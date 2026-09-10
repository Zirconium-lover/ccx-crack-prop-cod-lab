/*     CalculiX - damage/fracture extension                              */
/*     convstate.c: what counts as CONVERGED, made measurable first.     */

/* Why this module exists
   ----------------------
   06-TARGET.md lists four judgements with no owner.  This is one of them,
   and it is the one the second s3rad wall keeps pointing at.

   The history is short and worth stating.  At that wall the peak force
   residual sits on a node that cannot be equilibrated by any displacement a
   specimen of this size could undergo: node 1246 is a fragment held by ONE
   live bulk element with all six of its cohesive facets failed and open, so
   equilibrating it needs a motion of order the grip travel.  The response so
   far has been CCX_DAMAGE_AUTOSPC_FORCE, which removes such nodes from the
   force norm using a STIFFNESS mask.  Two measurements say that is the wrong
   handle:

     - it has become a fiction by 02-DIAGNOSTICS.md 4's own criterion: at the
       stall it excludes 4.4835e-02 against a tolerance of 1.8355e-03, 24.4x,
       larger than the residual it retains, and not returning to zero;
     - it is not robust: the recorded 554 -> 930 gain holds on two threads
       and not on four, where the run stalls at 604.

   A stiffness ratio also needs a tuned threshold, and 02-DIAGNOSTICS.md 2
   already records why that is a dead end - 76 nodes below 1e-3, 166 below
   1e-2, 381 below 1e-1, and picking among them proves nothing.

   The quantity that needs no threshold is in the same section:

       need_du = |R| / k     against a physical length

   It is a length.  Compared with the grip travel it is dimensionless, and it
   states the thing directly: can this node be equilibrated by a displacement
   the specimen could actually undergo?

   THIS MODULE DECIDES NOTHING YET, ON PURPOSE.  There is no deck in the tree
   that reproduces the class - 03-TEST-LADDER.md calls that the missing rung
   and the highest-value test work available - so changing the convergence
   criterion now would be changing a judgement nobody can put a gate case
   against.  What is missing first is the measurement, so this module
   computes and reports it and takes no decision at all.  Consumers ask; for
   the moment nobody acts.

   That order is the one damstate.c and loadpath.c were built in, and both
   times the measurement changed what the fix should be.                  */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CalculiX.h"

/* The displacement this node would need to shed its residual, given the
   stiffness actually holding it.  A non-positive diagonal has no meaning
   here and returns 0, which the report prints as "unknown" rather than as
   "needs nothing". */
double convstate_needdu(double r,double k)
{
  if(k<=0.) return 0.;
  return fabs(r)/k;
}

/* The physical length to judge it against: how far the driven grip has
   actually moved.  Taken from the driven node set loadpath.c already owns,
   so the two modules cannot disagree about what a grip is. */
double convstate_griplength(const double *vold,ITG mt,
                            const ITG *nodes,ITG n,ITG idir)
{
  ITG i,node; double u,m=0.;
  if((vold==NULL)||(nodes==NULL)||(n<=0)) return 0.;
  if((idir<1)||(idir>3)) idir=1;
  for(i=0;i<n;i++){
    node=nodes[i];
    if(node<1) continue;
    u=fabs(vold[mt*(node-1)+idir]);
    if(u>m) m=u;
  }
  return m;
}

/* One line, in the form 02-DIAGNOSTICS.md 2 asks it to be read in. */
void convstate_report(ITG inc,ITG node,double r,double k,double grip,
                      double tol)
{
  double need=convstate_needdu(r,k);
  printf("[CONVSTATE] inc=%" ITGFORMAT " peak residual node %" ITGFORMAT
         " |R|=%.6e k=%.6e tol=%.6e",inc,node,fabs(r),k,tol);
  if(k<=0.){
    printf(" need_du=unknown (non-positive diagonal)\n");
  }else if(grip>0.){
    printf(" need_du=%.6e grip=%.6e need_du/grip=%.4e\n",need,grip,need/grip);
  }else{
    printf(" need_du=%.6e grip=unknown\n",need);
  }
  fflush(stdout);
}

/* ---------------------------------------------------------------- tests */

static ITG cs_chk(const char *name,double got,double want,double tol,
                  ITG *nbad)
{
  ITG ok=(fabs(got-want)<=tol*(1.+fabs(want)));
  printf("   %-40s got=%-12.5g want=%-12.5g %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

ITG convstate_selftest(void)
{
  ITG nbad=0;

  printf("[CONVSTATE] self test%s","\n");

  /* A: need_du is a length, |R|/k */
  cs_chk("A need_du = |R|/k",convstate_needdu(2.,4.),0.5,1e-12,&nbad);
  cs_chk("A sign of R does not matter",convstate_needdu(-2.,4.),0.5,1e-12,
         &nbad);

  /* B: a non-positive diagonal is not "needs nothing".  The historical
     defect this guards is reporting 0 and having a consumer read it as
     healthy. */
  cs_chk("B non-positive diagonal -> 0, not huge",
         convstate_needdu(2.,0.),0.,1e-12,&nbad);
  cs_chk("B negative diagonal likewise",convstate_needdu(2.,-4.),0.,1e-12,
         &nbad);

  /* C: the grip length is the largest driven displacement, in the driven
     direction only. */
  {
    ITG mt=4,nodes[3]; double vold[16];
    ITG i; for(i=0;i<16;i++) vold[i]=0.;
    nodes[0]=1; nodes[1]=2; nodes[2]=3;
    vold[mt*0+1]=0.10; vold[mt*1+1]=0.25; vold[mt*2+1]=0.05;
    vold[mt*0+2]=9.99;                       /* another direction: ignored */
    cs_chk("C grip length is the largest driven |u|",
           convstate_griplength(vold,mt,nodes,3,1),0.25,1e-12,&nbad);
    cs_chk("C a different direction is read separately",
           convstate_griplength(vold,mt,nodes,3,2),9.99,1e-12,&nbad);
    cs_chk("C no nodes -> unknown, not zero-by-accident",
           convstate_griplength(vold,mt,NULL,0,1),0.,1e-12,&nbad);
  }

  /* D: the reading 02-DIAGNOSTICS.md 2 is built on.  A healthy node needs a
     microscopic motion; the second-wall node needs one of order the grip. */
  {
    double healthy=convstate_needdu(1.e-3,1.e3);      /* 1e-6 */
    double stuck  =convstate_needdu(4.4835e-2,1.e-1); /* 0.448 */
    double grip=0.2550;
    cs_chk("D healthy node: need_du/grip is microscopic",
           healthy/grip,3.92e-6,1e-2,&nbad);
    cs_chk("D stranded node: need_du/grip is of order one",
           stuck/grip,1.758,1e-2,&nbad);
  }

  printf("[CONVSTATE] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
