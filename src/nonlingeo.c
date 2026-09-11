/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2025 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include "CalculiX.h"
#include "mortar.h"
#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif
#ifdef PASTIX
#include "pastix.h"
#endif

#define max(a,b) ((a) >= (b) ? (a) : (b))

/* Damage Patch A3 adaptive-fast controller.
   The values below only affect trial load placement.  No trial damage is
   committed unless the post-deletion Newton/active-set solve converges. */
#define DAMAGE_FAST_BATCH_MAX 24
#define DAMAGE_FAST_MIN_FRACTION 0.25
#define DAMAGE_FAST_DIRECT_ALPHA 0.35

/* DE1.1 fixed-point controller.
   DE1 damage is coupled to equilibrium by same-load staggered iterations.
   Requiring bitwise/no-change closure is unnecessarily strict and can make
   one physical increment repeat indefinitely.  Use an engineering tolerance
   on the largest element damage correction, allow a slightly relaxed closure
   at the pass limit, and otherwise cut the physical increment back. */
#define DAMAGE_DE1_FP_TOL 1.e-4
#define DAMAGE_DE1_FP_RELAX_TOL 5.e-4
#define DAMAGE_DE1_MAX_PASSES 8
#define DAMAGE_DE1_CUTBACK_FACTOR 0.50

/* DE1.3 terminal-failure controller.
   Continuous DE1.2 degradation remains fully Newton-integrated.  Topology
   changes are allowed only very near complete material failure.  Limit one
   topology update to a modest batch; additional terminal elements are picked
   up by the existing same-load active-set closure after redistribution. */
/* bounds on how fast the dissipation-based step control may change
   the increment: shrinking hard is safe, growing back has to be
   gentle or the controller oscillates across the event */
#define DAMAGE_DISS_GROW 1.25
#define DAMAGE_DISS_SHRINK 0.20
/* largest lambda move one iteration may make, in units of the
   nominal increment: a bad secant must not throw the load factor
   across the whole step */
#define DAMAGE_DISS_DLAM 0.50
/* fraction of the target dissipation an accepted increment has to
   reach before the load factor is handed to the constraint */
#define DAMAGE_DISS_ENGAGE 0.10

#define DAMAGE_DE13_DELETE_D 0.999
#define DAMAGE_DE13_BATCH_MAX 64

/* a node whose entire live support is below this fraction of its original
   stiffness is treated as free for stabilisation purposes (E-61) */
#define DAMAGE_STAB_GDEAD 1.e-2

/* Damage-aware slow-Newton controller.
   The stock convergence controller starts estimating convergence at ir and
   cuts the increment at ic (normally 16), even when the residual and the
   displacement correction are both contracting monotonically.  DE1.2 uses
   an approximate damaged tangent, so linear convergence is expected while D
   evolves.  Permit a bounded extension only for that state.  The ordinary
   checkconvergence() residual/correction criteria and all of its divergence
   paths remain unchanged. */
#define DAMAGE_SLOW_NEWTON_MAX_ITERS 40
#define DAMAGE_SLOW_NEWTON_MAX_EXTRA 20
#define DAMAGE_SLOW_NEWTON_D_TOL 1.e-12
#define DAMAGE_SLOW_NEWTON_INVALID_EST 1000000000

/* Damage-specific Newton globalization.
   A full Newton correction is retained whenever it contracts the force
   residual.  Only a genuine residual increase in an active DE1.2/DM2.0
   softening trial triggers one safeguarded secant line-search correction.
   This adds a constitutive/residual evaluation, but no PARDISO
   factorization.  Final convergence tolerances remain the stock ones. */
#define DAMAGE_LINESEARCH_GROWTH 1.10
#define DAMAGE_LINESEARCH_MIN 0.10
#define DAMAGE_LINESEARCH_MAX 0.80
#define DAMAGE_LINESEARCH_FALLBACK 0.50
#define DAMAGE_LINESEARCH_MAX_TRIALS 3

/* AUTOSPC.  Nodes whose assembled diagonal has collapsed relative to their
   own intact value are excluded from the DISPLACEMENT convergence norm.
   Read by resultsini.c, which is where the Newton correction is still in
   hand - by the time control returns to nonlingeo, b holds residual forces
   and cam[0] cannot be recomputed.
   NULL, or a zero length, means the feature is off.

   WHY THIS IS SAFE.  The excluded DOF is still solved, still updated, and
   still measured by the FORCE residual ram[0], which is untouched.  Only its
   veto over cam[0] is removed.  A node that still carries load therefore
   still blocks convergence through the force criterion; what is dropped is a
   correction whose size is set by conditioning rather than by physics.
   Nothing is deleted, which is what separates this from E-22, E-57 and E-72.

   WHY IT IS NEEDED.  checkconvergence.c gates on cam[0] <= c2*uam[0] with
   c2 = 0.01, and on m12_epsf50 the run died with a force residual 338x
   inside tolerance and a correction/increment ratio of 0.29, both reported
   at a node with no stiffness left (E-67).  62-68% of the blocking reports
   in the final phase name a topologically defective node. */
ITG *damage_spc_mask=NULL;
ITG damage_spc_nk=0;
ITG damage_spc_count=0;

/* [DAMSTATE] the single owner of the load-path judgement; the three
   globals above are views onto it once it has been updated. */
damstate damage_dstate={0,NULL,NULL,NULL,NULL,0,0.,0};
/* [STALLSTATE] the run's own opinion of whether it is still progressing */
static stallstate damage_stall;
static ITG damage_stall_ok=0;

/* [LOADPATH] the single owner of the judgement one level up: whether a
   load path still links the grips at all, i.e. whether what is being
   solved is still a specimen.  Asked once per converged increment. */
loadpath damage_lp;
ITG damage_lp_ready=0;
ITG damage_lp_stop=1;
ITG damage_lp_armed_once=0;
ITG damage_lp_testok=0;

/* Number of active damage integration points for the standard 3-D
   continuum elements used by calcdamage.  For uncommon/composite
   formulations fall back to mi[0], i.e. the allocated damage stride. */
static void damage_aba_cmp(const char *name,const double *a,
                           const double *b,ITG n,ITG *nbad)
{
  ITG k,first=-1,ndiff=0;
  double mx=0.,d;
  if((a==NULL)||(b==NULL)||(n<=0)){
    printf("   %-14s SKIPPED (null or empty)%s",name,"\n");
    return;
  }
  for(k=0;k<n;k++){
    if(a[k]!=b[k]){
      ndiff++;
      if(first<0) first=k;
      d=fabs(a[k]-b[k]);
      if(d>mx) mx=d;
    }
  }
  if(ndiff==0){
    printf("   %-14s BITWISE IDENTICAL  (n=%" ITGFORMAT ")%s",
           name,n,"\n");
  }else{
    (*nbad)++;
    printf("   %-14s *** DIFFERS *** %" ITGFORMAT " of %" ITGFORMAT
           " entries, first k=%" ITGFORMAT " (%.17e vs %.17e) maxabs=%.6e%s",
           name,ndiff,n,first,a[first],b[first],mx,"\n");
  }
}

static ITG damage_history_nip(const char *lakonel,ITG mi0);

/* Project f_hat (equation space) onto a nodal displacement difference.
   Only free dofs contribute, and j starts at 1 because j=0 is the thermal
   dof, exactly as resultsini.c applies the correction.  This is how the
   path following reads the model instead of shadowing it. */

static double damage_fd_dot(const double *a,const double *b,ITG n){
  ITG i;
  double s=0.;
  for(i=0;i<n;i++) s+=a[i]*b[i];
  return s;
}

static double pf_project(const double *fh,const double *a,const double *c,
                         const ITG *nactdof,ITG nk,ITG mt){
  ITG i,j,k;
  double p=0.;
  for(i=0;i<nk;i++){
    for(j=1;j<mt;j++){
      k=nactdof[mt*i+j];
      if(k>0) p+=fh[k-1]*(a[mt*i+j]-c[mt*i+j]);
    }
  }
  return p;
}

/* ---- discrete-branch census -------------------------------------------
   J-14 measured a sharp loss of local linearity between alpha=0.0625 and
   alpha=0.125 and called it a "discrete switch".  That was a hypothesis, not
   a measurement: nothing had been counted.  These two helpers count it.

   One bitmask per integration point.  The element type decides which bits are
   meaningful, and getting that wrong is not hypothetical - xstate slot 1 is
   the equivalent plastic strain for a bulk C3D4 (calcdamage.f:649) and dmax
   for a UC6 facet (resultsmech_uc6.f:43).  Reading one as the other yields a
   plausible, meaningless census.

   Layout: dam(mi(1),*)          -> dam[mi0*i+j]
           xstate(nstate_,mi(1),*) -> xstate[nstate_*mi0*i + nstate_*j + k] */
ITG ccx_rescue_active=0,ccx_rescue_arm=0,ccx_rescue_req=0;

#define DAMCAT_PLAST   1   /* bulk: accumulated plastic strain this increment */
#define DAMCAT_DINIT   2   /* bulk: damage initiated (dam >= 1)               */
#define DAMCAT_DGROW   4   /* bulk: damage grew from the committed baseline   */
#define DAMCAT_USOFT   8   /* UC6 : dmax > 0, law has left the elastic branch */
#define DAMCAT_UVISC  16   /* UC6 : viscous damage active                     */
#define DAMCAT_UFAIL  32   /* UC6 : fully failed flag set                     */
#define DAMCAT_UADV   64   /* UC6 : dmax ADVANCED this trial, i.e. deff>dmax0  */
#define DAMCAT_UCOMP 128   /* UC6 : traction(1)<0, the compression branch     */

/* ==================================================================
   [DAMAGE CT] bounded experimental coupled local continuation (SPEC
   FREEZE v1).  Opt-in, arms only as rescue LEVEL 4 after Rescue2 and the
   dogleg have both failed on one wall.  Everything below is inert unless
   CCX_DAMAGE_CONTINUATION is set.
   ================================================================== */

/* Exact UC6 kinematics, mirroring cohesive_uc6.f:96-160.  rmat rows map a
   global vector to (normal, shear-1, shear-2) and are built from the
   REFERENCE geometry co alone, so they are a mesh constant.  shape is the
   integration-point specific weight set: 1/6 off-point, 2/3 on-point. */

static void damage_ct_kin(const double *co,const ITG *kon,ITG indexe,
                          const double *v,ITG mt,ITG mint,
                          double *dl,double *rmat,double *shape)
{
  double e1[3],e2[3],cv[3],jump[3],n1,nc;
  ITG i,k,nm,np;

  for(k=0;k<3;k++){
    e1[k]=co[3*(kon[indexe+1]-1)+k]-co[3*(kon[indexe]-1)+k];
    e2[k]=co[3*(kon[indexe+2]-1)+k]-co[3*(kon[indexe]-1)+k];
  }
  cv[0]=e1[1]*e2[2]-e1[2]*e2[1];
  cv[1]=e1[2]*e2[0]-e1[0]*e2[2];
  cv[2]=e1[0]*e2[1]-e1[1]*e2[0];
  n1=sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
  nc=sqrt(cv[0]*cv[0]+cv[1]*cv[1]+cv[2]*cv[2]);
  if((n1<=1.e-30)||(nc<=1.e-30)){
    for(i=0;i<3;i++){dl[i]=0.;for(k=0;k<3;k++) rmat[3*i+k]=0.;}
    return;
  }
  for(k=0;k<3;k++){rmat[k]=cv[k]/nc;rmat[3+k]=e1[k]/n1;}
  rmat[6]=rmat[1]*rmat[5]-rmat[2]*rmat[4];
  rmat[7]=rmat[2]*rmat[3]-rmat[0]*rmat[5];
  rmat[8]=rmat[0]*rmat[4]-rmat[1]*rmat[3];

  for(i=0;i<3;i++) shape[i]=1./6.;
  shape[mint]=2./3.;
  for(k=0;k<3;k++){
    jump[k]=0.;
    for(i=0;i<3;i++){
      nm=kon[indexe+i]-1;
      np=kon[indexe+i+3]-1;
      jump[k]+=shape[i]*(v[mt*np+k+1]-v[mt*nm+k+1]);
    }
  }
  for(i=0;i<3;i++){
    dl[i]=0.;
    for(k=0;k<3;k++) dl[i]+=rmat[3*i+k]*jump[k];
  }
}

/* One committed endpoint of the ring: the three local separations of every
   live UC6 integration point, plus its failed flag and its compression
   category.  Both flags are needed because the historical admissibility
   guards test them at BOTH endpoints of an interval. */

static void damage_ct_snap(const double *co,const ITG *kon,const ITG *ipkon,
                           const char *lakon,const double *v,
                           const double *stx,const double *xstate,
                           ITG ne0,ITG mi0,ITG nstate,ITG mt,
                           double *ring,ITG *fl)
{
  double dl[3],rmat[9],shape[3];
  ITG i,j,np,idx;

  np=(mi0<3)?mi0:3;
  for(i=0;i<ne0;i++){
    for(j=0;j<mi0;j++){
      idx=mi0*i+j;
      ring[3*idx]=0.;ring[3*idx+1]=0.;ring[3*idx+2]=0.;
      fl[idx]=-1;                                  /* -1: not a live UC6 ip */
    }
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    for(j=0;j<np;j++){
      idx=mi0*i+j;
      damage_ct_kin(co,kon,ipkon[i],v,mt,j,dl,rmat,shape);
      ring[3*idx]=dl[0];ring[3*idx+1]=dl[1];ring[3*idx+2]=dl[2];
      fl[idx]=((xstate[nstate*idx+3]>=0.5)?1:0)
             |((stx[6*idx]<0.)?2:0);              /* bit0 failed, bit1 comp */
    }
  }
}

/* Predicted reduction is not used here; the continuation needs only the
   bordered algebra.  Isolated so the self-test exercises the same code.

   den = c_lambda + c_u^T y ,  dlambda = (-c - c_u^T z)/den
   Returns 0 on a degenerate or non-finite input - a REFUSAL, never a step
   of zero length.  Every test is written positively so NaN falls through
   to the refusal branch. */

static ITG damage_ct_bordered(double clam,double cuz,double cuy,double c,
                              double *den,double *dlam)
{
  *den=clam+cuy;
  if(!(*den==*den)) return 0;
  if(!(fabs(*den)>0.)) return 0;
  *dlam=(-c-cuz)/(*den);
  if(!(*dlam==*dlam)) return 0;
  return 1;
}

/* Conditioning of the constraint row, measured on the constraint SUPPORT
   only: the raw |den| is not comparable between candidates. */

static double damage_ct_rhoden(double clam,double cunorm,double ysupp,
                               double den)
{
  double d=fabs(clam)+cunorm*ysupp;
  if(!(d>0.)) return 0.;
  return fabs(den)/d;
}

/* [DAMAGE CT] candidate scan.  Uses ONLY committed ring data, so it can be
   evaluated BEFORE the level-4 retry is granted: a refusal then costs the
   stock path nothing at all and the wall goes to the ORIGINAL stock stop
   byte-identically.  Returns 1 and fills the frozen quantities on success,
   0 on refusal.  reason: 1 no candidate, 2 kappa unstable. */

static ITG damage_ct_select(const double *ring,const ITG *fl,
                            const double *dtr,const ITG *ipkon,
                            const char *lakon,const ITG *ielprop,
                            const double *prop,ITG ne0,ITG mi0,ITG head,
                            ITG *belem,ITG *bip,double *m,double *kappa,
                            double *ds0,double *tau,ITG *ncand,ITG *reason,
                            double kaptol,const ITG *bl,ITG nbl)
{
  ITG i,j,k,n,s0,s1,f0,f1,ok,cnt=0;
  double tn0,ts0,gc,beta,df,atau,d0[3],d1[3],mm[3],deff,adv,best=-1.;
  double q[5],kp[5],sq[5],sk[5],t;

  n=mi0*ne0;*ncand=0;*reason=1;*belem=-1;*bip=-1;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    if(ielprop[i]<0) continue;
    tn0=prop[ielprop[i]+1];ts0=prop[ielprop[i]+2];gc=prop[ielprop[i]+3];
    if((tn0<=0.)||(ts0<=0.)||(gc<=0.)) continue;
    beta=(ts0/tn0)*(ts0/tn0);
    df=2.*gc/tn0;
    atau=1.e-9*df;
    for(j=0;j<3;j++){
      if(j>=mi0) break;
      f1=fl[n*head+mi0*i+j];
      if(f1<0) continue;
      if(f1&1) continue;
      {ITG bk,skip=0;
       for(bk=0;bk<nbl;bk++) if(bl[bk]==4*i+j){skip=1;break;}
       if(skip) continue;}
      for(k=0;k<3;k++) d1[k]=ring[3*(n*head+mi0*i+j)+k];
      deff=(d1[0]>0.?d1[0]*d1[0]:0.)+beta*(d1[1]*d1[1]+d1[2]*d1[2]);
      deff=(deff>0.)?sqrt(deff):0.;
      if(!(deff>0.)) continue;
      if(fabs(d1[0])<0.05*deff) continue;
      mm[0]=(d1[0]>0.?d1[0]:0.)/deff;
      mm[1]=beta*d1[1]/deff;
      mm[2]=beta*d1[2]/deff;
      ok=1;
      for(k=0;k<5;k++){
        s1=(head-k+6)%6;s0=(head-k-1+6)%6;
        f1=fl[n*s1+mi0*i+j];f0=fl[n*s0+mi0*i+j];
        if((f1<0)||(f0<0)){ok=0;break;}
        if((f1&1)||(f0&1)){ok=0;break;}
        if((f1&2)!=(f0&2)){ok=0;break;}
        if((ring[3*(n*s1+mi0*i+j)]>0.)!=(ring[3*(n*s0+mi0*i+j)]>0.)){
          ok=0;break;}
        adv=mm[0]*(ring[3*(n*s1+mi0*i+j)]-ring[3*(n*s0+mi0*i+j)])
           +mm[1]*(ring[3*(n*s1+mi0*i+j)+1]-ring[3*(n*s0+mi0*i+j)+1])
           +mm[2]*(ring[3*(n*s1+mi0*i+j)+2]-ring[3*(n*s0+mi0*i+j)+2]);
        if(!(adv>atau)){ok=0;break;}
        if(!(dtr[s1]>0.)){ok=0;break;}
        q[k]=adv;kp[k]=dtr[s1]/adv;
        if(!(kp[k]>0.)){ok=0;break;}
      }
      if(ok==0) continue;
      cnt++;
      if(q[0]>best){
        ITG b1,b2;
        best=q[0];*belem=i;*bip=j;
        for(k=0;k<3;k++) m[k]=mm[k];
        for(k=0;k<5;k++){sq[k]=q[k];sk[k]=kp[k];}
        for(b1=1;b1<5;b1++){
          t=sq[b1];for(b2=b1;(b2>0)&&(sq[b2-1]>t);b2--) sq[b2]=sq[b2-1];
          sq[b2]=t;
          t=sk[b1];for(b2=b1;(b2>0)&&(sk[b2-1]>t);b2--) sk[b2]=sk[b2-1];
          sk[b2]=t;
        }
        *ds0=sq[2];*kappa=sk[2];
        *tau=(sk[0]>0.)?sk[4]/sk[0]:1.e30;
      }
    }
  }
  *ncand=cnt;
  if((cnt==0)||(*belem<0)){*reason=1;return 0;}
  if(!(*tau<=kaptol)){*reason=2;return 0;}
  *reason=0;
  return 1;
}

/* [DAMAGE CT SELFTEST] bordered algebra on a case with a closed-form
   answer.  Returns the number of failures; a non-zero result stops the run
   before the first increment. */

static ITG damage_ct_selftest(void)
{
  double den,dlam,r;
  ITG k,nbad=0;

  /* c_lambda=0, c_u^T y = 4, c_u^T z = 1, c = 2  ->  den=4, dlam=-3/4 */
  k=damage_ct_bordered(0.,1.,4.,2.,&den,&dlam);
  printf("[DAMAGE CT SELFTEST] den=%.12e (want 4) dlam=%.12e (want -0.75) "
         "rc=%" ITGFORMAT " %s\n",den,dlam,k,
         ((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam+0.75)<1.e-14))?
         "PASS":"FAIL");
  if(!((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam+0.75)<1.e-14))) nbad++;

  /* c_lambda non-zero must enter den */
  k=damage_ct_bordered(3.,0.,1.,-8.,&den,&dlam);
  printf("[DAMAGE CT SELFTEST] den=%.12e (want 4) dlam=%.12e (want 2) "
         "rc=%" ITGFORMAT " %s\n",den,dlam,k,
         ((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam-2.)<1.e-14))?
         "PASS":"FAIL");
  if(!((k==1)&&(fabs(den-4.)<1.e-14)&&(fabs(dlam-2.)<1.e-14))) nbad++;

  /* c=0 at a converged constraint must give dlam = -c_u^T z / den */
  k=damage_ct_bordered(0.,2.,2.,0.,&den,&dlam);
  printf("[DAMAGE CT SELFTEST] c=0 -> dlam=%.12e (want -1) %s\n",dlam,
         ((k==1)&&(fabs(dlam+1.)<1.e-14))?"PASS":"FAIL");
  if(!((k==1)&&(fabs(dlam+1.)<1.e-14))) nbad++;

  /* degenerate and non-finite input must REFUSE */
  {
    double z=0.,dn=z/z;
    ITG k1,k2,k3;
    k1=damage_ct_bordered(0.,1.,0.,1.,&den,&dlam);      /* den=0     */
    k2=damage_ct_bordered(dn,1.,1.,1.,&den,&dlam);      /* clam=NaN  */
    k3=damage_ct_bordered(0.,dn,1.,1.,&den,&dlam);      /* cuz=NaN   */
    printf("[DAMAGE CT SELFTEST] refusals: den=0 -> %" ITGFORMAT
           ", clam=NaN -> %" ITGFORMAT ", cuz=NaN -> %" ITGFORMAT
           " (all want 0) %s\n",k1,k2,k3,
           ((k1==0)&&(k2==0)&&(k3==0))?"PASS":"FAIL");
    if(!((k1==0)&&(k2==0)&&(k3==0))) nbad++;
  }

  /* rho_den is a cancellation measure in [0,1] and 0 on a null row */
  r=damage_ct_rhoden(0.,1.,4.,4.);
  printf("[DAMAGE CT SELFTEST] rho_den=%.12e (want 1) ; null row -> %.12e "
         "(want 0) %s\n",r,damage_ct_rhoden(0.,0.,0.,0.),
         ((fabs(r-1.)<1.e-14)&&(damage_ct_rhoden(0.,0.,0.,0.)==0.))?
         "PASS":"FAIL");
  if(!((fabs(r-1.)<1.e-14)&&(damage_ct_rhoden(0.,0.,0.,0.)==0.))) nbad++;

  printf("[DAMAGE CT SELFTEST] %" ITGFORMAT " failure(s)\n",nbad);
  fflush(stdout);
  return nbad;
}

/* ---- [DAMAGE TR] the dogleg step, isolated ------------------------------

   Given the trust radius dl and the five scalars the trust region works in
   (nd2=|d|^2, nw2=|Jd|^2, npn2=|p_N|^2, dtpn=dot(d,p_N), and tc=nd2/nw2 for
   the Cauchy point), return the step as p = (*pa)*d + (*pb)*p_N and its norm.

   Return value: 1 NEWTON (the full step is inside the radius), 2 CAUCHY (the
   Cauchy point is already outside it), 3 DOGLEG (the blend that meets the
   boundary), 0 REFUSE - a degenerate or non-finite input, which the caller
   must treat as "no step", never as a step of zero length.

   Every test is written as a POSITIVE comparison (!(x>0.) rather than x<=0.)
   so that a NaN falls into the refusal branch instead of silently passing. */

static ITG damage_dl_pick(double dl,double nd2,double nw2,double npn2,
                          double dtpn,double *pa,double *pb,double *nrm)
{
  double tc,pcn,pnn,uu,uv,vv,disc,tau,q;

  *pa=0.;*pb=0.;*nrm=0.;
  if(!(dl>0.)) return 0;
  if(!(nw2>0.)) return 0;
  if(!(nd2>0.)) return 0;
  if(!(npn2>0.)) return 0;
  if(!(dtpn==dtpn)) return 0;

  tc=nd2/nw2;
  if(!(tc>0.)) return 0;
  pnn=sqrt(npn2);
  pcn=tc*sqrt(nd2);
  if(!(pnn>0.)||!(pcn>0.)) return 0;

  if(pnn<=dl){
    *pa=0.;*pb=1.;
  }else if(pcn>=dl){
    *pa=dl/sqrt(nd2);*pb=0.;
  }else{
    uu=pcn*pcn;
    uv=tc*dtpn-uu;
    vv=npn2-2.*tc*dtpn+uu;
    if(!(vv>0.)){
      *pa=0.;*pb=1.;
      q=npn2;
      *nrm=sqrt(q);
      return 1;
    }
    disc=uv*uv-vv*(uu-dl*dl);
    if(!(disc>=0.)) disc=0.;
    tau=(-uv+sqrt(disc))/vv;
    if(!(tau==tau)) return 0;
    if(tau<0.) tau=0.;
    if(tau>1.) tau=1.;
    *pa=(1.-tau)*tc;*pb=tau;
    q=(*pa)*(*pa)*nd2+2.*(*pa)*(*pb)*dtpn+(*pb)*(*pb)*npn2;
    *nrm=(q>0.)?sqrt(q):0.;
    if(!((*nrm)>0.)) return 0;
    return 3;
  }
  q=(*pa)*(*pa)*nd2+2.*(*pa)*(*pb)*dtpn+(*pb)*(*pb)*npn2;
  *nrm=(q>0.)?sqrt(q):0.;
  if(!((*nrm)>0.)) return 0;
  return (*pb>0.)?1:2;
}

/* predicted reduction of phi=1/2|R|^2 for p = pa*d + pb*p_N, using
   J p = pa*w + pb*r0 and R = -r0:  phi(u) - 1/2|R+Jp|^2 */

static double damage_dl_pred(double pa,double pb,double nb2,double nd2,
                             double nw2)
{
  return 0.5*nb2-0.5*((pb-1.)*(pb-1.)*nb2+2.*pa*(pb-1.)*nd2+pa*pa*nw2);
}

/* ---- [DAMAGE TR] geometry self-test ------------------------------------

   Runs on a J whose answer is known in closed form, so a wrong branch is a
   failure of arithmetic and not of the model:  J=diag(1,10), r0=(1,1).
       p_N = (1, 0.1)      |p_N|^2 = 1.01
       d   = J^T r0 = (1,10)   |d|^2   = 101
       w   = J d    = (1,100)  |w|^2   = 10001
       |r0|^2 = 2  and  dot(d,p_N) = 2   - the same identity the run checks.
   Prints PASS/FAIL per case and returns the number of failures. */

static ITG damage_dl_selftest(void)
{
  const double nb2=2.,nd2=101.,nw2=10001.,npn2=1.01,dtpn=2.;
  double pa,pb,nr,pr,pnn,pcn,tc;
  ITG k,nbad=0;

  tc=nd2/nw2;pnn=sqrt(npn2);pcn=tc*sqrt(nd2);
  printf("[DAMAGE TR SELFTEST] closed-form case: |p_N|=%.12e |p_C|=%.12e\n",
         pnn,pcn);

  /* 1. the full Newton step lies inside the radius */
  k=damage_dl_pick(2.,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
  pr=damage_dl_pred(pa,pb,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] Delta=2 -> kind=%" ITGFORMAT
         " (want 1 NEWTON) pa=%.6e pb=%.6e |p|=%.12e (want %.12e) "
         "pred=%.12e %s\n",k,pa,pb,nr,pnn,pr,
         ((k==1)&&(fabs(nr-pnn)<=1.e-12*pnn)&&(pr>0.))?"PASS":"FAIL");
  if(!((k==1)&&(fabs(nr-pnn)<=1.e-12*pnn)&&(pr>0.))) nbad++;

  /* 2. the Cauchy point is already outside the radius */
  k=damage_dl_pick(0.05,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
  pr=damage_dl_pred(pa,pb,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] Delta=0.05 -> kind=%" ITGFORMAT
         " (want 2 CAUCHY) pa=%.6e pb=%.6e |p|=%.12e (want 5.0e-02) "
         "pred=%.12e %s\n",k,pa,pb,nr,pr,
         ((k==2)&&(pb==0.)&&(fabs(nr-0.05)<=1.e-12*0.05)&&(pr>0.))?
         "PASS":"FAIL");
  if(!((k==2)&&(pb==0.)&&(fabs(nr-0.05)<=1.e-12*0.05)&&(pr>0.))) nbad++;

  /* 3. the intermediate step must sit exactly ON the boundary */
  k=damage_dl_pick(0.5,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
  pr=damage_dl_pred(pa,pb,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] Delta=0.5 -> kind=%" ITGFORMAT
         " (want 3 DOGLEG) pa=%.6e pb=%.6e |p|=%.12e (want 5.0e-01) "
         "pred=%.12e %s\n",k,pa,pb,nr,pr,
         ((k==3)&&(fabs(nr-0.5)<=1.e-10*0.5)&&(pr>0.)&&(pa>0.)&&(pb>0.))?
         "PASS":"FAIL");
  if(!((k==3)&&(fabs(nr-0.5)<=1.e-10*0.5)&&(pr>0.)&&(pa>0.)&&(pb>0.))) nbad++;

  /* 4. the model prediction is exact for the full Newton step */
  pr=damage_dl_pred(0.,1.,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] pred(full Newton)=%.12e (want %.12e, i.e. the "
         "linear model predicts phi=0) %s\n",pr,0.5*nb2,
         (fabs(pr-0.5*nb2)<=1.e-14*nb2)?"PASS":"FAIL");
  if(!(fabs(pr-0.5*nb2)<=1.e-14*nb2)) nbad++;

  /* 5. the Cauchy point is the exact minimiser along d */
  pr=damage_dl_pred(tc,0.,nb2,nd2,nw2);
  printf("[DAMAGE TR SELFTEST] pred(Cauchy point)=%.12e (want %.12e = "
         "|d|^4/(2|Jd|^2)) %s\n",pr,0.5*nd2*nd2/nw2,
         (fabs(pr-0.5*nd2*nd2/nw2)<=1.e-12*fabs(pr))?"PASS":"FAIL");
  if(!(fabs(pr-0.5*nd2*nd2/nw2)<=1.e-12*fabs(pr))) nbad++;

  /* 6-9. degenerate and non-finite input must REFUSE, not return a step */
  {
    double dnan=0.,dzero=0.;
    ITG k6,k7,k8,k9;
    dnan=dzero/dzero;                       /* NaN without a literal */
    k6=damage_dl_pick(0.,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
    k7=damage_dl_pick(0.5,nd2,0.,npn2,dtpn,&pa,&pb,&nr);
    k8=damage_dl_pick(dnan,nd2,nw2,npn2,dtpn,&pa,&pb,&nr);
    k9=damage_dl_pick(0.5,nd2,nw2,npn2,dnan,&pa,&pb,&nr);
    printf("[DAMAGE TR SELFTEST] refusals: Delta=0 -> %" ITGFORMAT
           ", |Jd|^2=0 -> %" ITGFORMAT ", Delta=NaN -> %" ITGFORMAT
           ", dot(d,p_N)=NaN -> %" ITGFORMAT " (all want 0) %s\n",
           k6,k7,k8,k9,
           ((k6==0)&&(k7==0)&&(k8==0)&&(k9==0))?"PASS":"FAIL");
    if(!((k6==0)&&(k7==0)&&(k8==0)&&(k9==0))) nbad++;
  }

  printf("[DAMAGE TR SELFTEST] %" ITGFORMAT " failure(s)\n",nbad);
  fflush(stdout);
  return nbad;
}

static ITG damage_ray_catof(const double *xstate,const double *xstateini,
                            const double *dam,const double *dambase,
                            const double *visc,const double *stx,
                            const char *lakonel,
                            ITG i,ITG j,ITG mi0,ITG nstate)
{
  ITG c=0,ix,is;
  double d,db;
  ix=nstate*mi0*i+nstate*j;
  is=6*mi0*i+6*j;
  if(lakonel[0]=='C'){
    if(nstate>0){
      if(xstate[ix]-xstateini[ix]>1.e-14) c|=DAMCAT_PLAST;
    }
    /* NO real bulk branch flag.  Exporting mattyp into a spare xstate slot
       was attempted and REVERTED: incplas_lin.f:217 owns slots 2..7 for the
       plastic strain tensor, so slot 2 is epl(1), and writing there destroys
       the return map ("no convergence in incplas" at the first increment).
       With *Depvar 4 there is no spare slot at all.  PLAST above stays as
       the observable: peeq grows if and only if the return map ran, so it is
       equivalent by construction, but it IS an inference, not the flag. */
    d=dam[mi0*i+j];
    if(d>=1.) c|=DAMCAT_DINIT;
    if(dambase!=NULL){
      db=dambase[mi0*i+j];
      if(d-db>1.e-14) c|=DAMCAT_DGROW;
    }
  }else if(lakonel[0]=='U'){
    if(nstate>0){ if(xstate[ix]>1.e-14)   c|=DAMCAT_USOFT; }
    if(nstate>1){ if(xstate[ix+1]>1.e-14) c|=DAMCAT_UVISC; }
    if(nstate>3){ if(xstate[ix+3]>0.5)    c|=DAMCAT_UFAIL; }
    /* LOADING/UNLOADING.  cohesive_uc6.f stores dmax=max(dmax0,deff) in slot
       1 and takes the softening tangent only while deff>=dmax0, so
       xstate>xstateini in that slot IS "this point advanced its maximum on
       this trial" - the one active-set indicator of the cohesive law that
       the census did not carry.  A point sitting exactly on the kink has
       xstate==xstateini and is counted as NOT advancing, which is what makes
       a branch crossing along a search direction visible as a transition. */
    if(nstate>0){ if(xstate[ix]>xstateini[ix]) c|=DAMCAT_UADV; }
    /* tension/compression: resultsmech_uc6.f:55 stores traction(1), and both
       branches multiply deltal(1) by a positive factor (kn or g*kn), so the
       sign of the stored traction IS the sign of the normal opening.  This
       switch changes the facet stiffness by 1/g ~ 1e4 and was entirely absent
       from the earlier census. */
    if(stx!=NULL){ if(stx[is]<0.) c|=DAMCAT_UCOMP; }
  }
  return c;
}

static void damage_ray_census(ITG *cat,const double *xstate,
                              const double *xstateini,const double *dam,
                              const double *dambase,const double *visc,
                              const double *stx,
                              const ITG *ipkon,const char *lakon,
                              ITG ne0,ITG mi0,ITG nstate)
{
  ITG i,j,nip;
  for(i=0;i<ne0;i++){
    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    for(j=0;j<mi0;j++) cat[mi0*i+j]=0;
    if(ipkon[i]<0) continue;
    for(j=0;j<nip;j++)
      cat[mi0*i+j]=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                                    &lakon[8*i],i,j,mi0,nstate);
  }
}

/* [DAMAGE RAY] the two branch indicators, counted over the live mesh at
   whatever state the caller has just built.  Same iteration as the census
   above, so the two always speak about the same integration points. */

static void damage_ray_tally(const double *xstate,const double *xstateini,
                             const double *dam,const double *dambase,
                             const double *visc,const double *stx,
                             const ITG *ipkon,const char *lakon,
                             ITG ne0,ITG mi0,ITG nstate,
                             ITG *nplast,ITG *nucomp)
{
  ITG i,j,nip,c;
  *nplast=0;*nucomp=0;
  for(i=0;i<ne0;i++){
    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    if(ipkon[i]<0) continue;
    for(j=0;j<nip;j++){
      c=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                         &lakon[8*i],i,j,mi0,nstate);
      if(c&DAMCAT_PLAST) (*nplast)++;
      if(c&DAMCAT_UCOMP) (*nucomp)++;
    }
  }
}

/* [DAMAGE EVT] The UC6 tension/compression set.  traction(1) is kn*deltal(1)
   in compression and g*kn*deltal(1) in tension, and both kn and g are
   positive, so sign(stx(1)) IS sign(deltal(1)) and the branch is readable
   without storing deltal anywhere.  UC6 carries exactly three integration
   points. */

static void damage_evt_sign(const double *stx,const ITG *ipkon,
                            const char *lakon,ITG ne0,ITG mi0,ITG *sgn)
{
  ITG i,j,np;
  np=(mi0<3)?mi0:3;
  for(i=0;i<ne0;i++){
    for(j=0;j<mi0;j++) sgn[mi0*i+j]=0;
    if(ipkon[i]<0) continue;
    /* UC6 is labelled 'U'; 'C' is the C3D4 bulk (damage_ray_catof splits on
       exactly this).  Getting the letter wrong built the map over the bulk
       and tracked nothing - measured, s3rad inc=231 reported "no ladder alpha
       changes the set" while the ray saw the crossing. */
    if(lakon[8*i]!='U') continue;
    for(j=0;j<np;j++)
      sgn[mi0*i+j]=(stx[6*mi0*i+6*j]<0.)?-1:1;
  }
}

static ITG damage_evt_flips(const double *stx,const ITG *ipkon,
                            const char *lakon,ITG ne0,ITG mi0,const ITG *sgn,
                            ITG *firste,ITG *firstip)
{
  ITG i,j,np,n=0,sg;
  np=(mi0<3)?mi0:3;
  *firste=0;*firstip=0;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    for(j=0;j<np;j++){
      if(sgn[mi0*i+j]==0) continue;
      sg=(stx[6*mi0*i+6*j]<0.)?-1:1;
      if(sg!=sgn[mi0*i+j]){
        if(n==0){ *firste=i+1; *firstip=j+1; }
        n++;
      }
    }
  }
  return n;
}

static ITG damage_ray_census_diff(const ITG *cat,const double *xstate,
                                  const double *xstateini,const double *dam,
                                  const double *dambase,const double *visc,
                                  const double *stx,
                                  const ITG *ipkon,const char *lakon,
                                  ITG ne0,ITG mi0,ITG nstate,
                                  ITG *firste,ITG *firstip,
                                  ITG *firsta,ITG *firstb)
{
  ITG i,j,nip,c,n=0;
  *firste=-1;*firstip=-1;*firsta=0;*firstb=0;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    for(j=0;j<nip;j++){
      c=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                         &lakon[8*i],i,j,mi0,nstate);
      if(c!=cat[mi0*i+j]){
        n++;
        if(*firste<0){
          *firste=i+1;*firstip=j+1;*firsta=cat[mi0*i+j];*firstb=c;
        }
      }
    }
  }
  return n;
}

/* [WALLDIAG] The same census difference, but resolved per category bit and
   split between the two element families, because "the active set moved" and
   "8000 cohesive points crossed their loading/unloading kink" are different
   findings and the aggregate count cannot tell them apart.

   nb[k] counts the integration points whose bit k differs from the reference
   census cat.  Bits are the DAMCAT_ masks; the caller names them. */

/* [WALLDIAG] WHERE the residual and the Newton correction live.

   A norm says how big they are and nothing about what they touch, and the
   two competing explanations of the second wall differ precisely in what
   they touch: a globalisation defect spreads the correction over the mesh,
   a branch-switching one concentrates it on the fracture front.  So the
   probe reports the largest entries by NODE, and for each such node the
   local state of the front - how much of its support is left, how soft its
   cohesive facets are and how many of them are carrying compression.

   nactdof is the post-SPC numbering (1-based, mt-strided) that b, ad and au
   share, so inverting it is what turns an equation index back into a node
   and a direction.  It is inverted here rather than assumed. */

static void damage_wall_where(const char *tag,const double *x,ITG neq1,
                              const ITG *nactdof,ITG mt,ITG nk,ITG ntop,
                              const ITG *ipkon,const ITG *kon,
                              const char *lakon,const double *xstate,
                              const double *stx,const double *dam,
                              ITG ne,ITG ne0,ITG mi0,ITG nstate)
{
  ITG *inode=NULL,*idir=NULL,i,j,k,t,ip,np,nb,nu,ncomp,ibest,idx;
  double gmn,gmx,dv,a;

  if((ntop<1)||(neq1<1)) return;
  NNEW(inode,ITG,neq1);
  NNEW(idir,ITG,neq1);
  for(i=0;i<neq1;i++){inode[i]=-1;idir[i]=0;}
  for(i=0;i<nk;i++){
    for(j=1;j<mt;j++){
      k=nactdof[mt*i+j];
      if((k>0)&&(k<=neq1)){inode[k-1]=i;idir[k-1]=j;}
    }
  }

  for(t=0;t<ntop;t++){
    ibest=-1;a=-1.;
    for(i=0;i<neq1;i++){
      if(inode[i]<0) continue;          /* already reported, or unmapped */
      if(fabs(x[i])>a){a=fabs(x[i]);ibest=i;}
    }
    if(ibest<0) break;
    i=inode[ibest];
    nb=0;nu=0;ncomp=0;gmn=2.;gmx=-1.;
    for(j=0;j<ne;j++){
      if(ipkon[j]<0) continue;
      if(lakon[8*j]=='U') np=6; else if(lakon[8*j]=='C') np=4; else continue;
      ip=0;
      for(k=0;k<np;k++) if(kon[ipkon[j]+k]-1==i) ip=1;
      if(ip==0) continue;
      if(lakon[8*j]=='C'){nb++;continue;}
      nu++;
      if(j>=ne0) continue;
      for(k=0;k<3;k++){
        idx=mi0*j+k;
        if(nstate>1){
          dv=1.-xstate[nstate*idx+1];
          if(dv<gmn) gmn=dv;
          if(dv>gmx) gmx=dv;
        }
        if(stx[6*idx]<0.) ncomp++;
      }
    }
    printf("[WALLDIAG]   %s #%" ITGFORMAT ": node %" ITGFORMAT " dir %"
           ITGFORMAT " value %.6e   live bulk %" ITGFORMAT " live UC6 %"
           ITGFORMAT " facet g in [%.3e,%.3e] UC6 ips in compression %"
           ITGFORMAT "%s",tag,t+1,i+1,idir[ibest],x[ibest],nb,nu,
           (gmn<=1.)?gmn:0.,(gmx>=0.)?gmx:0.,ncomp,"\n");
    inode[ibest]=-1;                       /* do not report it twice */
  }
  SFREE(inode);SFREE(idir);
}

/* [WALLDIAG] WHERE a vector lives, split by what its node touches.

   The linear-model defect at small eps IS eps*(J - dR/du)p, so localising it
   localises the missing term: a defect that sits on nodes with cohesive
   facets accuses the interface tangent, one that sits on nodes whose
   elements are softening accuses the damage rank-1 term, and one spread over
   ordinary plastic bulk accuses the return map.  Shares of |x|_2^2, so they
   add to 1. */

static void damage_wall_split(const char *tag,const double *x,ITG neq1,
                              const ITG *nactdof,ITG mt,ITG nk,
                              const ITG *ipkon,const ITG *kon,
                              const char *lakon,const double *dam,
                              const double *xstate,ITG ne,ITG ne0,ITG mi0,
                              ITG nstate)
{
  ITG *touch=NULL,i,j,k,np,idx,nu=0,ns=0,np2=0;
  double tot=0.,su=0.,ss=0.,sp=0.,v;

  NNEW(touch,ITG,nk);
  for(i=0;i<nk;i++) touch[i]=0;
  for(i=0;i<ne;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]=='U') np=6;
    else if(lakon[8*i]=='C') np=4;
    else continue;
    if(lakon[8*i]=='U'){
      for(k=0;k<np;k++){
        j=kon[ipkon[i]+k]-1;
        if((j>=0)&&(j<nk)) touch[j]|=1;            /* bit0: a live facet   */
      }
    }else{
      idx=mi0*i;
      if((i<ne0)&&(dam[idx]>1.+1.e-12)){
        for(k=0;k<np;k++){
          j=kon[ipkon[i]+k]-1;
          if((j>=0)&&(j<nk)) touch[j]|=2;          /* bit1: softening bulk */
        }
      }
      if((i<ne0)&&(nstate>0)&&(xstate[nstate*idx]>1.e-14)){
        for(k=0;k<np;k++){
          j=kon[ipkon[i]+k]-1;
          if((j>=0)&&(j<nk)) touch[j]|=4;          /* bit2: plastified     */
        }
      }
    }
  }
  for(i=0;i<nk;i++){
    if(touch[i]&1) nu++;
    if(touch[i]&2) ns++;
    if(touch[i]&4) np2++;
    for(j=1;j<mt;j++){
      k=nactdof[mt*i+j];
      if((k<=0)||(k>neq1)) continue;
      v=x[k-1]*x[k-1];
      tot+=v;
      if(touch[i]&1) su+=v;
      if(touch[i]&2) ss+=v;
      if(touch[i]&4) sp+=v;
    }
  }
  if(tot<=0.) tot=1.;
  printf("[WALLDIAG]   %s lives on: nodes with a live UC6 facet %.6f (%"
         ITGFORMAT " nodes), nodes on softening bulk %.6f (%" ITGFORMAT
         "), nodes on plastified bulk %.6f (%" ITGFORMAT
         "); |x|2=%.6e%s",tag,su/tot,nu,ss/tot,ns,sp/tot,np2,sqrt(tot),"\n");
  SFREE(touch);
}

static ITG damage_wall_setdiff(const ITG *cat,const double *xstate,
                               const double *xstateini,const double *dam,
                               const double *dambase,const double *visc,
                               const double *stx,
                               const ITG *ipkon,const char *lakon,
                               ITG ne0,ITG mi0,ITG nstate,ITG *nb)
{
  ITG i,j,nip,c,d,k,n=0;

  for(k=0;k<8;k++) nb[k]=0;
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    for(j=0;j<nip;j++){
      c=damage_ray_catof(xstate,xstateini,dam,dambase,visc,stx,
                         &lakon[8*i],i,j,mi0,nstate);
      d=c^cat[mi0*i+j];
      if(d==0) continue;
      n++;
      for(k=0;k<8;k++) if(d&(1<<k)) nb[k]++;
    }
  }
  return n;
}

static ITG damage_history_nip(const char *lakonel,ITG mi0)
{
  if((lakonel[6]=='L')&&(lakonel[7]=='C')) return mi0;
  if(strncmp(lakonel,"C3D20RB",7)==0) return mi0;
  if(strncmp(lakonel,"C3D8R",5)==0) return 1;
  if(strncmp(lakonel,"C3D8I",5)==0) return 8;
  if(strncmp(lakonel,"C3D20R",6)==0) return 8;
  if(strncmp(lakonel,"C3D20",5)==0) return 27;
  if(strncmp(lakonel,"C3D10",5)==0) return 4;
  if(strncmp(lakonel,"C3D4",4)==0) return 1;
  if(strncmp(lakonel,"C3D15",5)==0) return 9;
  if(strncmp(lakonel,"C3D6",4)==0) return 2;
  if(strncmp(lakonel,"C3D8",4)==0) return 8;
  return mi0;
}

/* Progressive damage material classifier shared by DE1 and DM2.0.
   Rice-Tracey + Evolution=Displacement keeps the historical four-constant
   signature.  DM2.0 is identified by model type 3 and a variable-length
   constant count 3+2*NPOINTS (NPOINTS>=2). */
static ITG damage_progressive_material(ITG imat,const ITG *ndmcon,
                                       const double *dmcon,ITG ndmat,
                                       ITG ntmat)
{
  ITG nconst,type,off;

  if(imat<1) return 0;
  nconst=ndmcon[2*(imat-1)];
  if((dmcon==NULL)||(ndmat<1)||(ntmat<1)) return 0;

  off=1+(ndmat+1)*ntmat*(imat-1);
  type=(ITG)dmcon[off];
  if((type==1)&&(nconst==4)) return 1;
  if((type==3)&&(nconst>=7)&&(((nconst-3)%2)==0)) return 1;

  return 0;
}

/* Detect actual progressive softening in the present Newton trial.  Merely
   having a DE1/DM2.0 material in the model is not enough: at least one active
   integration point must have D_trial>D_committed.  Comparing degradation D
   rather than the overloaded raw dam value also handles initiation crossing
   (omega<1 -> dam=1+D) without a false large jump. */
static ITG damage_de12_trial_softening(const double *dam,
                                      const double *dambase,
                                      const ITG *ipkon,const char *lakon,
                                      const ITG *ielmat,ITG mi2,
                                      const ITG *ndmcon,const double *dmcon,
                                      ITG ndmat,ITG ntmat,ITG ne0,ITG mi0,
                                      ITG *nsoft,double *maxdd)
{
  ITG i,j,nip,imat,elementsoft;
  double dtrial,dbase,dd;

  *nsoft=0;
  *maxdd=0.;
  if((dam==NULL)||(dambase==NULL)) return 0;

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='C') continue;

    imat=ielmat[mi2*i];
    if(!damage_progressive_material(imat,ndmcon,dmcon,ndmat,ntmat))
      continue;

    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    elementsoft=0;

    for(j=0;j<nip;j++){
      dtrial=dam[mi0*i+j]-1.;
      dbase=dambase[mi0*i+j]-1.;
      if(dtrial<0.) dtrial=0.;
      if(dbase<0.) dbase=0.;
      dd=dtrial-dbase;
      if(dd>DAMAGE_SLOW_NEWTON_D_TOL){
        elementsoft=1;
        if(dd>*maxdd) *maxdd=dd;
      }
    }
    if(elementsoft) (*nsoft)++;
  }

  return (*nsoft>0)?1:0;
}

/* Geometric estimate of the total iteration number required for value to
   reach target at the contraction measured over the latest Newton step.
   Invalid/non-contracting data deliberately returns a large sentinel so the
   stock cutback is not postponed. */
static ITG damage_slow_newton_estimate(ITG iit,double value,
                                       double previous,double target)
{
  double ratio,nextra,iest;

  if((!isfinite(value))||(!isfinite(previous))||(!isfinite(target))||
     (value<0.)||(previous<=0.)||(target<=0.))
    return DAMAGE_SLOW_NEWTON_INVALID_EST;
  if(value<=target) return iit;
  if(value>=previous) return DAMAGE_SLOW_NEWTON_INVALID_EST;

  ratio=value/previous;
  if((ratio<=0.)||(ratio>=1.)) return DAMAGE_SLOW_NEWTON_INVALID_EST;
  nextra=ceil(log(target/value)/log(ratio));
  if((!isfinite(nextra))||(nextra<0.)||
     (nextra>(double)DAMAGE_SLOW_NEWTON_INVALID_EST-iit))
    return DAMAGE_SLOW_NEWTON_INVALID_EST;

  iest=(double)iit+nextra;
  if(iest>(double)DAMAGE_SLOW_NEWTON_INVALID_EST)
    return DAMAGE_SLOW_NEWTON_INVALID_EST;
  return (ITG)iest;
}

/* Decide whether it is worth paying for more PARDISO factorizations.  The
   direct stock force and correction targets are used.  Requiring two
   consecutive contractions rejects residual oscillation and one-step noise;
   the hard cap and remaining-iteration cap bound the cost. */
static ITG damage_slow_newton_allow(ITG iit,const double *ram,
                                    const double *ram1,const double *ram2,
                                    const double *cam,const double *uam,
                                    double camprev1,double camprev2,
                                    const double *qa,const double *qam,
                                    const double *ctrl,ITG maxiters,
                                    ITG *iestres,ITG *iestcorr,
                                    ITG *iesttotal,double *rratio,
                                    double *cratio)
{
  ITG ip;
  double ea,c1,c2,targetres,targetcorr;

  *iestres=DAMAGE_SLOW_NEWTON_INVALID_EST;
  *iestcorr=DAMAGE_SLOW_NEWTON_INVALID_EST;
  *iesttotal=DAMAGE_SLOW_NEWTON_INVALID_EST;
  *rratio=0.;
  *cratio=0.;

  if((iit<3)||(iit>=maxiters)) return 0;
  if((ram[0]<=0.)||(ram1[0]<=0.)||(ram2[0]<=0.)||
     (cam[0]<=0.)||(camprev1<=0.)||(camprev2<=0.)||(uam[0]<=0.))
    return 0;

  /* Strict two-step monotonicity.  ram2 still contains the genuine
     two-iterations-old residual here; checkconvergence() has not yet folded
     it into its running minimum. */
  if(!((ram[0]<ram1[0])&&(ram1[0]<ram2[0])&&
       (cam[0]<camprev1)&&(camprev1<camprev2))) return 0;

  *rratio=ram[0]/ram1[0];
  *cratio=cam[0]/camprev1;

  ea=ctrl[23];
  ip=(ITG)ctrl[2];
  if(qa[0]>ea*qam[0]){
    c1=(iit<=ip)?ctrl[18]:ctrl[22];
    c2=ctrl[19];
  }else{
    c1=ea;
    c2=ctrl[24];
  }

  targetres=c1*qam[0];
  targetcorr=c2*uam[0];
  *iestres=damage_slow_newton_estimate(iit,ram[0],ram1[0],targetres);
  *iestcorr=damage_slow_newton_estimate(iit,cam[0],camprev1,targetcorr);
  *iesttotal=(*iestres>*iestcorr)?*iestres:*iestcorr;

  if(*iesttotal>maxiters) return 0;
  if(*iesttotal-iit>DAMAGE_SLOW_NEWTON_MAX_EXTRA) return 0;
  return 1;
}

/* Mark a bounded batch of DE1.2 C3D4 elements for terminal deletion.
   DE1 uses dam = 1 + D after initiation, so D is recovered locally without
   changing the public history layout.  Both Rice-Tracey DE1 and DM2.0
   tabulated ductile progressive materials are eligible.  The element is only
   marked here; remastruct and transactional
   commit/rollback remain under nonlingeo's existing A3 machinery. */
/* Terminal deletion can be restricted to named materials.

   The TP1 ladder needs its rungs to differ ONLY in which phase is allowed
   to erode, with the constitutive routines untouched:

     CCX_DAMAGE_DELETE_MAT=NONE   continuous damage, no topology change
     CCX_DAMAGE_DELETE_MAT=ZrH    only the hydride erodes
     CCX_DAMAGE_DELETE_MAT=ALL    stock behaviour (default when unset)

   The list is comma separated and matched case-insensitively against the
   *MATERIAL names.  A rejected element is removed from the terminal scan
   only; its damage keeps evolving exactly as before, so the comparison
   isolates the topology change and nothing else. */

static ITG damage_delete_allowed(ITG imat,const char *matname,
                                 const char *filter)
{
  const char *p,*q;
  char nm[81];
  ITG i,n;

  if(filter==NULL) return 1;
  if((strcmp(filter,"ALL")==0)||(strcmp(filter,"all")==0)) return 1;
  if((strcmp(filter,"NONE")==0)||(strcmp(filter,"none")==0)) return 0;
  if((matname==NULL)||(imat<1)) return 1;

  n=0;
  for(i=0;i<80;i++){
    if(matname[80*(imat-1)+i]==' ') break;
    nm[i]=matname[80*(imat-1)+i];
    n++;
  }
  nm[n]=0;
  if(n==0) return 1;

  p=filter;
  while(*p!=0){
    q=strchr(p,',');
    if(q==NULL) q=p+strlen(p);
    if((ITG)(q-p)==n){
      for(i=0;i<n;i++){
        if(tolower((unsigned char)p[i])!=tolower((unsigned char)nm[i])) break;
      }
      if(i==n) return 1;
    }
    p=(*q==0)?q:q+1;
  }
  return 0;
}

/* Read one coefficient back out of the assembled sparse operator.

   Convention taken from add_sm_st_as.f, which is what mafillsm and
   mafilldamas write through: ad(i) holds the diagonal; a coefficient with
   i>j lives in column j of au; one with i<j lives in column i, offset by
   nzs(3).  Row indices inside a column are sorted, so the search is a
   bisection, exactly as nident does on the way in.

   Returns 0 for a structurally absent coefficient, which is the correct
   value - the sparsity pattern is a superset of the assembled entries and
   a missing slot means the two degrees of freedom share no element. */

static double damage_fd_coeff(const double *ad,const double *au,
                              const ITG *jq,const ITG *irow,const ITG *nzs,
                              ITG i,ITG j)
{
  ITG col,want,lo,hi,mid,off;

  if(i==j) return ad[i-1];

  if(i>j){col=j;want=i;off=0;}
  else   {col=i;want=j;off=nzs[2];}

  lo=jq[col-1];hi=jq[col]-1;
  while(lo<=hi){
    mid=(lo+hi)/2;
    if(irow[mid-1]==want) return au[mid-1+off];
    if(irow[mid-1]<want) lo=mid+1; else hi=mid-1;
  }
  return 0.;
}


/* Deletes a DEAD element that is the SOLE support of a node.

   Measured configuration (E-61): on m14_fine node 776 had exactly one live
   element, that element stood at D = 1.0000, and the node travelled 9.63 mm
   while the other three nodes of the same element moved 0.44.  The element's
   longest edge went from 0.0974 to 9.809 - a stretch of 101x on a 4 mm
   specimen.  An element below a per cent of its stiffness carries almost
   nothing, so a node whose whole support is one such element is very nearly
   free and Newton solves for it.

   Holding that node with a diagonal term was tried and does not cure it
   (E-61): the added stiffness would have to be non-perturbative to compete
   with the element itself.  Removing the dead element instead is physically
   near-free - it was carrying under 1% - and it hands the node to the
   existing BK4 / damfloat machinery, which is built for a node that has lost
   all its bulk.

   This is NOT damdangle.  damdangle counted support without ever reading
   degradation, so it deleted healthy load-bearing material and emptied small
   meshes (E-22).  Here the element must itself be dead, and it must be the
   only thing a node has.  Both conditions are necessary and both are
   measured, not assumed.

   The batch is bounded exactly like the terminal batch, and marking uses the
   same ipkon -> -ipkon-2 convention. */
static ITG damage_de13_mark_deadsole(const double *dam,const double *visc,
                                     ITG usevisc,ITG *ipkon,const char *lakon,
                                     const ITG *kon,ITG nk,ITG ne0,ITG mi0,
                                     double gdead,ITG batchmax)
{
  ITG i,j,n,nip,nnew=0,usedam,*nlive=NULL;
  const double *src;
  double dmx,g;

  if((nk<=0)||(batchmax<=0)) return 0;
  usedam=((usevisc&&(visc!=NULL))?0:1);
  src=usedam?dam:visc;

  NNEW(nlive,ITG,nk);
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    /* kon is 0-based here: the nodes of element i are kon[ipkon[i]+0] to
       kon[ipkon[i]+nope-1].  frd.c:1684 takes the last node as
       kon[ipkon[i]+nope-1], and the VTK writer below (which uses j=0..3)
       reproduces the deck connectivity exactly - 27359 of 27359 cells on
       m12_eta15.  An earlier version of this loop ran j=1..4, which skips
       node 1 and picks up node 1 of the NEXT element instead. */
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n>=0)&&(n<nk)) nlive[n]++;
    }
  }

  for(i=0;i<ne0;i++){
    if(nnew>=batchmax) break;
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    dmx=0.;
    for(j=0;j<nip;j++){
      /* dam holds 1+D (the integer part flags initiation); visc holds D
         directly.  Reading dam raw gives 1.0024 for D=0.0024 and deletes a
         healthy element - measured, and it is exactly how damdangle went
         wrong.  mark_terminal decodes it the same way. */
      double dd=usedam?(src[mi0*i+j]-1.):src[mi0*i+j];
      if(dd<0.) dd=0.;
      if(dd>1.) dd=1.;
      if(dd>dmx) dmx=dd;
    }
    g=1.-dmx; if(g<0.) g=0.;
    if(g>=gdead) continue;
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n<0)||(n>=nk)) continue;
      if(nlive[n]==1){
        printf("[DAMAGE DEADSOLE]   element %" ITGFORMAT " g=%.6e sole "
               "support of node %" ITGFORMAT "\n",
               i+1,g,n+1);
        ipkon[i]=-ipkon[i]-2;
        nnew++;
        /* the node counts are now stale for this element's nodes; drop them
           so a second element at the same node cannot also be taken in the
           same pass */
        for(j=0;j<4;j++){
          n=kon[(-ipkon[i]-2)+j]-1;
          if((n>=0)&&(n<nk)) nlive[n]=0;
        }
        break;
      }
    }
  }
  SFREE(nlive);
  return nnew;
}

/* DEADALL - the Class A artefact (E-69).
 *
 * DEADSOLE above asks whether a dead element is the SOLE support of a node.
 * The three runs that died by divergence died on a node with TWO, not one:
 *
 *   uf50   node 2471  2 live tets, 0 facets, both D=1.0000, edge stretch 22.2
 *   eta1e4 node 3354  2 live tets, 0 facets, both D=1.0000, edge stretch 49.8
 *   eta15  node 1581  2 live tets, 0 facets, both D=1.0000, edge stretch 35.5
 *
 * In each of those models exactly ONE node in the whole mesh had its entire
 * live support dead - out of 17, 145 and 65 low-support nodes respectively -
 * and it was the node the solver threw.  Three out of three, with a
 * selectivity of one in seventeen to one in a hundred and forty-five.
 *
 * The node is a free swinging point: everything holding it carries only gmin,
 * so Newton solves it almost unconstrained and draws the elements out into a
 * needle.  Note that a needle PRESERVES VOLUME - it stretches along one
 * direction and collapses across the other two - so V/V0 does not see it
 * (1.2 on eta15 while an edge went 0.120 -> 4.151 mm).  Do not look for this
 * with a volume ratio.
 *
 * SAFETY.  A node qualifies only when every live element at it is dead, so
 * every element this routine deletes is itself dead: element i is live and
 * touches a qualifying node n, and n qualifies only if all its live elements
 * are dead.  Nothing load-bearing can be removed.  That is the property
 * damdangle did not have (E-22), and it is what makes this different.
 *
 * The facet guard is the second half.  A node still tied to the other side of
 * an interface is not free, whatever its bulk looks like, and the
 * cohesive-only node is Class B - a CONDITIONING defect (E-67, E-68) that must
 * not be answered by deleting material.  286 of 318 cohesive-only nodes on
 * m12_epsf50 are perfectly well supported.
 */
static ITG damage_de13_mark_deadall(const double *dam,const double *visc,
                                    ITG usevisc,ITG *ipkon,const char *lakon,
                                    const ITG *kon,ITG nk,ITG ne,ITG ne0,
                                    ITG mi0,double gdead,ITG batchmax,
                                    ITG *nnodes)
{
  ITG i,j,n,nip,nope,nnew=0,usedam;
  ITG *nlive=NULL,*ndead=NULL,*nfac=NULL,*take=NULL;
  const double *src;
  double dmx,g;

  if(nnodes!=NULL) *nnodes=0;
  if((nk<=0)||(batchmax<=0)) return 0;
  usedam=((usevisc&&(visc!=NULL))?0:1);
  src=usedam?dam:visc;

  NNEW(nlive,ITG,nk);
  NNEW(ndead,ITG,nk);
  NNEW(nfac,ITG,nk);

  /* live cohesive support, over the whole element range: user elements carry
     their node count in byte 7 of lakon, the idiom the node dump uses. */
  for(i=0;i<ne;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='U') continue;
    nope=(ITG)((unsigned char)lakon[8*i+7]);
    if((nope<1)||(nope>20)) continue;
    for(j=0;j<nope;j++){
      n=kon[ipkon[i]+j]-1;
      if((n>=0)&&(n<nk)) nfac[n]++;
    }
  }

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;
    dmx=0.;
    for(j=0;j<nip;j++){
      /* dam holds 1+D, visc holds D - same decode as mark_deadsole. */
      double dd=usedam?(src[mi0*i+j]-1.):src[mi0*i+j];
      if(dd<0.) dd=0.;
      if(dd>1.) dd=1.;
      if(dd>dmx) dmx=dd;
    }
    g=1.-dmx; if(g<0.) g=0.;
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n<0)||(n>=nk)) continue;
      nlive[n]++;
      if(g<gdead) ndead[n]++;
    }
  }

  NNEW(take,ITG,nk);
  for(n=0;n<nk;n++){
    if(nfac[n]>0) continue;
    if(nlive[n]<1) continue;
    if(ndead[n]!=nlive[n]) continue;
    take[n]=1;
    if(nnodes!=NULL) (*nnodes)++;
  }

  /* The qualifying set is a SNAPSHOT.  Deleting an element also takes support
     away from its other three nodes, and resolving that here would need a
     fixed point; the next increment picks it up instead, which is how every
     other batch in this file behaves. */
  for(i=0;i<ne0;i++){
    if(nnew>=batchmax) break;
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
    for(j=0;j<4;j++){
      n=kon[ipkon[i]+j]-1;
      if((n<0)||(n>=nk)) continue;
      if(take[n]==0) continue;
      printf("[DAMAGE DEADALL]   element %" ITGFORMAT " deleted: node %"
             ITGFORMAT " has %" ITGFORMAT " live element(s), all dead, "
             "no cohesive facet\n",i+1,n+1,nlive[n]);
      ipkon[i]=-ipkon[i]-2;
      nnew++;
      break;
    }
  }

  SFREE(take);SFREE(nfac);SFREE(ndead);SFREE(nlive);
  return nnew;
}

static ITG damage_de13_mark_terminal(double *dam,ITG *ipkon,
                                     const char *lakon,const ITG *ielmat,
                                     ITG mi2,const ITG *ndmcon,
                                     const double *dmcon,ITG ndmat,ITG ntmat,
                                     ITG ne0,ITG mi0,double ddelete,
                                     ITG batchmax,double *batch_dmax,
                                     double *trigger_value,ITG *trigger_ip,
                                     const char *matname,
                                     const char *delfilter,
                                     const double *visc,ITG usevisc,
                                     double *batch_vmin)
{
  ITG i,j,nip,imat,nnew=0;
  double de,demax,dv,dvmax,dtrig;

  *batch_dmax=0.;
  if(batch_vmin!=NULL) *batch_vmin=1.;

  for(i=0;i<ne0;i++){
    if(nnew>=batchmax) break;
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;

    imat=ielmat[mi2*i];
    if(imat<1) continue;
    if(!damage_progressive_material(imat,ndmcon,dmcon,ndmat,ntmat)) continue;
    if(!damage_delete_allowed(imat,matname,delfilter)) continue;

    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;

    demax=0.;
    dvmax=0.;
    for(j=0;j<nip;j++){
      de=dam[mi0*i+j]-1.;
      if(de<0.) de=0.;
      if(de>1.) de=1.;
      if(de>demax) demax=de;
      if(visc!=NULL){
        dv=visc[mi0*i+j];
        if(dv<0.) dv=0.;
        if(dv>1.) dv=1.;
        if(dv>dvmax) dvmax=dv;
      }
    }

    /* The stress is scaled by 1-Dvis, not by 1-D: resultsmech.f uses
       damvisc whenever the viscosity is on.  Triggering deletion on D
       therefore removes an element that is still carrying 1-Dvis of its
       effective stress, and releases that force in one increment at
       constant load.  It also explains why cutting the step never helped
       and sometimes hurt: beta=dt/(eta+dt), so a smaller step makes Dvis
       lag further behind D and the deleted element carries more.
       usevisc makes the trigger read the same variable the stress does. */
    dtrig=demax;
    if((usevisc==1)&&(visc!=NULL)) dtrig=dvmax;

    if(dtrig>=ddelete){
      if(batch_vmin!=NULL){
        if(dvmax<*batch_vmin) *batch_vmin=dvmax;
      }
      /* DE1.3.1: capture the terminal trigger at the instant the element
         first changes topology.  Later same-load redistribution can alter
         dam for already deleted elements, so .damage must not reconstruct
         this value from the current constitutive state. */
      if((trigger_value!=NULL)&&(trigger_ip!=NULL)){
        trigger_value[i]=demax;
        trigger_ip[i]=1;
      }
      ipkon[i]=-ipkon[i]-2;
      nnew++;
      if(demax>*batch_dmax) *batch_dmax=demax;
    }
  }

  return nnew;
}


/* Exact DE1 element statistics.
   For each element the maximum degradation over its active integration
   points is used.  In the present DE1 implementation only C3D4 is enabled,
   therefore this is exactly the single integration-point value. */
static void damage_de1_stats(const double *dam,const double *damold,
                             const ITG *ipkon,const char *lakon,
                             ITG ne0,ITG mi0,
                             ITG *nactive,ITG *ngt01,ITG *ngt05,
                             ITG *ngt09,ITG *nfull,ITG *nchanged,
                             double *dmax,double *maxdelta)
{
  ITG i,j,nip;
  double de,dold,delta,demax,delmax;

  *nactive=0;*ngt01=0;*ngt05=0;*ngt09=0;*nfull=0;*nchanged=0;
  *dmax=0.;*maxdelta=0.;

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(lakon[8*i]!='C') continue;

    nip=damage_history_nip(&lakon[8*i],mi0);
    if(nip<1) nip=1;
    if(nip>mi0) nip=mi0;

    demax=0.;delmax=0.;
    for(j=0;j<nip;j++){
      de=dam[mi0*i+j]-1.;
      if(de<0.) de=0.;
      if(de>1.) de=1.;
      if(de>demax) demax=de;

      if(damold!=NULL){
        dold=damold[mi0*i+j]-1.;
        if(dold<0.) dold=0.;
        if(dold>1.) dold=1.;
        delta=fabs(de-dold);
        if(delta>delmax) delmax=delta;
      }
    }

    if(demax>0.){
      (*nactive)++;
      if(demax>0.1) (*ngt01)++;
      if(demax>0.5) (*ngt05)++;
      if(demax>0.9) (*ngt09)++;
      if(demax>=0.999) (*nfull)++;
      if(demax>*dmax) *dmax=demax;
    }

    if(delmax>DAMAGE_DE1_FP_TOL) (*nchanged)++;
    if(delmax>*maxdelta) *maxdelta=delmax;
  }
}

/* Append one compact accepted-state record.  This file is intentionally
   independent of the legacy .damage hard-deletion history. */
static void damage_de1_append_stats(const char *jobnamec,ITG istep,ITG iinc,
                                    double steptime,double totaltime,
                                    ITG passes,ITG nactive,ITG ngt01,
                                    ITG ngt05,ITG ngt09,ITG nfull,
                                    double dmax,double maxdelta)
{
  char fname[200]="";
  FILE *f=NULL;
  static ITG de1stats_initialized=0;

  strcpy2(fname,jobnamec,132);
  strcat(fname,".de1stats");

  /* DE1.3.1: one solver process == one fresh diagnostic history.
     The first accepted damage state truncates stale data from an older run;
     subsequent accepted states append normally. */
  if(de1stats_initialized==0){
    f=fopen(fname,"w");
  }else{
    f=fopen(fname,"a");
  }
  if(f==NULL) return;

  if(de1stats_initialized==0){
    fprintf(f,"# CalculiX DE1.3.1 accepted damage states\n");
    fprintf(f,"# step increment step_time total_time passes active "
              "Dgt0.1 Dgt0.5 Dgt0.9 Dfull Dmax max_dD\n");
    de1stats_initialized=1;
  }

  fprintf(f,"%" ITGFORMAT " %" ITGFORMAT " %.15e %.15e "
            "%" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT " "
            "%" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT " "
            "%.15e %.15e\n",
          istep,iinc,steptime,totaltime,passes,nactive,ngt01,ngt05,
          ngt09,nfull,dmax,maxdelta);
  fclose(f);
}

/* Write the latest accepted DE1 state as an exact element-cell VTK snapshot.
   The file is overwritten, so long calculations do not accumulate large
   post-processing files.  All active C3D4 cells are written.  DE1_D and
   DUCT_IP are CELL_DATA taken directly from the solver integration-point
   history; no extrapolation or nodal averaging is involved.  Coordinates
   are written in the current deformed configuration. */
static void damage_de1_write_vtk(const char *jobnamec,
                                 const double *co,const double *vold,
                                 ITG nk,ITG mt,const ITG *kon,
                                 const ITG *ipkon,const char *lakon,
                                 const ITG *ielmat,ITG mi2,
                                 const double *dam,ITG mi0,ITG ne0,
                                 ITG istep,ITG iinc,double steptime)
{
  char fname[200]="",seq[32]="";
  FILE *f=NULL;
  ITG i,j,indexe,node,ncell=0;
  double d,duct,x,y,z;
  static ITG vtkseries=-1,vtkcount=0;

  /* CCX_DAMAGE_VTK_SERIES keeps every accepted state as
     <job>.de1.NNNNN.vtk instead of overwriting a single frame.  Needed to
     show a sequence of damage rather than one final picture. */

  if(vtkseries<0){
    vtkseries=(getenv("CCX_DAMAGE_VTK_SERIES")!=NULL)?1:0;
  }

  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)==0) ncell++;
  }

  strcpy2(fname,jobnamec,132);
  if(vtkseries==1){
    sprintf(seq,".de1.%5.5d.vtk",(int)vtkcount);
    strcat(fname,seq);
    vtkcount++;
  }else{
    strcat(fname,".de1.vtk");
  }
  f=fopen(fname,"w");
  if(f==NULL) return;

  fprintf(f,"# vtk DataFile Version 3.0\n");
  fprintf(f,"CalculiX DE1.3.1 exact cell damage step=%" ITGFORMAT
            " inc=%" ITGFORMAT " time=%.12e\n",istep,iinc,steptime);
  fprintf(f,"ASCII\n");
  fprintf(f,"DATASET UNSTRUCTURED_GRID\n");

  fprintf(f,"POINTS %" ITGFORMAT " double\n",nk);
  for(i=0;i<nk;i++){
    x=co[3*i];
    y=co[3*i+1];
    z=co[3*i+2];
    if(vold!=NULL){
      x+=vold[mt*i+1];
      y+=vold[mt*i+2];
      z+=vold[mt*i+3];
    }
    fprintf(f,"%.15e %.15e %.15e\n",x,y,z);
  }

  fprintf(f,"CELLS %" ITGFORMAT " %" ITGFORMAT "\n",ncell,5*ncell);
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    indexe=ipkon[i];
    fprintf(f,"4");
    for(j=0;j<4;j++){
      node=kon[indexe+j]-1;
      fprintf(f," %" ITGFORMAT,node);
    }
    fprintf(f,"\n");
  }

  fprintf(f,"CELL_TYPES %" ITGFORMAT "\n",ncell);
  for(i=0;i<ncell;i++) fprintf(f,"10\n");

  fprintf(f,"CELL_DATA %" ITGFORMAT "\n",ncell);

  fprintf(f,"SCALARS DE1_D double 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    d=dam[mi0*i]-1.;
    if(d<0.) d=0.;
    if(d>1.) d=1.;
    fprintf(f,"%.15e\n",d);
  }

  fprintf(f,"SCALARS DUCT_IP double 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    duct=dam[mi0*i];
    fprintf(f,"%.15e\n",duct);
  }

  fprintf(f,"SCALARS ELEMENT_ID int 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    fprintf(f,"%" ITGFORMAT "\n",i+1);
  }

  fprintf(f,"SCALARS MATERIAL_ID int 1\nLOOKUP_TABLE default\n");
  for(i=0;i<ne0;i++){
    if(ipkon[i]<0) continue;
    if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
    fprintf(f,"%" ITGFORMAT "\n",ielmat[mi2*i]);
  }

  fclose(f);
}

void nonlingeo(double **cop,ITG *nk,ITG **konp,ITG **ipkonp,char **lakonp,
	       ITG *ne,
	       ITG *nodeboun,ITG *ndirboun,double *xboun,ITG *nboun,
	       ITG **ipompcp,ITG **nodempcp,double **coefmpcp,char **labmpcp,
	       ITG *nmpc,
	       ITG *nodeforc,ITG *ndirforc,double *xforc,ITG *nforc,
	       ITG **nelemloadp,char **sideloadp,double *xload,ITG *nload,
	       ITG *nactdof,
	       ITG **icolp,ITG *jq,ITG **irowp,ITG *neq,ITG *nzl,
	       ITG *nmethod,ITG **ikmpcp,ITG **ilmpcp,ITG *ikboun,
	       ITG *ilboun,
	       double *elcon,ITG *nelcon,double *rhcon,ITG *nrhcon,
	       double *alcon,ITG *nalcon,double *alzero,ITG **ielmatp,
	       ITG **ielorienp,ITG *norien,double *orab,ITG *ntmat_,
	       double *t0,double *t1,double *t1old,
	       ITG *ithermal,double *prestr,ITG *iprestr,
	       double **voldp,ITG *iperturb,double *sti,ITG *nzs, 
	       ITG *kode,char *filab,
	       ITG *idrct,ITG *jmax,ITG *jout,double *timepar,
	       double *eme,
	       double *xbounold,double *xforcold,double *xloadold,
	       double *veold,double *accold,
	       char *amname,double *amta,ITG *namta,ITG *nam,
	       ITG *iamforc,ITG **iamloadp,
	       ITG *iamt1,double *alpha,ITG *iexpl,
	       ITG *iamboun,double *plicon,ITG *nplicon,double *plkcon,
	       ITG *nplkcon,
	       double **xstatep,ITG *npmat_,ITG *istep,double *ttime,
	       char *matname,double *qaold,ITG *mi,
	       ITG *isolver,ITG *ncmat_,ITG *nstate_,
	       double *cs,ITG *mcs,ITG *nkon,double **enerp,ITG *mpcinfo,
	       char *output,
	       double *shcon,ITG *nshcon,double *cocon,ITG *ncocon,
	       double *physcon,ITG *nflow,double *ctrl,
	       char *set,ITG *nset,ITG *istartset,
	       ITG *iendset,ITG *ialset,ITG *nprint,char *prlab,
	       char *prset,ITG *nener,ITG *ikforc,ITG *ilforc,double *trab,
	       ITG *inotr,ITG *ntrans,double **fmpcp,char *cbody,
	       ITG *ibody,double *xbody,ITG *nbody,double *xbodyold,
	       ITG *ielprop,double *prop,ITG *ntie,char *tieset,
	       ITG *itpamp,ITG *iviewfile,char *jobnamec,double *tietol,
	       ITG *nslavs,double *thicke,ITG *ics,
	       ITG *nintpoint,ITG *mortar,ITG *ifacecount,char *typeboun,
	       ITG **islavsurfp,double **pslavsurfp,double **clearinip,
	       ITG *nmat,double *xmodal,ITG *iaxial,ITG *inext,ITG *nprop,
	       ITG *network,char *orname,double *vel,ITG *nef,
	       double *velo,double *veloo,double *energy,ITG *itempuser,
	       ITG *ipobody,ITG *inewton,double *t0g,double *t1g,
	       ITG *ifreebody,ITG *nlabel,ITG *ndmat_,ITG *ndmcon,
	       double *dmcon,double *dam){

  char description[13]="            ",*lakon=NULL,jobnamef[396]="",
    *sideface=NULL,*labmpc=NULL,*lakonf=NULL,*env,*envsys,fneig[132]="",
    *sideloadref=NULL,*sideload=NULL,stiffmatrix[132]="",
    damagefilename[160]="",*sideloadf=NULL,cflag[1]=" ",
    *damage_tangent_env=NULL,*damage_reeq_scale_env=NULL,
    *damage_linesearch_env=NULL,*damage_topology_env=NULL,
    *damage_de13_env=NULL,*damage_delete_filter=NULL,
    *damage_visc_env=NULL,*damage_diss_env=NULL,
    *damage_fracture_env=NULL,
    *damage_fracture_seta=NULL,*damage_fracture_setb=NULL;
  char damage_fracture_a[81],damage_fracture_b[81]; 
 
  ITG *inum=NULL,k,l,iout=0,icntrl,iinc=0,jprint=0,iit=-1,jnz=0,
    icutb=0,istab=0,uncoupled,n1,n2,itruecontact=1,iclean=0,
    iperturb_sav[2],iforbou,*icol=NULL,*irow=NULL,ielas=0,icmd=0,
    memmpc_,mpcfree,icascade,maxlenmpc,*nodempc=NULL,*iaux=NULL,
    *nodempcref=NULL,memmpcref_,mpcfreeref,*itg=NULL,*ineighe=NULL,
    *ieg=NULL,ntg=0,ntr,*kontri=NULL,*nloadtr=NULL,idamping=0,
    *ipiv=NULL,ntri,newstep,mode=-1,noddiam=-1,nasym=0,im,
    ntrit,*inocs=NULL,*nacteq=NULL,*ipface=NULL,masslesslinear=0,
    *nactdog=NULL,nteq,*itietri=NULL,*koncont=NULL,istrainfree=0,
    ncont,ne0,nkon0,*ipkon=NULL,*kon=NULL,*ielorien=NULL,
    *ielmat=NULL,itp=0,symmetryflag=0,inputformat=0,kscale=1,
    *iruc=NULL,iitterm=0,iturbulent,ngraph=1,ismallsliding=0,
    *ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,i0ref,irref,icref,
    *itiefac=NULL,*islavsurf=NULL,*islavnode=NULL,*imastnode=NULL,
    *nslavnode=NULL,*nmastnode=NULL,*imastop=NULL,imat,
    *iponoels=NULL,*inoels=NULL,*islavsurfold=NULL,maxlenmpcref,
    *islavact=NULL,mt=mi[1]+1,*nactdofinv=NULL,*ipe=NULL, 
    *ime=NULL,*ikactmech=NULL,nactmech,inode,idir,neold,neini,
    iemchange=0,nzsrad,*mast1rad=NULL,*irowrad=NULL,*icolrad=NULL,
    *jqrad=NULL,*ipointerrad=NULL,*integerglob=NULL,negpres=0,
    mass[2]={0,0},stiffness=1, buckling=0, rhsi=1, intscheme=0,idiscon=0,
    coriolis=0,*ipneigh=NULL,*neigh=NULL,maxprevcontel,nslavs_prev_step,
    *nelemface=NULL,*ipoface=NULL,*nodface=NULL,*ifreestream=NULL,
    *isolidsurf=NULL,*neighsolidsurf=NULL,*iponoeln=NULL,*inoeln=NULL,
    nface,nfreestream,nsolidsurf,i,icfd=0,id,nslavquadel=0,
    node,networknode,iflagact=0,*nodorig=NULL,*ipivr=NULL,iglob=0,
    *inomat=NULL,ntrimax,*nx=NULL,*ny=NULL,*nz=NULL,nforcrhs,nloadrhs,
    idampingwithoutcontact=0,*nactdoh=NULL,*nactdohinv=NULL,*ipkonf=NULL,
    *ielmatf=NULL,*ielorienf=NULL,ialeatoric=0,nloadref,isym,
    *nelemloadref=NULL,*iamloadref=NULL,*idefload=NULL,nload_,
    *nelemload=NULL,*iamload=NULL,ncontacts=0,inccontact=0,nrhs=1,
    j=0,inoelnsize=0,isensitivity=0,*konf=NULL,nbodyrhs,
    *iwork=NULL,nelt,lrgw,*igwk=NULL,itol,itmax,iter,ierr,iunit,ligw,
    mei[4]={0,0,0,0},*itreated=NULL,mscalmethod=-1,inoelfree,
    isiz=0,num_cpus,sys_cpus,ne1d2d=0,kchdep,nkftot,
    ifreesurface=0,*iponoelf=NULL,*inoelf=NULL,*iponoel=NULL,
    *damage_iponoel_trial=NULL,*damage_orphan_seen=NULL,*damage_addok=NULL,
    mortartrafoflag=0,*nelold=NULL,*nelnew=NULL,*nkold=NULL,*nknew=NULL,
    *ipompcf=NULL,*nodempcf=NULL,*nodebounf=NULL,*ndirbounf=NULL,
    *nelemloadf=NULL,*ipobodyf=NULL,nkf,nkonf,memmpcf,nbounf,nloadf,nmpcf,
    *ikbounf=NULL,*ilbounf=NULL,*ikmpcf=NULL,*ilmpcf=NULL,*iambounf=NULL,
    *iamloadf=NULL,*inotrf=NULL,*jqw=NULL,*iroww=NULL,nzsw,*jqtherm=NULL,
    *kslav=NULL,*lslav=NULL,*ktot=NULL,*ltot=NULL,nmasts,neqtot,
    intpointvarm,calcul_fn,calcul_f,calcul_qa,calcul_cauchy,ikin,
    intpointvart,*jqbi=NULL,*irowbi=NULL,*jqib=NULL,*irowib=NULL,
    idispfrdonly,*inumcp=NULL,nmethodold=*nmethod,
    idamage=0,iitsav=0,idamagereeq=0,ilocalsubstep=0,*ipkondamageini=NULL,
    *damage_tent_elem=NULL,*damage_tent_mat=NULL,*damage_tent_ip=NULL,
    *damage_de13_trigger_ip=NULL,*damage_ract=NULL,
    damage_ray_probe=0,damage_ray_shots=0,damage_ray_max=8,
    *damage_ray_cat=NULL,damage_bt_mode=0,damage_bt_ntrial=0,
    damage_rescue_mode=0,damage_rescue_bt_on=0,damage_rescue_used=0,
    damage_rescue_nfired=0,damage_rescue_nok=0,damage_rescue_maxlevel=1,
    damage_evt_on=0,damage_evt_nstep=0,*damage_evt_sgn=NULL,
    damage_evt_nsw[8]={0,0,0,0,0,0,0,0},
    damage_evt_fe[8]={0,0,0,0,0,0,0,0},
    damage_evt_fp[8]={0,0,0,0,0,0,0,0},
    damage_reg_on=0,damage_reg_nlam=0,damage_reg_napply=0,
    damage_reg_level=0,damage_rec_window=5,damage_rec_maxunrec=3,
    damage_corr_mode=0,damage_corr_on=0,damage_corr_ninc=0,
    damage_corr_nint=0,damage_corr_nfact=0,damage_corr_clean=0,
    damage_corr_try=0,damage_corr_maxinc=60,damage_corr_exit=5,
    damage_corr_tryevery=5,damage_corr_grace=10,damage_corr_dtn=0,
    damage_corr_trial=0,damage_corr_nsince=0,damage_corr_stableneed=3,
    damage_corr_nwallstab=0,damage_corr_maxesc=3,damage_corr_nwalltot=0,
    damage_corr_maxwall=30,damage_dth_n=0,damage_dth_i=0,
    damage_rec_healthy=0,damage_rec_unrec=0,damage_rec_used_in_inc=0,
    damage_rec_disarmed=0,
    damage_ct_mode=0,damage_ct_on=0,damage_ct_have=0,damage_ct_nring=0,
    damage_ct_elem=-1,damage_ct_ip=-1,damage_ct_nsupp=0,
    damage_ct_supp[18],damage_ct_node[18],damage_ct_dir[18],
    damage_ct_ncommit=0,damage_ct_nstep=0,damage_ct_ncorr=0,
    damage_ct_nfact=0,damage_ct_neval=0,damage_ct_nretry=0,
    damage_ct_nkink=0,damage_ct_partial=0,damage_ct_nbl=0,
    *damage_ct_bl=NULL,*damage_ct_fl=NULL,damage_ct_maxstep=60,
    damage_ct_maxcorr=15,damage_ct_maxfact=700,damage_ct_maxeval=1600,
    damage_ct_maxretry=5,damage_ct_maxkink=2,damage_ct_lastminc=0,
    damage_ct_refused=0,damage_ct_p1=0,damage_ct_p2=0,damage_ct_p3=0,
    damage_ct_head=0,damage_ct_alloc=0,damage_ct_arm=0,damage_ct_step=0,
    damage_ct_it=0,damage_ct_pred=0,*damage_ct_sgn=NULL,
    damage_ct_newstep=1,damage_ct_epsok=0,damage_ct_p1r=0,
    damage_ct_used=0,
    damage_ct_p2r=0,damage_ct_p3r=0,damage_ct_nprog=0,
    damage_dl_mode=0,damage_dl_on=0,damage_dl_have=0,damage_dl_used=0,
    damage_dl_maxtrial=6,damage_dl_maxeval=600,damage_dl_maxfact=250,
    damage_dl_neval=0,damage_dl_nfact=0,damage_dl_narm=0,
    damage_dl_maxarm=12,damage_dl_nacc=0,damage_dl_nrej=0,
    damage_dl_nnewt=0,damage_dl_ncau=0,damage_dl_ndog=0,
    damage_dl_nfail=0,damage_dl_banner=0,damage_dl_incarm=0,
    damage_dl_selfrec=0,damage_dl_lincheck=0,damage_dl_lc_due=0,
    damage_dl_lc_it=1,damage_dl_lc_nit=1,damage_ls_probe=0,
    damage_wall_armed=0,*damage_wall_cat=NULL,damage_wall_nb[8],
    damage_wall_null=0,
    damage_dl_lasthelp=0,damage_dl_recdone=-1,
    damage_aba_mode=0,damage_aba_done=0,damage_bt_nring=0,
    damage_aba_ninc=0,damage_aba_hit=0,damage_aba_inc[4]={0,0,0,0},
    damage_ray_ninc=0,damage_ray_incok=1,damage_ray_inc[4]={0,0,0,0},
    damage_bt_window=1,
    damage_release_probe=0,damage_release_armed=0,damage_release_pass=0,
    damage_de13_term_only=0,damage_release_rebuild=0,
    damage_release_nterm=0,damage_release_nother=0,
    damage_release_nisl=0,damage_release_ncoh=0,damage_release_iforbou=0,
    damage_tent_count=0,damage_tent_step=0,damage_tent_increment=0,
    damage_batch=0,damage_scan_count=0,damage_nip_local=0,
    damage_mode=0,damage_predict_count=0,damage_event_cut=0,
    damage_active_pass=0,damage_soft_reeq=0,damage_fast_retry=0,damage_fast_used=0,
    damage_fast_recover=0,damage_fast_failures=0,
    damage_de1_nactive=0,damage_de1_gt01=0,damage_de1_gt05=0,
    damage_de1_gt09=0,damage_de1_nfull=0,damage_de1_nchanged=0,
    damage_de12_enabled=0,damage_de12_matcount=0,damage_dm20_matcount=0,
    damage_tangent_mode=0,damage_reeq_scale_mode=0,
    damage_linesearch_mode=0,damage_linesearch_active=0,
    damage_linesearch_applied=0,damage_linesearch_nsoft=0,
    damage_linesearch_trial=0,damage_linesearch_contracted=0,
    damage_topology_deferred_mode=0,damage_topology_rebuild=1,
    damage_unsym_active=0,damage_unsym_elems=0,damage_unsym_report=0,
    damage_unsym_hole=0,damage_unsym_floor=0,damage_unsym_holerep=0,
    damage_unsym_census=0,damage_unsym_tanfull=0,
    damage_snap_elem=0,damage_snap_bad=0,
    damage_diss_report=0,damage_diss_init=0,damage_diss_ctrl=0,
    damage_diss_have=0,damage_diss_ok=0,damage_diss_engaged=0,
    damage_diss_probe=0,damage_batch_list=0,
    damage_float_new=0,damage_float_reach=0,damage_float_total=0,
    damage_float_coh=0,damage_float_isl=0,
    damage_fracture_complete=0,
    /* bounded exactly like the DE1.3 terminal batch: once the hydride is
       gone every facet on its surface loses its plus side at the same
       instant, and removing all of them in one topology transaction is a
       far larger redistribution than the same-load solve can absorb */
    damage_float_batch=DAMAGE_DE13_BATCH_MAX,
    damage_dangle_new=0,damage_dangle_weak=0,damage_dangle_total=0,
    damage_dangle_max=0,damage_stiff_probe=0,damage_delete_visc=1,damage_path_on=0,damage_path_retry=0,damage_path_desc=0,
    damage_path_nstep=20,damage_path_arm=3,damage_path_used=0,
    damage_path_att=0,damage_null_inc=0,damage_null_it=0,
    damage_null_cnt=0,damage_null_seed=987654321,damage_null_nit=4,
    damage_stab_maxdof=0,damage_stab_maxdead=0,*damage_stab_node=NULL,
    damage_deadsole_total=0,damage_deadall_total=0,damage_deadall_nodes=0,
    damage_fracture_link=0,damage_facetdel=0,
    damage_facetdel_new=0,damage_facetdel_total=0,damage_arc=0,damage_diss_step=1,damage_spc_neg=0,damage_stiff_nneg=0,
    damage_ls_trials=DAMAGE_LINESEARCH_MAX_TRIALS,
    damage_bare=0,damage_bare_rep=-1,damage_free_probe=0,
    damage_free_cnt=0,damage_free_rep=-1,*damage_free_nb=NULL,
    damage_free_worst=-1,damage_free_raw=0,damage_free_rawrep=-1,
    damage_dump_node=0,damage_dump_inc=1,damage_dump_n=0,
    damage_dump_nb=0,damage_dump_nu=0,damage_dump_alive=0,
    damage_dump_idx=0,damage_dump_np=0,damage_dump_hit=0,
    damage_fd_inc=0,damage_fd_it=1,damage_fd_ncol=0,damage_fd_el=-1,
    damage_fd_j=0,damage_fd_s=0,damage_fd_d=0,damage_fd_node=0,
    damage_fd_col=0,damage_fd_row=0,damage_fd_worst=-1,
    damage_fd_nbad=0,damage_fd_wd=0,damage_unsym_skip=0,
    damage_unsym_skiprep=-1,damage_unsym_adv=0,
    damage_unsym_advrep=-1,
    *damage_stiff_haz=NULL,
    damage_stiff_n1=0,damage_stiff_n2=0,damage_stiff_n3=0,
    damage_stiff_worst=-1,
    damage_topology_orphans=0,damage_indexe=0,
    damage_ls_bestused=0,damage_ls_legacy=0,damage_lsr=0,
    damage_de13_new=0,damage_de13_transaction=0,
    damage_slow_active=0,damage_slow_extended=0,damage_slow_nsoft=0,
    damage_slow_allow=0,damage_slow_estres=0,damage_slow_estcorr=0,
    damage_slow_esttotal=0,damage_slow_maxiters=0;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*epn=NULL,*cg=NULL,
    *cdn=NULL,*pslavsurfold=NULL,*fextload=NULL,
    *f=NULL,*fn=NULL,qa[4]={0.,0.,-1.,0.},qam[2]={0.,0.},dtheta,theta,
    err,ram[6]={0.,0.,0.,0.,0.,0.},*areaslav=NULL,
    *springarea=NULL,ram1[6]={0.,0.,0.,0.,0.,0.},
    ram2[6]={0.,0.,0.,0.,0.,0.},deltmx,ptime,smaxls,sminls,
    uam[2]={0.,0.},*vini=NULL,*ac=NULL,qa0,qau,ea,*straight=NULL,
    *t1act=NULL,qamold[2],*xbounact=NULL,*bc=NULL,
    *xforcact=NULL,*xloadact=NULL,*fext=NULL,*clearini=NULL,
    reltime,time,bet=0.,gam=0.,*aux2=NULL,dtime,*fini=NULL,
    *fextini=NULL,*veini=NULL,*accini=NULL,*xstateini=NULL,
    *ampli=NULL,scal1,*eei=NULL,*t1ini=NULL,pressureratio,
    *xbounini=NULL,dev,*xstiff=NULL,*stx=NULL,*stiini=NULL,
    *enern=NULL,*coefmpc=NULL,*aux=NULL,*xstaten=NULL,
    *coefmpcref=NULL,*enerini=NULL,*emn=NULL,alpham,betam,
    *tarea=NULL,*tenv=NULL,*erad=NULL,*fnr=NULL,*fni=NULL,
    *adview=NULL,*auview=NULL,*qfx=NULL,*cvini=NULL,*cv=NULL,
    *qfn=NULL,*co=NULL,*vold=NULL,*fenv=NULL,sigma=0.,
    *xbodyact=NULL,*cgr=NULL,dthetaref,dthetadamage,thetadamage,
    dthetarefdamage,theta_goal=0.,theta_local_start=0.,
    dtheta_restore=0.,dtheta_remaining, *vr=NULL,*vi=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*fmpc=NULL,*ener=NULL,
    *f_cm=NULL, *f_cs=NULL,*adc=NULL,*auc=NULL,*res=NULL,
    *xstate=NULL,*eenmax=NULL,*adrad=NULL,*aurad=NULL,*bcr=NULL,
    *xmastnor=NULL,*emeini=NULL,*tinc,*tper,*tmin,*tmax,*tincf,
    *doubleglob=NULL,*xnoels=NULL,*au=NULL,*resold=NULL,
    *damage_linesearch_step=NULL,
    *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,
    *x=NULL,*y=NULL,*z=NULL,*xo=NULL,sum1,sum2,flinesearch,
    *yo=NULL,*zo=NULL,*cdnr=NULL,*cdni=NULL,*fnext=NULL,*fnextini=NULL,
    allwk=0.,allwkini,energyini[4]={0.,0.,0.,0.},*cof=NULL,
    energyref,dtcont,dtvol,wavespeed[*nmat],emax,r_abs,
    enetoll,dampwk=0.,dampwkini=0.,temax,*tmp=NULL,energystartstep[4],
    sizemaxinc,*adblump=NULL,*adcpy=NULL,*aucpy=NULL,*rwork=NULL,
    *sol=NULL,*rgwk=NULL,tol,*sb=NULL,*sx=NULL,delcon,alea,
    *smscale=NULL,dtset,energym=0.,energymold=0.,*voldf=NULL,
    *coefmpcf=NULL,*xbounf=NULL,*xloadf=NULL,*xbounoldf=NULL,
    *xbounactf=NULL,*xloadoldf=NULL,*xloadactf=NULL,*auw=NULL,*volddof=NULL,
    *qb=NULL,*aloc=NULL,dtmin,*fric=NULL,*aubi=NULL,*auib=NULL,
    *fullgmatrix=NULL,*fullr=NULL,*alglob=NULL,*damn=NULL,*errn=NULL,
    *damdamageini=NULL,*damde1prev=NULL,*veolddamageini=NULL,
    *damage_tent_value=NULL,*damage_de13_trigger_value=NULL,*damagebase=NULL,
    *damage_damjac=NULL,damage_snap_ratio=0.,
    *damage_damvisc=NULL,*damage_damviscini=NULL,damage_visc_eta=0.,
    *damage_frel=NULL,*damage_ray_p=NULL,*damage_ray_res=NULL,
    *damage_ray_r0=NULL,*damage_bt_dam=NULL,*damage_bt_visc=NULL,
    *damage_bt_xs=NULL,damage_bt_r0=0.,damage_aba_a=0.25,
    damage_rescue_dtheta_last=0.,damage_rescue_dthetaref_last=0.,
    damage_reg_lambda=0.,
    damage_reg_lam[5]={1.e-2,1.e-1,1.e0,4.e0,1.6e1},
    damage_corr_lam=0.,damage_corr_lamstable=0.,damage_corr_dtref=0.,
    damage_corr_dtsum=0.,damage_corr_theta0=0.,damage_corr_minfrac=0.05,
    damage_dtheta_healthy=0.,
    damage_dth_ring[20]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                         0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
    damage_bt_growth=1.,damage_bt_floor=0.015625,
    damage_bt_ring[8]={0.,0.,0.,0.,0.,0.,0.,0.},
    *daba_res=NULL,*daba_v=NULL,*daba_stx=NULL,*daba_fn=NULL,
    *daba_f=NULL,*daba_dam=NULL,*daba_visc=NULL,*daba_xs=NULL,
    *daba_eme=NULL,*daba_stiff=NULL,*daba_qa=NULL,*daba_cam=NULL,
    damage_ray_growth=1.10,
    *damage_ct_ring=NULL,*damage_ct_w=NULL,*damage_ct_r0=NULL,
    *damage_ct_beps=NULL,*damage_ct_y=NULL,*damage_ct_z=NULL,
    *damage_ct_qh=NULL,
    *damage_ct_dam=NULL,
    *damage_ct_visc=NULL,*damage_ct_xs=NULL,*damage_ct_jac=NULL,
    damage_ct_dt[6]={0.,0.,0.,0.,0.,0.},
    damage_ct_dlam[6]={0.,0.,0.,0.,0.,0.},
    damage_ct_duinf[6]={0.,0.,0.,0.,0.,0.},
    damage_ct_m[3]={0.,0.,0.},damage_ct_g[3]={0.,0.,0.},
    damage_ct_dc[3]={0.,0.,0.},
    damage_ct_lam=0.,damage_ct_lamc=0.,damage_ct_ds=0.,damage_ct_ds0=0.,
    damage_ct_dsmin=0.,damage_ct_dsmax=0.,damage_ct_kappa=0.,
    damage_ct_clam=0.,damage_ct_tau=0.,damage_ct_tolc=0.,
    damage_ct_rhomin=1.e-4,damage_ct_clim=20.,damage_ct_ulim=20.,
    damage_ct_lamref=0.,damage_ct_duref=0.,damage_ct_taupr=0.,
    damage_ct_den=0.,damage_ct_rho=0.,damage_ct_eps=1.e-6,
    damage_ct_kaptol=1.5,damage_ct_qa[4]={0.,0.,0.,0.},
    damage_ct_cam[5]={0.,0.,0.,0.,0.},damage_ct_uam[2]={0.,0.},
    damage_ct_lamsnap=0.,damage_ct_cprev=0.,
    *damage_dl_r0=NULL,*damage_dl_d=NULL,*damage_dl_w=NULL,
    *damage_dl_pn=NULL,*damage_dl_res=NULL,*damage_dl_dam=NULL,
    *damage_dl_pm=NULL,*damage_dl_wm=NULL,damage_dl_npm2=0.,
    *damage_dl_visc=NULL,*damage_dl_xs=NULL,
    damage_dl_delta=0.,damage_dl_nb2=0.,damage_dl_nd2=0.,
    damage_dl_nw2=0.,damage_dl_tc=0.,damage_dl_npn2=0.,
    damage_dl_dtpn=0.,damage_dl_ident=0.,damage_dl_asym=0.,
    damage_dl_d0fac=1.,damage_dl_dmax=0.,damage_dl_phi0=0.,
    damage_release_qa=0.,damage_release_qam=0.,
    damage_release_dt=0.,
    *damage_diss_fhat=NULL,*damage_diss_uf=NULL,
    *damage_addiag=NULL,*damage_addiag0=NULL,
    damage_diss_p=0.,damage_diss_pprev=0.,damage_diss_lprev=0.,
    damage_diss_dg=0.,damage_diss_total=0.,damage_diss_target=0.,damage_ls_min=DAMAGE_LINESEARCH_MIN,damage_arc_lam=0.,damage_diss_engage_t=-1.,damage_arc_theta0=0.,damage_diss_lamnow=0.,
    damage_diss_scale=1.,damage_diss_dtheta=0.,
    damage_diss_lamcur=0.,damage_diss_lamold=0.,damage_diss_dgcur=0.,
    damage_diss_dgold=0.,damage_diss_g=0.,damage_diss_slope=0.,
    damage_diss_dlam=0.,damage_diss_kpp=0.,damage_diss_fr=0.,
    damage_diss_ff=0.,damage_diss_den=0.,
    damage_addmin=0.,damage_addrat=0.,damage_de13_batch_vmin=1.,
    damage_stiff_min=0.,damage_path_lam=0.,damage_path_dev=0.,
    damage_fd_h=1.e-7,damage_fd_dmax=0.,damage_fd_num=0.,
    damage_fd_asm=0.,damage_fd_amax=0.,damage_fd_emax=0.,
    *damage_fd_vsav=NULL,*damage_fd_fp=NULL,*damage_fd_fm=NULL,
    *damage_fd_ad=NULL,*damage_fd_au=NULL,
    damage_fd_wa=0.,damage_fd_wf=0.,
    damage_path_devmax=-1.,damage_path_ref=0.,
    damage_path_lamcom=0.,damage_path_drop=0.25,
    *damage_null_x=NULL,damage_null_nb=0.,damage_null_nx=0.,
    damage_null_amax=0.,damage_null_nn=0.,
    *damage_free_g=NULL,damage_free_gm=0.,damage_free_dv=0.,
    damage_dump_dv=0.,
    damage_de13_delete_d=DAMAGE_DE13_DELETE_D,
    damage_tent_step_time=0.,damage_tent_total_time=0.,damage_dmax=0.,
    damage_alphaevent=2.,damage_event_dtheta=0.,
    damage_event_raw=0.,damage_event_floor=0.,
    damage_de1_dmax=0.,damage_de1_maxdelta=0.,
    damage_de13_batch_dmax=0.,damage_slow_maxdd=0.,
    damage_slow_camprev1=1.e300,damage_slow_camprev2=1.e300,
    damage_slow_rratio=0.,damage_slow_cratio=0.,
    damage_reeq_uam_ref[2]={0.,0.},damage_reeq_uam_actual[2]={0.,0.},
    damage_reeq_uam_floor=0.,damage_reeq_uam_peak[2]={0.,0.},
    damage_linesearch_oldnorm=0.,damage_linesearch_fullnorm=0.,
    damage_linesearch_dampednorm=0.,damage_linesearch_maxdd=0.,
    damage_wall_theta=-1.,*damage_wall_def=NULL;
  ITG damage_wall_maskstep=0;
  ITG damage_spc_force=0,damage_spc_fnode=0,damage_spc_fcount=0;
  double damage_spc_fmax=0.;
  double damage_nl_ell=0.;
  ITG damage_nl_mode=0;
  double damage_qam_floor=0.,damage_qam_peak=0.;
  double damage_stab_alpha=0.,damage_deadsole_g=0.,damage_deadall_g=0.,
    damage_spc_g=0.;
  char *damage_stab_env=NULL,*damage_deadsole_env=NULL,
    *damage_deadall_env=NULL;

  /* ---- CCX_PATHFOLLOW: consistent dissipation path following ----------
     Everything below is inert unless CCX_PATHFOLLOW is set to a positive
     dissipation increment.  The state proper lives in pathfollow.c; these
     are only the handles the Newton loop needs. */

  /* ---- topology / deletion diagnostics (topodiag.c) -----------------
     Measurement only.  td_trace prints a content hash of the COMMITTED
     state at the top of every attempt and the sorted composition of every
     deletion batch, which is what turns "the same batch repeats from the
     same state" from an inference into a measurement.  td_from runs the
     connectivity/rank/residual report from that increment on. */

  ITG td_trace=0,td_from=0,*td_comp=NULL,*td_sort=NULL,td_armed=0;
  ITG td_nmode=0,td_qneq=0;
  double *td_qkeep=NULL;
  unsigned long long td_state=0ULL,td_batch=0ULL;
  topodiag_report td_rep;

  lsladder damage_lsl;

  char *pf_env=NULL;
  ITG pf_on=0,pf_engaged=0,pf_pending=0,pf_reason=0,pf_applied=0,
    pf_neqarm=0,pf_nstep=0,pf_icutbprev=0,pf_ncut=0;
  double pf_lam=0.,pf_lamprev=0.,pf_dlampred=0.,pf_dlamjump=0.,pf_tauv=0.,
    pf_dg=0.,pf_g=0.,pf_dlam=0.,pf_dgc=0.,pf_clip=0.05,pf_taucur=0.,
    pf_dtheta_eng=1.e-3,pf_pdu=0.,*pf_uf=NULL,*pf_rhs0=NULL,*pf_y=NULL,
    *pf_uref=NULL,*pf_p=NULL,*pf_r0=NULL,*pf_r1=NULL,*pf_q=NULL,
    *pf_sv=NULL,*pf_sxs=NULL,*pf_sdam=NULL,*pf_sf=NULL,*pf_sfn=NULL,
    *pf_sstx=NULL,*pf_sxb=NULL,*pf_sxst=NULL,*pf_lhs=NULL,*pf_rhs=NULL;
  double pf_sqa[4],pf_scam[5];
  ITG pf_lincheck=0,pf_codmode=0;
  double pf_dphi=0.,pf_cu=0.,*pf_cvec=NULL;

  /* ---- mixed-mode crack control (pf_codmode==2) ---------------------
     The control functional is rebuilt from the COMMITTED state at the top
     of every attempt and the constraint is written incrementally, so a
     redefinition between increments carries nothing over.  See
     crackcontrol.c. */

  ITG pf_ccmode=0,pf_ccnw=0,pf_ccengage=0,pf_ccarmed=0,pf_ncutfloor=0;
  double pf_dphicur=0.,pf_lam0it=0.,pf_dlamit=0.,pf_gacc=0.,pf_phiacc=0.,
    pf_fhcos=0.,pf_fhrat=0.,pf_eps=1.e-5,pf_ccgrow=1.1;
  crackcontrol_census pf_cs;
	 
  FILE *f1,*fdamage=NULL;

#ifdef SGI
  ITG token;
#endif
  
  /* declarations for mortar contact */

  ITG *nslavspc=NULL,*islavspc=NULL,*nslavmpc=NULL,*islavmpc=NULL,
    *nmastspc=NULL,*imastspc=NULL,*nmastmpc=NULL,*imastmpc=NULL,
    *islavactdof=NULL,*islavactini=NULL,*islavtie=NULL,
    *irowt=NULL,*jqt=NULL,*irowtinv=NULL,*jqtinv=NULL,
    *irowb=NULL,*jqb=NULL,*irowd=NULL,*jqd=NULL,*irowdtil=NULL,*jqdtil=NULL,
    *irowbtil=NULL,*jqbtil=NULL,*irowbhelp=NULL,*jqbhelp=NULL,
    *islavnodeinv=NULL,*islavquadel=NULL,*irowc2=NULL,*jqc2=NULL,nzsc2,
    *icolc2=NULL,
    *jqbd=NULL,*irowbd=NULL,*jqbdtil=NULL,*irowbdtil=NULL,*jqbdtil2=NULL,
    *irowbdtil2=NULL,
    *jqdd=NULL,*irowdd=NULL,*jqddtil=NULL,*irowddtil=NULL,*jqddtil2=NULL,
    *irowddtil2=NULL,
    *jqddinv=NULL,*irowddinv=NULL,*jqtemp=NULL,*irowtemp=NULL,*icoltemp=NULL,
    nzstemp[3];
  
  double *bp=NULL,*gap=NULL,*slavnor=NULL,
    *slavtan=NULL,*cdisp=NULL,*cstress=NULL,*cfs=NULL,*cfm=NULL,*cfsini=NULL,
    *cfstil=NULL,*bpini=NULL,
    *cstressini=NULL,*pslavdual=NULL,*aut=NULL,
    *autinv=NULL,*Bd=NULL,*Bdhelp=NULL,
    *Dd=NULL,*Ddtil=NULL,*Bdtil=NULL,*auc2=NULL,*adc2=NULL,*aubd=NULL,
    *audd=NULL,*auddtil=NULL,*auddtil2=NULL,*auddinv=NULL,*bhat=NULL,
    *aubdtil=NULL,*aubdtil2=NULL;

  /* end of declarations for mortar contact */

  icol=*icolp;irow=*irowp;co=*cop;vold=*voldp;
  ipkon=*ipkonp;lakon=*lakonp;kon=*konp;ielorien=*ielorienp;
  ielmat=*ielmatp;ener=*enerp;xstate=*xstatep;
  
  ipompc=*ipompcp;labmpc=*labmpcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  fmpc=*fmpcp;nodempc=*nodempcp;coefmpc=*coefmpcp;nelemload=*nelemloadp;
  iamload=*iamloadp;sideload=*sideloadp;

  islavsurf=*islavsurfp;pslavsurf=*pslavsurfp;clearini=*clearinip;

  /* determining whether a node belongs to at least one element
     (needed in resultsforc.c) */
  
  NNEW(iponoel,ITG,*nk);
  FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));

  if(filab[4]!=' ') ne1d2d=1;

  num_cpus=0;
  sys_cpus=0;
  
  /* explicit user declaration prevails */
  
  envsys=getenv("NUMBER_OF_CPUS");
  if(envsys){
    sys_cpus=atoi(envsys);
    if(sys_cpus<0) sys_cpus=0;
  }
  
  /* automatic detection of available number of processors */
  
  if(sys_cpus==0){
    sys_cpus=getSystemCPUs();
    if(sys_cpus<1) sys_cpus=1;
  }
  
  /* else global declaration, if any, applies */
  
  env = getenv("OMP_NUM_THREADS");
  if(num_cpus==0){
    if(env)
      num_cpus=atoi(env);
    if(num_cpus<1) {
      num_cpus=1;
    }else if(num_cpus>sys_cpus){
      num_cpus=sys_cpus;
    }
  }
  
  // MPADD: initialize rmin to the tolerance
  enetoll=0.02;
  r_abs=0.0;
  emax=0.0;
  // MPADD end

  delcon=ctrl[53];alea=ctrl[54];

  tinc=&timepar[0];
  tper=&timepar[1];
  tmin=&timepar[2];
  tmax=&timepar[3];
  tincf=&timepar[4];

  if(*ithermal==4){
    uncoupled=1;
    *ithermal=3;
  }else{
    uncoupled=0;
  }

  /* for massless explicit dynamics any other "nonlingeo" step in the same
     calculation (e.g. a static step or an implicit dynamics step) 
     is performed with node-to-face contact */
  
  if(*mortar!=1){
    if(*nintpoint!=0){
      maxprevcontel=0;
      *nintpoint=0;
      SFREE(pslavsurf);SFREE(clearini);
    }else{
      maxprevcontel=*nslavs;
    }
  }else if(*mortar==1){
    maxprevcontel=*nintpoint;
    if(*nstate_!=0){
      if(maxprevcontel!=0){
	MNEW(islavsurfold,ITG,2**ifacecount+2);
	MNEW(pslavsurfold,double,3**nintpoint);
	isiz=2**ifacecount+2;cpyparitg(islavsurfold,islavsurf,&isiz,&num_cpus);
	isiz=3**nintpoint;cpypardou(pslavsurfold,pslavsurf,&isiz,&num_cpus);
      }
    }
    nslavs_prev_step=*nslavs;
  }

  /* turbulence model 
     iturbulent==0: laminar
     iturbulent==1: k-epsilon
     iturbulent==2: q-omega
     iturbulent==3: BSL
     iturbulent==4: SST */
  
  iturbulent=(ITG)physcon[8];
  
  for(k=0;k<3;k++){
    strcpy1(&jobnamef[k*132],&jobnamec[k*132],132);
  }
  
  qa0=ctrl[20];qau=ctrl[21];ea=ctrl[23];deltmx=ctrl[26];
  i0ref=ctrl[0];irref=ctrl[1];icref=ctrl[3];

  sminls=ctrl[28];smaxls=ctrl[29];
  
  memmpc_=mpcinfo[0];mpcfree=mpcinfo[1];icascade=mpcinfo[2];
  maxlenmpc=mpcinfo[3];

  alpham=xmodal[0];
  betam=xmodal[1];

  /* check whether, for a dynamic calculation, damping is involved */
  
  if(*nmethod==4){
    if(*iexpl<=1){
	  
      /* implicit dynamics */
	  
      if((fabs(alpham)>1.e-30)||(fabs(betam)>1.e-30)){
	idamping=1;idampingwithoutcontact=1;
      }else{
	for(i=0;i<*ne;i++){
	  if(ipkon[i]<0) continue;
	  if(strcmp1(&lakon[i*8],"ED")==0){
	    idamping=1;idampingwithoutcontact=1;break;
	  }
	}
      }
    }else{
	  
      /* explicit dynamics */
	  
      if((fabs(alpham)>1.e-30)||((fabs(betam)>1.e-30)&&(*mortar==-1))){
	idamping=1;idampingwithoutcontact=1;
      }
      if((fabs(betam)>1.e-30)&&(*mortar!=-1)){
	printf
	  (" *ERROR in nonlingeo: in explicit dynamic calculations\n");
	printf
	  ("         without massless contact the damping is only\n");
	printf
	  ("         allowed to be mass proportional: the coefficient beta\n");
	printf("         of the stiffness proportional term must be zero\n");
	FORTRAN(stop,());
      }
    }
  }
  
  /* check whether a sensitivity step may follow (whether design variables
     were defined */

  for(i=0;i<*ntie;i++){
    if(strcmp1(&tieset[i*243+80],"D")==0){
      isensitivity=1;
      NNEW(adcpy,double,neq[1]);
      /* no asymmetric matrices allowed for sensitivity */
      NNEW(aucpy,double,nzs[1]);
      break;
    }
  }

  if((icascade==2)&&(*iexpl>1)){
    printf
      (" *ERROR in nonlingeo: linear and nonlinear MPC's depend on each other\n");
    printf("        This is not allowed in a explicit dynamic calculation\n");
    FORTRAN(stop,());
  }
      
  /* determining the global values to be used as boundary conditions
     for a submodel */

  ITG irefine=0;
  getglobalresults(&jobnamec[396],&integerglob,&doubleglob,nboun,iamboun,xboun,
		   nload,sideload,iamload,&iglob,nforc,iamforc,xforc,
		   ithermal,nk,t1,iamt1,&sigma,&irefine);
  
  if(iglob<0){
    printf(" *ERROR in nonlingeo: a submodel calculation for which\n");
    printf("        the global model results from a *FREQUENCY\n");
    printf("        calculation must be geometrically linear\n");
    FORTRAN(stop,());
  }

  /* reading temperatures from frd-file */
  
  if((itempuser[0]==2)&&(itempuser[1]!=itempuser[2])) {
    utempread(t1,&itempuser[2],jobnamec);
  }      
  
  /* invert nactdof */
  
  /*  NNEW(nactdofinv,ITG,mt**nk);
  MNEW(nodorig,ITG,*nk);
  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			 ipkon,lakon,kon,ne));
			 SFREE(nodorig);*/
  
  /* allocating a field for the stiffness matrix */
  
  NNEW(xstiff,double,(long long)27*mi[0]**ne);
  
  /* allocating force fields */
  
  NNEW(f,double,neq[1]);
  NNEW(fext,double,neq[1]);
  
  NNEW(b,double,neq[1]);
  NNEW(vini,double,mt**nk);
  
  NNEW(aux,double,7*maxlenmpc);
  NNEW(iaux,ITG,2*maxlenmpc);
  
  /* allocating fields for the actual external loading */
  
  NNEW(xbounact,double,*nboun);
  NNEW(xbounini,double,*nboun);
  for(k=0;k<*nboun;++k){
    xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){
    xbodyact[k]=xbody[k];}
  
  /* assigning the body forces to the elements */ 
  
  if(*nbody>0){

    /* check whether there the previous step was in the relative
       system and a change to the absolute system was requested */
      
    if((*nmethod==4)&&(alpha[1]>1.)){
      NNEW(itreated,ITG,*nk);
      FORTRAN(velinireltoabs,(ibody,xbody,cbody,nbody,set,
			      istartset,iendset,ialset,nset,veold,mi,
			      ipkon,kon,lakon,co,itreated));
      SFREE(itreated);
    }
    if(*inewton==1){NNEW(cgr,double,4**ne);}
  }
  
  /* for mechanical calculations: updating boundary conditions
     calculated in a previous thermal step */
  
  if(*ithermal<2) FORTRAN(gasmechbc,(vold,nload,sideload,
				     nelemload,xload,mi));
  
  /* for thermal calculations: forced convection and cavity
     radiation*/
  
  if(*ithermal>1){
    NNEW(itg,ITG,*nload+3**nflow);
    NNEW(ieg,ITG,*nflow);
    /* max 6 triangles per face, 4 entries per triangle */
    NNEW(kontri,ITG,24**nload);
    NNEW(nloadtr,ITG,*nload);
    NNEW(nacteq,ITG,4**nk);
    NNEW(nactdog,ITG,4**nk);
    NNEW(v,double,mt**nk);
    FORTRAN(envtemp,(itg,ieg,&ntg,&ntr,sideload,nelemload,
		     ipkon,kon,lakon,ielmat,ne,nload,
		     kontri,&ntri,nloadtr,nflow,ndirboun,nactdog,
		     nodeboun,nacteq,nboun,ielprop,prop,&nteq,
		     v,network,physcon,shcon,ntmat_,co,
		     vold,set,nshcon,rhcon,nrhcon,mi,nmpc,nodempc,
		     ipompc,labmpc,ikboun,&nasym,ttime,&time,
		     iaxial));
    SFREE(v);
      
    if((*mcs>0)&&(ntr>0)){
      NNEW(inocs,ITG,*nk);
      radcyc(nk,kon,ipkon,lakon,ne,cs,mcs,nkon,ialset,istartset,
	     iendset,&kontri,&ntri,&co,&vold,&ntrit,inocs,mi);
    }
    else{ntrit=ntri;}
      
    nzsrad=100*ntr;
    NNEW(mast1rad,ITG,nzsrad);
    NNEW(irowrad,ITG,nzsrad);
    NNEW(icolrad,ITG,ntr);
    NNEW(jqrad,ITG,ntr+1);
    NNEW(ipointerrad,ITG,ntr);
      
    if(ntr>0){
      mastructrad(&ntr,nloadtr,sideload,ipointerrad,
		  &mast1rad,&irowrad,&nzsrad,
		  jqrad,icolrad);
    }
      
    /* determine the network elements belonging to a given node (for usage
       in user subroutine film */

    if((*network>0)||(ntg>0)){
      NNEW(iponoeln,ITG,*nk);
      NNEW(inoeln,ITG,2**nkon);
      if(*network>0){
	FORTRAN(networkelementpernode,(iponoeln,inoeln,lakon,ipkon,kon,
				       &inoelnsize,nflow,ieg,ne,network));
	FORTRAN(checkforhomnet,(ieg,nflow,lakon,ipkon,kon,itg,&ntg,
				iponoeln,inoeln));
      }
      RENEW(inoeln,ITG,2*inoelnsize);
    }

    SFREE(ipointerrad);SFREE(mast1rad);
    RENEW(irowrad,ITG,nzsrad);
      
    RENEW(itg,ITG,ntg);
    NNEW(ineighe,ITG,ntg);
    RENEW(kontri,ITG,4*ntrit);
    RENEW(nloadtr,ITG,ntr);
      
    NNEW(adview,double,ntr);
    NNEW(auview,double,2*nzsrad);
    NNEW(tarea,double,ntr);
    NNEW(tenv,double,ntr);
    NNEW(fenv,double,ntr);
    NNEW(erad,double,ntr);
      
    NNEW(ac,double,nteq*nteq);
    NNEW(bc,double,nteq);
    NNEW(ipiv,ITG,nteq);
    NNEW(adrad,double,ntr);
    NNEW(aurad,double,2*nzsrad);
    NNEW(bcr,double,ntr);
    NNEW(ipivr,ITG,ntr);
  }
  
  /* check for fluid elements
     check for strain-less elements */
  
  NNEW(nactdoh,ITG,*ne);
  NNEW(nactdohinv,ITG,*ne);
  *nef=0;
  for(i=0;i<*ne;++i){
    if(ipkon[i]<0) continue;
    if(strcmp1(&lakon[8*i],"F")==0){
      icfd=1;nactdohinv[*nef]=i+1;++*nef;nactdoh[i]=*nef;}
    if(istrainfree==0){
      if(ielmat[i]<0){istrainfree=1;}
    }
  }

  if(icfd==1){
    if(iturbulent>=20){

      /* CBS method for shallow water equations */
      
      iturbulent=iturbulent-20;
      ifreesurface=1;
      icfd=2;
    }else if(iturbulent>=10){

      /* CBS method for all other applications */
      
      iturbulent=iturbulent-10;
      icfd=2;
    }
  }
  
  if(icfd==1){
  }else if(icfd==2){
      SFREE(nactdoh);SFREE(nactdohinv);

      /* rearranging the fluid nodes and elements such that
	 no gaps occur */

      NNEW(ipkonf,ITG,*nef);
      NNEW(lakonf,char,8**nef);
      NNEW(ielmatf,ITG,mi[2]**nef);
      if(*norien>0) NNEW(ielorienf,ITG,mi[2]**nef);
      NNEW(nelold,ITG,*nef);
      NNEW(nelnew,ITG,*ne);
      NNEW(cof,double,3**nk);
      NNEW(voldf,double,mt**nk);
      NNEW(nkold,ITG,*nk);
      NNEW(nknew,ITG,*nk);
      NNEW(inotrf,ITG,2**nk);
      NNEW(konf,ITG,*nkon);
      NNEW(ipompcf,ITG,*nmpc);
      NNEW(ikmpcf,ITG,*nmpc);
      NNEW(ilmpcf,ITG,*nmpc);
      NNEW(nodempcf,ITG,3*memmpc_);
      NNEW(coefmpcf,double,memmpc_);
      NNEW(nodebounf,ITG,*nboun);
      NNEW(ndirbounf,ITG,*nboun);
      NNEW(ikbounf,ITG,*nboun);
      NNEW(ilbounf,ITG,*nboun);
      if(*nam>0) NNEW(iambounf,ITG,*nboun);
      NNEW(xbounf,double,*nboun);
      NNEW(xbounoldf,double,*nboun);
      NNEW(xbounactf,double,*nboun);
      NNEW(nelemloadf,ITG,2**nload);
      if(*nam>0) NNEW(iamloadf,ITG,2**nload);
      NNEW(xloadf,double,2**nload);
      NNEW(xloadoldf,double,2**nload);
      NNEW(xloadactf,double,2**nload);
      NNEW(sideloadf,char,20**nload);
      if(*nbody>0) NNEW(ipobodyf,ITG,2*(*ifreebody-1));

      FORTRAN(rearrangecfd,(ne,ipkon,lakon,ielmat,ielorien,norien,nef,ipkonf,
			    lakonf,ielmatf,ielorienf,mi,nelold,nelnew,nkold,
			    nknew,nk,&nkf,konf,&nkonf,nmpc,ipompc,nodempc,
			    coefmpc,&memmpc_,&nmpcf,ipompcf,nodempcf,coefmpcf,
			    &memmpcf,nboun,nodeboun,ndirboun,xboun,&nbounf,
			    nodebounf,ndirbounf,xbounf,nload,nelemload,
			    sideload,xload,&nloadf,nelemloadf,sideloadf,xloadf,
			    ipobody,ipobodyf,kon,&nkftot,co,cof,vold,voldf,
			    ikbounf,ilbounf,ikmpcf,ilmpcf,iambounf,iamloadf,
			    iamboun,iamload,xbounold,xbounoldf,xbounact,
			    xbounactf,xloadold,xloadoldf,xloadact,xloadactf,
			    inotr,inotrf,nam,ntrans,nbody));

      /* call rearrangecfd */
      
      RENEW(nkold,ITG,nkftot);
      RENEW(cof,double,3*nkftot);
      RENEW(voldf,double,mt*nkftot);
      NNEW(inotrf,ITG,2*nkftot);
      RENEW(konf,ITG,nkonf);
      RENEW(ipompcf,ITG,nmpcf);
      RENEW(ikmpcf,ITG,nmpcf);
      RENEW(ilmpcf,ITG,nmpcf);
      RENEW(nodempcf,ITG,3*memmpcf);
      RENEW(coefmpcf,double,memmpcf);
      RENEW(nodebounf,ITG,nbounf);
      RENEW(ndirbounf,ITG,nbounf);
      RENEW(ikbounf,ITG,nbounf);
      RENEW(ilbounf,ITG,nbounf);
      if(*nam>0) RENEW(iambounf,ITG,nbounf);
      RENEW(xbounf,double,nbounf);
      RENEW(xbounoldf,double,nbounf);
      RENEW(xbounactf,double,nbounf);
      RENEW(nelemloadf,ITG,2*nloadf);
      if(*nam>0) NNEW(iamloadf,ITG,2*nloadf);
      RENEW(xloadf,double,2*nloadf);
      RENEW(xloadoldf,double,2*nloadf);
      RENEW(xloadactf,double,2*nloadf);
      RENEW(sideloadf,char,20*nloadf);
      
      /* calculating topological properties for CFD */
      
      NNEW(sideface,char,6**nef);
      NNEW(nelemface,ITG,6**nef);
      NNEW(ipface,ITG,*nef);
      NNEW(ipoface,ITG,nkf);
      NNEW(nodface,ITG,5*6**nef);
      NNEW(ifreestream,ITG,nkf);
      NNEW(isolidsurf,ITG,nkf);
      NNEW(neighsolidsurf,ITG,nkf);
      NNEW(iponoelf,ITG,nkf);
      NNEW(inoelf,ITG,2*8**nef);
      NNEW(inomat,ITG,nkftot);
      FORTRAN(topocfdfem,(nelemface,sideface,&nface,ipoface,nodface,nef,ipkonf,
			  konf,lakonf,&nkf,isolidsurf,&nsolidsurf,ifreestream,
			  &nfreestream,neighsolidsurf,iponoelf,inoelf,
			  &inoelfree,cof,set,istartset,iendset,ialset,nset,
			  &iturbulent,inomat,ielmatf,ipface,nknew));
      RENEW(sideface,char,nface);
      RENEW(nelemface,ITG,nface);
      SFREE(ipoface);SFREE(nodface);
      RENEW(ifreestream,ITG,nfreestream);
      RENEW(isolidsurf,ITG,nsolidsurf);
      RENEW(neighsolidsurf,ITG,nsolidsurf);
      RENEW(inoelf,ITG,2*inoelfree);
      if(*ithermal==1){
	NNEW(qfx,double,3*mi[0]**ne);}
  }else{
    SFREE(nactdoh);SFREE(nactdohinv);
  }
  
  if(*ithermal>1){
    NNEW(qfx,double,3*mi[0]**ne);}
  
  /* contact conditions */
  
  inicont(nk,&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,&itietri,
	  lakon,ipkon,kon,&koncont,nslavs,tietol,&ismallsliding,&itiefac,
          &islavsurf,&islavnode,&imastnode,&nslavnode,&nmastnode,
          mortar,&imastop,nkon,&iponoels,&inoels,&ipe,&ime,ne,ifacecount,
	  iperturb,ikboun,nboun,co,istep,&xnoels);
  
  if(ncont!=0){
      
    NNEW(cg,double,3*ncont);
    NNEW(straight,double,16*ncont);
	  
    /* 11 instead of 10: last position is reserved for the
       local contact spring element number; needed as
       pointer into springarea */
      
    if(*mortar<=0){
      RENEW(kon,ITG,*nkon+11**nslavs);
      NNEW(springarea,double,2**nslavs);
      if((*nener==1)&&((maxprevcontel==0)&&(*nslavs!=0))){
	RENEW(ener,double,2*mi[0]*(*ne+*nslavs));
	DOUMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+*nslavs),0.);

	/* setting the entries for the friction contact energy to zero */

	/*	for(k=mi[0]*(2**ne+*nslavs);k<mi[0]*(*ne+*nslavs)*2;k++){
		ener[k]=0.;}*/
      }
      RENEW(ipkon,ITG,*ne+*nslavs);
      RENEW(lakon,char,8*(*ne+*nslavs));
	  
      if(*norien>0){
	RENEW(ielorien,ITG,mi[2]*(*ne+*nslavs));
	for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++){
	  ielorien[k]=0;}
      }

      RENEW(ielmat,ITG,mi[2]*(*ne+*nslavs));
      for(k=mi[2]**ne;k<mi[2]*(*ne+*nslavs);k++){
	ielmat[k]=1;}

      if((maxprevcontel==0)&&(*nslavs!=0)){
	RENEW(xstate,double,*nstate_*mi[0]*(*ne+*nslavs));
	for(k=*nstate_*mi[0]**ne;k<*nstate_*mi[0]*(*ne+*nslavs);k++){
	  xstate[k]=0.;
	}
      }
      maxprevcontel=*nslavs;

      NNEW(areaslav,double,*ifacecount);
    }else if(*mortar==1){
      NNEW(islavact,ITG,nslavnode[*ntie]);
      if((*istep==1)||(nslavs_prev_step==0))
	NNEW(clearini,double,3*9**ifacecount);

      /* check whether at least one contact definition involves true contact
	 and not just tied contact */

      FORTRAN(checktruecontact,(ntie,tieset,tietol,elcon,&itruecontact,
				ncmat_,ntmat_));
    }else if(*mortar>1){
      ismallsliding=1;
      NNEW(slavnor,double,3**nslavs);
      NNEW(slavtan,double,6**nslavs);
      inimortar(&ener,mi,ne,nslavs,nk,nener,&ipkon,&lakon,&kon,nkon,
		&maxprevcontel,&xstate,nstate_,&islavtie,&bp,&islavact,
		&gap,&cdisp,&cstress,&cfs,
		&bpini,&islavactini,&cstressini,ntie,
		tieset,nslavnode,islavnode,&islavnodeinv,&islavquadel,
		&pslavdual,&aut,&irowt,&jqt,&autinv,
		&irowtinv,&jqtinv,&Bd,&irowb,&jqb,&Bdhelp,&irowbhelp,
		&jqbhelp,&Dd,&irowd,&jqd,&Ddtil,&irowdtil,&jqdtil,&Bdtil,
		&irowbtil,&jqbtil,itiefac,islavsurf,nboun,
		nmpc,&nslavspc,&islavspc,&nslavmpc,&islavmpc,
		&nmastspc,&imastspc,&nmastmpc,&imastmpc,
		imastnode,nmastnode,&nasym,mortar,&ielmat,&ielorien,norien,
		ipompc,nodempc,ikboun,ilboun,ikmpc,ilmpc,jobnamef,set,
		co,vold,nset,&nslavquadel);
    }
    NNEW(xmastnor,double,3*nmastnode[*ntie]);
  }
  
  if(icascade==2){
    memmpcref_=memmpc_;mpcfreeref=mpcfree;maxlenmpcref=maxlenmpc;
    NNEW(nodempcref,ITG,3*memmpc_);
    for(k=0;k<3*memmpc_;k++){
      nodempcref[k]=nodempc[k];}
    NNEW(coefmpcref,double,memmpc_);
    for(k=0;k<memmpc_;k++){
      coefmpcref[k]=coefmpc[k];}
  }
  
  if((*ithermal==1)||(*ithermal>=3)){
    NNEW(t1ini,double,*nk);
    NNEW(t1act,double,*nk);
    for(k=0;k<*nk;++k){
      t1act[k]=t1old[k];}
  }
  
  /* allocating a field for the instantaneous amplitude */
  
  NNEW(ampli,double,*nam);
  
  /* fini is also needed in static calculations if iforbou=1
     to get correct values of f after a divergent increment */

  NNEW(fini,double,neq[1]);
  
  /* allocating fields for nonlinear dynamics */
  
  if(*nmethod==4){
    mass[0]=1;
    mass[1]=1;
    NNEW(aux2,double,neq[1]);
    NNEW(fextini,double,neq[1]);
    NNEW(fnext,double,mt**nk);
    NNEW(fnextini,double,mt**nk);
    NNEW(veini,double,mt**nk);
    NNEW(accini,double,mt**nk);
    NNEW(adb,double,neq[1]);
    NNEW(aub,double,nzs[1]);
    NNEW(cvini,double,neq[1]);
    NNEW(cv,double,neq[1]);
  }

  //    if((*nstate_!=0)&&((*mortar!=1)||(ncont==0))){
  if((*nstate_!=0)&&(*mortar!=1)){
    NNEW(xstateini,double,*nstate_*mi[0]*(*ne+*nslavs));
    isiz=*nstate_*mi[0]*(*ne+*nslavs);cpypardou(xstateini,xstate,&isiz,&num_cpus);
    //    FORTRAN(stop,());
  }

  /* next lines: change on 8th of July 2023: initial state values
     for dynamic plastic calculations */
  
  if((*nstate_!=0)&&(*mortar==1)){
    NNEW(xstateini,double,*nstate_*mi[0]**ne);
    isiz=*nstate_*mi[0]**ne;cpypardou(xstateini,xstate,&isiz,&num_cpus);
  }
  
  NNEW(eei,double,6*mi[0]**ne);
  NNEW(stiini,double,6*mi[0]**ne);
  NNEW(emeini,double,6*mi[0]**ne);
  
  if(*nener==1){
    if((*mortar!=1)||(ncont==0)){
      NNEW(enerini,double,2*mi[0]*(*ne+*nslavs));
    }else{
      NNEW(enerini,double,2*mi[0]*(*ne+*nintpoint));
    }
    isiz=2*mi[0]**ne;cpypardou(enerini,ener,&isiz,&num_cpus);
  }
  
  qa[0]=qaold[0];
  qa[1]=qaold[1];
  
  /* normalizing the time */
  
  FORTRAN(checktime,(itpamp,namta,tinc,ttime,amta,tmin,inext,&itp,istep,tper));
  dtheta=(*tinc)/(*tper);

  /* taking care of a small increment at the end of the step
     for face-to-face penalty contact */

  dthetaref=dtheta;
  if((dtheta<=1.e-6)&&(*iexpl<=1)){
    printf("\n *ERROR in nonlingeo\n");
    printf(" increment size smaller than one millionth of step size\n");
    printf(" increase increment size\n\n");
  }
  /* [SWITCHES] State the configuration of this run before anything acts on
     it, and name any CCX_* variable that is set but not read.  A flag that
     was never read is the difference between an A/B that is wrong and one
     that is uninformative, and the second kind costs a whole run to notice.
     Reports only; reads nothing, decides nothing. */
  damswitch_report();

  /* [DAMAGE TMIN] statics.f:234-247 silently raises the deck's minimum
     increment to min(tinc,1e-6*tper) under automatic incrementation.  With
     tinc=1e-3 and tper=1 every deck in this project has therefore been run
     with tmin=1e-6, not the 1e-9 it asks for, and "increment size smaller
     than minimum" has been a CalculiX floor rather than a physical limit.
     Opt-in, physical units, applied before the normalisation below. */

  if((damage_de13_env=getenv("CCX_DAMAGE_TMIN"))!=NULL){
    double tmnew=atof(damage_de13_env);
    if(tmnew>0.){
      printf("[DAMAGE TMIN] minimum increment overridden: %.12e -> %.12e "
             "(physical units).  The statics.f 1e-6*tper floor is bypassed; "
             "the stock cutback machinery is otherwise untouched.%s",
             *tmin,tmnew,"\n");
      fflush(stdout);
      *tmin=tmnew;
    }
  }

  *tmin=*tmin/(*tper);
  *tmax=*tmax/(*tper);
  theta=0.;
  
  /* calculating an initial flux norm */
  
  if(*ithermal!=2){
    if(qau>1.e-10){qam[0]=qau;}
    else if(qa0>1.e-10){qam[0]=qa0;}
    else if(qa[0]>1.e-10){qam[0]=qa[0];}
    else {qam[0]=1.e-2;}
  }
  if(*ithermal>1){
    if(qau>1.e-10){
      qam[1]=qau;}
    else if(qa0>1.e-10){
      qam[1]=qa0;}
    else if(qa[1]>1.e-10){
      qam[1]=qa[1];}
    else {qam[1]=1.e-2;}
  }
  
  /* storing the element and topology information before introducing 
     contact elements */
  
  ne0=*ne;nkon0=*nkon;neold=*ne;

  /* [LOADPATH] Arm the specimen-level judgement: is there still a load
     path between the grips at all.

     This sits OUTSIDE the `if(*ndmat_>0)` block below on purpose.  The
     block below is where every damage switch is read, and putting the
     judgement there would have repeated the defect it exists to remove:
     whether the specimen is still in one piece has nothing to do with
     whether the deck happens to declare a bulk damage material.  Measured:
     test/pathfollow/close.inp is 24 bulk elements in two halves joined by
     two cohesive facets and nothing else, both facets reach the failure
     flag, and with the arming inside that block it was never judged at all.

     Armed when the deck can actually come apart - it has a damage material
     or a cohesive element - and not otherwise, so a stock elastic run pays
     nothing. */
  {
    ITG lp_nu=0,lp_i;
    for(lp_i=0;lp_i<*ne;lp_i++) if(lakon[8*lp_i]=='U') lp_nu++;

    if((damage_de13_env=getenv("CCX_FRACTURE_LINK"))!=NULL){
      if((strcmp(damage_de13_env,"FACE")==0)||
         (strcmp(damage_de13_env,"face")==0)) damage_fracture_link=1;
    }
    if((damage_de13_env=getenv("CCX_FRACTURE_PAST_SEVERANCE"))!=NULL){
      damage_lp_stop=(strcmp(damage_de13_env,"0")==0)?1:0;
    }
    damage_fracture_env=getenv("CCX_FRACTURE_TERMINATION");
    if(damage_fracture_env!=NULL){
      damage_fracture_seta=strdup(damage_fracture_env);
      damage_fracture_setb=strchr(damage_fracture_seta,':');
      if(damage_fracture_setb!=NULL){
        *damage_fracture_setb=0;
        damage_fracture_setb++;
      }else{
        free(damage_fracture_seta);
        damage_fracture_seta=NULL;
        printf("[LOADPATH] *WARNING: expected "
               "CCX_FRACTURE_TERMINATION=SETA:SETB; ignored, and the "
               "endpoints come from the deck instead\n");
      }
    }

    if((*ndmat_>0)||(lp_nu>0)){
      /* nonlingeo() is entered once per STEP.  Initialise on the first
         entry only and merely release the endpoint lists on later ones:
         the severance latch has to survive a step boundary, or a specimen
         that came apart in step 1 would be judged intact again at the top
         of step 2.  close.inp is exactly that deck - it severs in tension
         at the end of step 1 and step 2 exists to close the crack. */
      if(damage_lp_armed_once==0){
        loadpath_init(&damage_lp);
        damage_lp_armed_once=1;
        damage_lp_testok=(loadpath_selftest()==0);
        convstate_selftest();
        /* [STALLSTATE] arm the self-arming diagnostics.  Refuses on a failed
           self test, like every other judgement here. */
        stallstate_init(&damage_stall);
        damage_stall_ok=(stallstate_selftest()==0);
        if(!damage_stall_ok){
          printf("[STALLSTATE] *ERROR: self test failed; the diagnostics will "
                 "NOT arm themselves and CCX_DAMAGE_WALL_THETA stays the only "
                 "way in.%s","\n");
        }else{
          printf("[STALLSTATE] armed: if dtheta falls below %.3g of the median "
                 "of this deck's own last %d accepted steps for %" ITGFORMAT
                 " consecutive increments, the wall diagnostics arm "
                 "themselves - no load factor has to be known in advance.%s",
                 damage_stall.ratio,STALL_RING,damage_stall.nconfirm,"\n");
        }
        fflush(stdout);
      }else{
        loadpath_free(&damage_lp);
      }
      if(!damage_lp_testok){
        printf("[LOADPATH] *ERROR: the specimen-level load-path self test "
               "failed; this run will make NO statement about whether the "
               "specimen stayed in one piece, rather than judging it by a "
               "rule that is not the one that was tested.\n");
        damage_lp_ready=0;
      }else{
        damage_lp_ready=loadpath_arm(&damage_lp,*nk,set,nset,istartset,
                                     iendset,ialset,
                                     damage_fracture_seta,
                                     damage_fracture_setb,
                                     nodeboun,ndirboun,xboun,*nboun,
                                     damage_fracture_link);
        if(damage_lp_ready){
          /* the verdict must survive a run that dies at a wall */
          loadpath_report_at_exit(&damage_lp);
          printf("[LOADPATH] grip-to-grip connectivity is checked every "
                 "converged increment; endpoints %s.  A cohesive facet "
                 "whose every integration point has failed does not count "
                 "as a load path.  On severance the run %s.\n",
                 damage_lp.origin,
                 damage_lp_stop?"stops":
                 "CONTINUES (CCX_FRACTURE_PAST_SEVERANCE), and every "
                 "increment after it is stamped as phantom");
        }
      }
      fflush(stdout);
    }
  }

  /* DE1.2/DM2.0 progressive damage is enabled for either the historical
     four-constant Rice-Tracey evolution record or the new type-3 tabulated
     ductile fracture-locus record.  Both use the same Newton-integrated
     evolution and DE1.3.1 terminal topology backend. */
  if((*ndmat_>0)&&(*iexpl<=1)){
    for(i=0;i<*nmat;i++){
      if(damage_progressive_material(i+1,ndmcon,dmcon,*ndmat_,*ntmat_)){
        damage_de12_matcount++;
        if((ITG)dmcon[1+(*ndmat_+1)*(*ntmat_)*i]==3)
          damage_dm20_matcount++;
      }
    }
    if(damage_de12_matcount>0) damage_de12_enabled=1;
    if(damage_de12_enabled){
      damage_tangent_env=getenv("CCX_DAMAGE_TANGENT");
      if((damage_tangent_env!=NULL)&&
         ((strcmp(damage_tangent_env,"FD_SYM")==0)||
          (strcmp(damage_tangent_env,"fd_sym")==0)||
          (strcmp(damage_tangent_env,"1")==0))){
        damage_tangent_mode=1;
      }else if((damage_tangent_env!=NULL)&&
               ((strcmp(damage_tangent_env,"UNSYM")==0)||
                (strcmp(damage_tangent_env,"unsym")==0)||
                (strcmp(damage_tangent_env,"2")==0))){

        /* Stage 1 of the consistent-tangent work: assemble and solve the
           bulk problem through the asymmetric path with the constitutive
           tangent left untouched.  resultsmech.f only ever tests
           de12tangent.eq.1, so mode 2 changes no material code. */

        damage_tangent_mode=2;
      }
      damage_reeq_scale_env=getenv("CCX_DAMAGE_REEQ_SCALE");
      if((damage_reeq_scale_env!=NULL)&&
         ((strcmp(damage_reeq_scale_env,"PHYSICAL")==0)||
          (strcmp(damage_reeq_scale_env,"physical")==0)||
          (strcmp(damage_reeq_scale_env,"1")==0))){
        damage_reeq_scale_mode=1;
      }
      /* J-09.  The same-load re-equilibration after a terminal deletion
         is tested with the STOCK RELATIVE correction criterion cam/uam, and
         REEQ_SCALE=PHYSICAL floors uam at the converged physical increment's
         displacement norm.  That floor SHRINKS WITH THE STEP, while the
         correction it is tested against does not: the perturbation is an
         element vanishing, not a load increment.  Measured on
         run_m12s3_rad - three consecutive failed attempts at dt = 2.64e-5,
         6.59e-6, 1.65e-6 (16x apart) give residuals agreeing to three
         digits, 0.548 / 0.550 / 0.551, while actual_uam holds at 1.43e-3
         and physical_ref collapses 2.37e-3 -> 1.48e-4.  Cutting the step is
         therefore the wrong instrument, and the loop can only end at tmin.

         This holds the reference at a fraction of the largest value it
         reached, exactly as CCX_DAMAGE_QAM_FLOOR does for the force side.
         The force side already had this; the displacement side did not.

         THIS CHANGES THE CONVERGENCE CRITERION AND THEREFORE THE ANSWER.
         A run that goes further with it is NOT thereby a success - that has
         to be shown on the physics (hydride consumption, failed facets,
         peakedness), and adopting it needs the full verify + ladder gate.
         Default 0 = off = bit-identical to the unpatched binary. */
      if((damage_de13_env=getenv("CCX_DAMAGE_REEQ_FLOOR"))!=NULL){
        damage_reeq_uam_floor=atof(damage_de13_env);
        if(damage_reeq_uam_floor<0.) damage_reeq_uam_floor=0.;
        if(damage_reeq_uam_floor>1.) damage_reeq_uam_floor=1.;
      }

      /* CCX_DAMAGE_RELEASE_PROBE - PURE DIAGNOSTIC, reads only.

         How much internal force does a topology event actually release?
         The run log cannot answer it: "largest residual force" is printed
         AFTER Newton has already made its first correction, so it is what
         survived the release, not the release.

         The probe differences f_int(u*) across the event with the
         displacement state held fixed, and - the reason it exists - splits
         the result by what happened to the equation:

           surv     DOF active BEFORE and AFTER.  The equation still exists,
                    so this is the only place a perturbation of the system
                    Newton solves can live.  dF_surv_max/qam is the number.
           removed  DOF active BEFORE, gone AFTER.  The equation does not
                    exist any more; Newton neither resolves it nor owes it
                    anything.  Diagnosis only, never a criterion.

         Reporting one number for both is what made the raw residual
         unreadable in the first place.

         =1 per batch, =2 adds one line per deleted element.  Unset = off and
         nothing is allocated.  This flag CHANGES NO BIT OF THE ANSWER, and
         that is gated both ways: s0_coarse_ts must give m.damage md5
         74212e957d7cd649 with the probe off AND with it on. */
      /* CCX_DAMAGE_RESIDUAL_RAY - PURE DIAGNOSTIC, reads only.

         J-13 established that the force RELEASED by a topology event does not
         order fatal against non-fatal: the fatal event ranked 487th of 563 on
         bandrad, and 486 larger releases were survived.  So the perturbation
         is not the right-hand side.  What is left is the STEP: the residual
         contracts three times (0.265 -> 0.019 -> 0.0022) and then flies up by
         265x on the next full Newton step, at frozen load and frozen topology.
         That is the classic signature of a full step leaving the basin, and
         the standard instrument for it is a scan of the residual along the
         Newton direction.

         There is no such scan in this tree, and there cannot be a line search
         either: BK3 is excluded from re-equilibration by TWO independent
         gates - its own conjunction (idamagereeq==0, below) and the resold
         store, which is guarded by the same conjunction, so resold is stale
         throughout an idamagereeq pass.  This probe needs neither: it
         evaluates alpha=0 itself and uses that as the reference.

         What it does: at a re-equilibration iteration it evaluates the
         residual at alpha=0 and alpha=1.  If the full step grew the norm by
         more than the growth factor - the fatal signature - it walks the
         whole ray, repeats one alpha to prove the evaluation is a pure
         function of the step length, and ends at alpha=1, which is exactly
         the state the unprobed code would have had.

         PURITY IS THE GATE, NOT AN ASIDE.  The repeated alpha must reproduce
         bitwise.  It can only do so for the residual vector, at fixed dtime,
         with CCX_DAMAGE_NONLOCAL unset and no contact: damjac/xstiff are not
         rebuilt from a baseline, and the nonlocal field (dpsave/ebar) is
         trial-derived, never snapshotted and CG-warm-started to 1e-10.  A
         mismatch under those conditions is a real impurity, not a bug in the
         probe, and it would void any line search built on top.

         Value = max number of rays to walk (default 8).  Unset = off,
         nothing allocated.  Changes no bit of the answer: b is saved and
         restored exactly, and the last evaluation is the alpha=1 state. */
      /* CCX_DAMAGE_REEQ_BACKTRACK - SOLVER CHANGE, not a diagnostic.
         Damps the Newton step during same-load re-equilibration, restoring
         the committed baseline before every probe and restoring the full
         step when nothing is acceptable.  THIS CHANGES THE ANSWER: a run
         that goes further with it is not thereby a success, and adopting it
         needs the full verify + ladder gate.  Default off = bit-identical. */
      /* CCX_DAMAGE_ABA=<alpha> - PURE DIAGNOSTIC.  Proves, or refutes, that a
         trial evaluation is a pure function of the step length.  The earlier
         A-B-A compared ONE scalar (|R|inf); one scalar agreeing proves
         nothing about the rest of the state, and the backtracking snapshot
         only held dam/damvisc/xstate while results() writes more than that.
         This evaluates A, snapshots EVERY array results()/calcresidual
         touch, evaluates B, evaluates A again, and compares byte for byte.
         Fires once, then the run continues from the full step. */
      if((damage_de13_env=getenv("CCX_DAMAGE_ABA"))!=NULL){
        damage_aba_mode=1;
        damage_aba_a=atof(damage_de13_env);
        if((damage_aba_a<=0.)||(damage_aba_a>=1.)) damage_aba_a=0.25;
        /* CCX_DAMAGE_ABA_INC="143,209" - fire at those increments instead of
           at the first opportunity.  Purity proved on one activated path does
           not prove it on another: a different increment reaches the same
           code through a different constitutive state, and that is exactly
           what has to be shown before an ensemble rests on it. */
        damage_aba_ninc=0;
        if((damage_de13_env=getenv("CCX_DAMAGE_ABA_INC"))!=NULL){
          char *acp=damage_de13_env;
          while((*acp!=0)&&(damage_aba_ninc<4)){
            while((*acp==' ')||(*acp==',')) acp++;
            if(*acp==0) break;
            damage_aba_inc[damage_aba_ninc++]=atoi(acp);
            while((*acp!=0)&&(*acp!=',')) acp++;
          }
        }
        printf("[DAMAGE ABA] DIAGNOSTIC: full-state A-B-A at alpha=%.6f;\n"
               "   every array written by results()/calcresidual is compared\n"
               "   byte for byte between two evaluations at the same alpha.%s",
               damage_aba_a,"\n");
        fflush(stdout);
      }

      if(getenv("CCX_DAMAGE_REEQ_BACKTRACK")!=NULL){
        damage_bt_mode=1;
        /* Three tunables, each aimed at a MEASURED failure of the
           first version (J-15 -> bandrad regressed 25%).
           _GROWTH : engage only when the full step makes the residual
                     worse by more than this factor.  Damping a step
                     that merely fails Armijo is what made the method
                     more aggressive than BK3 (which needs 1.10) and
                     is what stalled bandrad.  1.0 = old behaviour.
           _WINDOW : non-monotone reference (Grippo-Lampariello-
                     Lucidi).  Acceptance compares against the MAX of
                     the last WINDOW residuals, not the current one,
                     so Newton may worsen the residual briefly and
                     cross the kink - which is exactly what the
                     undamped control does.  1 = monotone = old.
           _FLOOR  : refuse to accept a step shorter than this.  The
                     measured death mode was a chain of accepts at
                     alpha=0.031 and 0.016 buying 1-3% each while the
                     iteration budget drained.  0.015625 = old. */
        if((damage_de13_env=getenv("CCX_DAMAGE_BT_GROWTH"))!=NULL){
          damage_bt_growth=atof(damage_de13_env);
          if(damage_bt_growth<1.) damage_bt_growth=1.;
        }
        if((damage_de13_env=getenv("CCX_DAMAGE_BT_WINDOW"))!=NULL){
          damage_bt_window=atoi(damage_de13_env);
          if(damage_bt_window<1) damage_bt_window=1;
          if(damage_bt_window>8) damage_bt_window=8;
        }
        if((damage_de13_env=getenv("CCX_DAMAGE_BT_FLOOR"))!=NULL){
          damage_bt_floor=atof(damage_de13_env);
          if(damage_bt_floor<0.015625) damage_bt_floor=0.015625;
          if(damage_bt_floor>1.) damage_bt_floor=1.;
        }
        printf("[DAMAGE BT] transactional backtracking ENABLED in "
               "idamagereeq: alpha 1, 1/2 ... 1/64, Armijo on |R|inf with "
               "c1=1e-4, committed baseline restored before every probe, "
               "full step restored and the increment handed to the standard "
               "cutback if no probe is acceptable.  THIS CHANGES THE "
               "ANSWER.  growth=%.3f window=%" ITGFORMAT
               " floor=%.6f%s",damage_bt_growth,damage_bt_window,
               damage_bt_floor,"\n");
        fflush(stdout);
      }

      /* ---- CCX_DAMAGE_REEQ_RESCUE ------------------------------------
         Emergency-only backtracking.  Always-on BT is EXPERIMENTAL and was
         measured to shorten solver survival on three placements of four and
         to destroy the bandrad severance the control reaches (J-17), so the
         two must never run together. */

      if((getenv("CCX_DAMAGE_REEQ_RESCUE")!=NULL)||
         (getenv("CCX_DAMAGE_REEQ_RESCUE2")!=NULL)||
         (getenv("CCX_DAMAGE_REEQ_RESCUE3")!=NULL)||
         (getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL)){
        damage_rescue_mode=1;
        ccx_rescue_active=1;
        if(getenv("CCX_DAMAGE_REEQ_RESCUE2")!=NULL){
          damage_rescue_maxlevel=2;
          damage_evt_nstep=0;
        }
        if((getenv("CCX_DAMAGE_REEQ_RESCUE3")!=NULL)||
           (getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL)){
          damage_evt_nstep=0;
          damage_reg_nlam=5;
          damage_rescue_maxlevel=2+damage_reg_nlam;
        }
        if(getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL){
          damage_corr_mode=1;
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_MAXINC"))!=NULL)
            damage_corr_maxinc=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_EXIT"))!=NULL)
            damage_corr_exit=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_TRY"))!=NULL)
            damage_corr_tryevery=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_GRACE"))!=NULL)
            damage_corr_grace=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_MAXWALL"))!=NULL)
            damage_corr_maxwall=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_MAXESC"))!=NULL)
            damage_corr_maxesc=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_STABLE"))!=NULL)
            damage_corr_stableneed=atoi(damage_de13_env);
          if((damage_de13_env=getenv("CCX_DAMAGE_CORR_MINFRAC"))!=NULL)
            damage_corr_minfrac=atof(damage_de13_env);
          if(damage_corr_maxinc<1) damage_corr_maxinc=1;
          if(damage_corr_exit<1) damage_corr_exit=1;
          if(damage_corr_tryevery<1) damage_corr_tryevery=1;
          if(damage_corr_maxwall<1) damage_corr_maxwall=1;
          if(damage_corr_maxesc<1) damage_corr_maxesc=1;
          if(damage_corr_stableneed<1) damage_corr_stableneed=1;
          printf("[DAMAGE CORR] bounded recovery CORRIDOR enabled.  On a "
                 "wall levels 1 and 2 cannot touch (idamagereeq=0) the "
                 "regularization that made the increment converge is HELD, "
                 "so the next increments start already regularized.  dtime "
                 "stays with the stock controller.  Three separate states "
                 "are kept: the lambda in use, the PROVEN lambda (one that "
                 "survived %" ITGFORMAT " converged increments) and a probe "
                 "flag.  Every %" ITGFORMAT " increments the help is probed "
                 "downwards (lambda/4, then 0).  A wall on a PROBE returns "
                 "to the proven lambda; a wall on the HELD lambda would be "
                 "an identical repeat, so lambda is escalated one ladder "
                 "step instead, at most %" ITGFORMAT " times, and the "
                 "corridor closes if the ladder runs out.  Wall counters "
                 "are evaluated AT THE WALL, so a chain of walls cannot "
                 "run unbounded.  Exit after %" ITGFORMAT " increments with "
                 "NO help.  Breakers: <=%" ITGFORMAT " walls, <=%" ITGFORMAT
                 " increments, and after a grace of %" ITGFORMAT " the mean "
                 "dtime inside must stay above %.3f of the dtime at the "
                 "last clean increment before entry.  On failure the "
                 "guarantee is t_end NOT LOWER than rescue-2; a byte-exact "
                 "rescue-2 result is impossible once corridor increments "
                 "have been accepted, since no entry snapshot is taken.%s",
                 damage_corr_stableneed,damage_corr_tryevery,
                 damage_corr_maxesc,damage_corr_exit,damage_corr_maxwall,
                 damage_corr_maxinc,damage_corr_grace,
                 damage_corr_minfrac,"\n");
          fflush(stdout);
        }
        if((damage_de13_env=getenv("CCX_DAMAGE_RESCUE_WINDOW"))!=NULL){
          damage_rec_window=atoi(damage_de13_env);
          if(damage_rec_window<1) damage_rec_window=1;
        }
        if((damage_de13_env=getenv("CCX_DAMAGE_RESCUE_MAXUNREC"))!=NULL){
          damage_rec_maxunrec=atoi(damage_de13_env);
          if(damage_rec_maxunrec<1) damage_rec_maxunrec=1;
        }
        printf("[DAMAGE RESCUE] bounded recovery window: a rescue counts "
               "as RECOVERED only after %" ITGFORMAT " consecutive "
               "increments converge with no intervention; after %"
               ITGFORMAT " consecutive un-recovered rescues the mechanism "
               "DISARMS itself and the wall goes to the original stock "
               "stop.  This exists because a run that needs rescuing at "
               "nearly every increment is crawling, not passing a wall: "
               "measured, 92 regularized rescues bought 2.2e-4 of step "
               "time on s3rad.%s",damage_rec_window,damage_rec_maxunrec,"\n");
        if(damage_bt_mode==1){
          printf("[DAMAGE RESCUE] CCX_DAMAGE_REEQ_BACKTRACK (always-on, "
                 "experimental) must not run together with rescue; it is "
                 "switched OFF for this run.%s","\n");
          damage_bt_mode=0;
        }
        printf("[DAMAGE RESCUE] emergency rescue backtracking ENABLED.  The "
               "trajectory, the stock Newton and every stock cutback are "
               "unchanged.  Only where the next stock cutback would put "
               "dtheta below tmin and the run would stop, the increment is "
               "rolled back by the STANDARD cutback path and retried ONCE at "
               "the last admissible dtheta with transactional BT active for "
               "that attempt alone.  One attempt per wall; re-armed after any "
               "increment that converges.%s","\n");
        if(damage_rescue_maxlevel==2){
          printf("[DAMAGE RESCUE2] second level ARMED.  A wall now gets two "
                 "attempts.  The first is the accepted level-one behaviour, "
                 "unchanged.  Only if it fails does the second run, and there "
                 "the single branch \"nothing accepted -> restore the full "
                 "step\" is replaced by an EVENT STEP: the smallest ladder "
                 "alpha at which the set of UC6 points in compression differs "
                 "from alpha=0, read from sign(stx(1)) over every live UC6 "
                 "point.  e_c3d_uc6.f:45 assembles the stiffness from vold, so "
                 "the next assembly picks up ctan(1,1)=kn on the crossed facet "
                 "by itself.  No constitutive law, no kn, no g and no material "
                 "parameter is touched, and no element, ip or increment is "
                 "named.  If no ladder alpha changes the set, this level does "
                 "nothing.%s","\n");
        }
        if(damage_reg_nlam>0){
          printf("[DAMAGE RESCUE3] third level ARMED with %" ITGFORMAT
                 " attempt(s): positive diagonal regularization K+lambda*D."
                 "  A wall whose failing solve has idamagereeq=0 goes "
                 "straight here - it never enters a same-load solve, so "
                 "levels 1 and 2 are gated out and would only repeat the "
                 "identical attempt.  ad[k] += lambda*D[k] with "
                 "D[k]=max(|ad[k]|,1e-6*mean|ad|) > 0, immediately before "
                 "the solver dispatch, where the existing stabiliser "
                 "already edits the same diagonal.  NOT ad*=(1+lambda), "
                 "which shifts only where the diagonal is positive; and D "
                 "is NOT plain |ad|, which is zero where the diagonal is "
                 "zero and, for ad<0, gives |ad|*(lambda-1) so lambda=1 "
                 "lands exactly on zero.  The per-attempt sign census is a "
                 "diagnostic of the diagonal, NOT evidence about the "
                 "definiteness of K.  Acceptance stays on the UNMODIFIED "
                 "residual: only ad is shifted.  If mean|ad| is zero, NaN "
                 "or infinite the shift is skipped and the attempt runs "
                 "stock.  lambda ladder:",damage_reg_nlam);
          for(i=0;i<damage_reg_nlam;i++) printf(" %.3e",damage_reg_lam[i]);
          printf(".  When it is exhausted the wall is left to the original "
                 "stock stop, so the run ends exactly where rescue-2 ends "
                 "it.%s","\n");
        }
        fflush(stdout);
      }

      if((damage_de13_env=getenv("CCX_DAMAGE_RESIDUAL_RAY"))!=NULL){
        damage_ray_probe=1;
        damage_ray_max=atoi(damage_de13_env);
        if(damage_ray_max<1) damage_ray_max=8;
        if(damage_ray_max>1000) damage_ray_max=1000;
        /* CCX_DAMAGE_RAY_INC="231" - walk only at these increments. */
        damage_ray_ninc=0;
        if((damage_de13_env=getenv("CCX_DAMAGE_RAY_INC"))!=NULL){
          char *rcp=damage_de13_env;
          while((*rcp!=0)&&(damage_ray_ninc<4)){
            while((*rcp==' ')||(*rcp==',')) rcp++;
            if(*rcp==0) break;
            damage_ray_inc[damage_ray_ninc++]=atoi(rcp);
            while((*rcp!=0)&&(*rcp!=',')) rcp++;
          }
        }
        printf("[DAMAGE RAY] DIAGNOSTIC: the residual is scanned along the "
               "Newton direction during same-load re-equilibration, at most "
               "%" ITGFORMAT " times.  Reads only; b is restored exactly and "
               "the scan ends on the full step, so no bit of the answer "
               "changes.%s",damage_ray_max,"\n");
        if(damage_ray_ninc>0){
          printf("[DAMAGE RAY] restricted to increments:");
          for(i=0;i<damage_ray_ninc;i++)
            printf(" %" ITGFORMAT,damage_ray_inc[i]);
          printf("%s","\n");
        }
        fflush(stdout);
      }

      if((damage_de13_env=getenv("CCX_DAMAGE_RELEASE_PROBE"))!=NULL){
        damage_release_probe=atoi(damage_de13_env);
        if(damage_release_probe<0) damage_release_probe=0;
        if(damage_release_probe>2) damage_release_probe=2;
        if(damage_release_probe>0){
          printf("[DAMAGE RELEASE] DIAGNOSTIC level %" ITGFORMAT
                 ": the internal force released by each topology event is "
                 "measured at frozen displacement and reported split into "
                 "surviving and removed degrees of freedom.  Reads only; "
                 "changes no bit of the answer.%s", damage_release_probe,
                 "\n");
          fflush(stdout);
        }
      }
      damage_linesearch_env=getenv("CCX_DAMAGE_LINESEARCH");
      if((damage_linesearch_env!=NULL)&&
         ((strcmp(damage_linesearch_env,"ADAPTIVE")==0)||
          (strcmp(damage_linesearch_env,"adaptive")==0)||
          (strcmp(damage_linesearch_env,"1")==0))){
        damage_linesearch_mode=1;
      }
      /* DE1.3 keeps an element until it is degraded to 0.1% of its
         stiffness.  Between D=0.95 and that point the element carries
         almost no load but still owns equations, and in a phase that is
         failing over a whole region there can be thousands of them at
         once.  Making the threshold settable lets that be measured
         instead of argued about. */

      /* ---- CCX_DAMAGE_TR_DOGLEG ---------------------------------------

         A root-finding TRUST REGION with a dogleg step, on the ORIGINAL
         equilibrium residual.  It is NOT another regularisation: no matrix
         is modified, no diagonal is shifted, no constitutive law, no Kn, no
         g and no deletion criterion is touched.  The only thing that changes
         is HOW LONG and IN WHICH DIRECTION the correction is, inside one
         armed increment attempt.

         Why here and not from the start: J-19 rejected the corridor because
         a held regularisation carried the solver instead of returning it to
         itself.  So this arms ONLY as rescue LEVEL 3, i.e. only after the
         Rescue2 levels 1 and 2 have both failed terminally on the same wall,
         and only on a wall with idamagereeq==0 - where the measured failure
         is a LINE SEARCH failure: at s3rad inc=569 BK3 reports lambda pinned
         at its floor 0.100000 with res_damped 2.001658e-03 ABOVE res_old
         1.942365e-03 on every one of the four identical attempts.  A floor
         of 0.1 on the Newton DIRECTION is exactly what a trust region does
         not have: it may go shorter, and it may leave that direction.

         Model, following PETSc SNESNEWTONTRDC:
             phi(u) = 1/2 |R(u)|^2 ,  R = f_int - f_ext = -b
             J p_N  = -R = b        (p_N is what PARDISO returns)
             g      = J^T R = -d    with d := J^T b
             p_C    = (|d|^2/|Jd|^2) d          (Cauchy point)
             p      = dogleg(p_C,p_N,Delta)
             rho    = [phi(u)-phi(u+p)] / [phi(u) - 1/2|R+Jp|^2]
         Acceptance of the STEP is rho; acceptance of the INCREMENT stays
         with checkconvergence on the unmodified residual.

         J^T IS FORMED AS J^T, and the run proves it: dot(d,p_N) must equal
         |b|^2 exactly, because dot(J^T b, J^-1 b) = b^T b.  That identity is
         printed as TRANSPOSE-CHECK on the first armed iteration and the
         mechanism REFUSES TO ARM if it is not 1 to 1e-8.  The asymmetry of
         the operator is measured at the same point, not assumed. */

      if(getenv("CCX_DAMAGE_TR_DOGLEG")!=NULL){
        if(damage_rescue_mode==0){
          printf("*ERROR: CCX_DAMAGE_TR_DOGLEG requires "
                 "CCX_DAMAGE_REEQ_RESCUE2; it is a level ON TOP of "
                 "Rescue2, not a replacement.  Stopping.%s","\n");
          fflush(stdout);FORTRAN(stop,());
        }
        if((damage_corr_mode==1)||(damage_reg_nlam>0)){
          printf("*ERROR: CCX_DAMAGE_TR_DOGLEG must not run together with "
                 "CCX_DAMAGE_REEQ_RESCUE3 or CCX_DAMAGE_RESCUE_CORRIDOR - "
                 "both were measured NEGATIVE (J-19) and both would occupy "
                 "the same rescue levels.  Stopping.%s","\n");
          fflush(stdout);FORTRAN(stop,());
        }
        if(damage_bt_mode==1){
          printf("*ERROR: CCX_DAMAGE_TR_DOGLEG must not run together with "
                 "always-on CCX_DAMAGE_REEQ_BACKTRACK (rejected, J-17).  "
                 "Stopping.%s","\n");
          fflush(stdout);FORTRAN(stop,());
        }
        damage_dl_mode=1;
        damage_rescue_maxlevel=3;
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_MAXTRIAL"))!=NULL)
          damage_dl_maxtrial=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_MAXEVAL"))!=NULL)
          damage_dl_maxeval=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_MAXFACT"))!=NULL)
          damage_dl_maxfact=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_MAXARM"))!=NULL)
          damage_dl_maxarm=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_D0"))!=NULL)
          damage_dl_d0fac=atof(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_LINCHECK"))!=NULL)
          damage_dl_lincheck=atoi(damage_de13_env);
        if(damage_dl_lincheck<0) damage_dl_lincheck=0;
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_LINCHECK_IT"))!=NULL)
          damage_dl_lc_it=atoi(damage_de13_env);
        if(damage_dl_lc_it<1) damage_dl_lc_it=1;
        if((damage_de13_env=getenv("CCX_DAMAGE_TR_LINCHECK_NIT"))!=NULL)
          damage_dl_lc_nit=atoi(damage_de13_env);
        if(damage_dl_lc_nit<1) damage_dl_lc_nit=1;
        /* the geometry of the step is proved before the first increment, on
           a J whose dogleg is known in closed form.  A failure here is
           arithmetic, so the run must not start. */
        if(damage_dl_selftest()!=0){
          printf("*ERROR: the trust-region geometry self-test FAILED.  "
                 "Stopping rather than running a method whose step "
                 "construction is wrong.%s","\n");
          fflush(stdout);FORTRAN(stop,());
        }
        if(damage_dl_maxtrial<1) damage_dl_maxtrial=1;
        if(damage_dl_maxtrial>12) damage_dl_maxtrial=12;
        if(damage_dl_maxeval<1) damage_dl_maxeval=1;
        if(damage_dl_maxfact<1) damage_dl_maxfact=1;
        if(damage_dl_maxarm<1) damage_dl_maxarm=1;
        if(damage_dl_d0fac<=0.) damage_dl_d0fac=1.;
        printf("[DAMAGE TR] trust-region DOGLEG armed as rescue LEVEL 3.  "
               "Levels 1 and 2 are the unchanged Rescue2 behaviour and run "
               "first; only when BOTH have failed on the same wall does the "
               "increment get one more attempt, and in that attempt every "
               "Newton correction is chosen by a dogleg trust region on "
               "phi=1/2|R|^2 instead of by BK3.  No matrix entry, no "
               "material constant and no deletion rule is touched, and "
               "convergence is still judged by checkconvergence on the "
               "UNMODIFIED residual.  Budget: <=%" ITGFORMAT " trial steps "
               "per iteration, <=%" ITGFORMAT " residual evaluations, <=%"
               ITGFORMAT " armed factorisations, <=%" ITGFORMAT " armed "
               "attempts in the whole run; on exhaustion the state is "
               "restored and the ORIGINAL stock stop runs.  Initial radius "
               "= %.3f * |p_Newton| at the first armed iteration.%s",
               damage_dl_maxtrial,damage_dl_maxeval,damage_dl_maxfact,
               damage_dl_maxarm,damage_dl_d0fac,"\n");
        fflush(stdout);
      }

      /* ---- CCX_DAMAGE_CONTINUATION (SPEC FREEZE v1) ------------------

         Bounded EXPERIMENTAL coupled local continuation.  Arms only as
         rescue LEVEL 4, i.e. only after Rescue2 levels 1 and 2 AND the
         dogleg have all failed on one wall.  Its single purpose is to find
         out whether coupled local continuation crosses the s3rad wall near
         inc=589 with physical front advance.

             R(u,lambda) = f(u,xbounact(lambda)) - fext = 0
             c(u,lambda) = m.(delta - delta_c) - ds    = 0

         lambda is a genuine unknown of a bordered system, not a corrected
         theta.  It owns the boundary for the rest of the step once armed.
         This is NOT a production continuation: there is no terminal landing
         at lambda=1, no return to stock control, no completed step and no
         restart.  Every ending is PARTIAL. */

      if(getenv("CCX_DAMAGE_CONTINUATION")!=NULL){
        if((damage_rescue_mode==0)||(damage_dl_mode==0)){
          printf("*ERROR: CCX_DAMAGE_CONTINUATION requires BOTH "
                 "CCX_DAMAGE_REEQ_RESCUE2 and CCX_DAMAGE_TR_DOGLEG; it is a "
                 "level ON TOP of them, never a replacement.  Stopping.%s",
                 "\n");
          fflush(stdout);FORTRAN(stop,());
        }
        if((damage_corr_mode==1)||(damage_reg_nlam>0)||(damage_bt_mode==1)||
           (damage_arc==1)||(damage_diss_ctrl>=1)||(damage_path_on>0)){
          printf("*ERROR: CCX_DAMAGE_CONTINUATION conflicts with "
                 "CCX_DAMAGE_ARCLENGTH, CCX_DISSIPATION_CONTROL, "
                 "CCX_DAMAGE_PATH, CCX_DAMAGE_REEQ_RESCUE3, "
                 "CCX_DAMAGE_RESCUE_CORRIDOR and always-on "
                 "CCX_DAMAGE_REEQ_BACKTRACK.  None of them is used as a "
                 "foundation and simultaneous operation is refused.  "
                 "Stopping.%s","\n");
          fflush(stdout);FORTRAN(stop,());
        }
        damage_ct_mode=1;
        damage_rescue_maxlevel=4;
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_RHOMIN"))!=NULL)
          damage_ct_rhomin=atof(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_CLIM"))!=NULL)
          damage_ct_clim=atof(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_ULIM"))!=NULL)
          damage_ct_ulim=atof(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_EPS"))!=NULL)
          damage_ct_eps=atof(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_KAPTOL"))!=NULL)
          damage_ct_kaptol=atof(damage_de13_env);
        if(damage_ct_kaptol<1.) damage_ct_kaptol=1.5;
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_MAXSTEP"))!=NULL)
          damage_ct_maxstep=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_MAXCORR"))!=NULL)
          damage_ct_maxcorr=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_MAXFACT"))!=NULL)
          damage_ct_maxfact=atoi(damage_de13_env);
        if((damage_de13_env=getenv("CCX_DAMAGE_CT_MAXEVAL"))!=NULL)
          damage_ct_maxeval=atoi(damage_de13_env);
        if(damage_ct_rhomin<=0.) damage_ct_rhomin=1.e-4;
        if(damage_ct_clim<=0.) damage_ct_clim=20.;
        if(damage_ct_ulim<=0.) damage_ct_ulim=20.;
        if(damage_ct_eps<=0.) damage_ct_eps=1.e-6;
        if(damage_ct_maxstep<1) damage_ct_maxstep=1;
        if(damage_ct_maxcorr<1) damage_ct_maxcorr=1;
        if(damage_ct_maxfact<1) damage_ct_maxfact=1;
        if(damage_ct_maxeval<1) damage_ct_maxeval=1;
        printf("[DAMAGE CT] bounded EXPERIMENTAL continuation armed as rescue "
               "LEVEL 4.  Levels 1-3 (Rescue2 and the dogleg) are unchanged "
               "and run first.  On a wall none of them takes, lambda becomes "
               "a genuine unknown of a bordered system with one FROZEN "
               "local mixed-mode UC6 constraint, and owns the boundary for "
               "the rest of the step.  There is NO terminal landing, NO "
               "return to stock control, NO completed step and NO restart: "
               "every ending is PARTIAL.  Parameters (all opt-in, printed as "
               "actually used): rho_den_min=%.3e C_lambda=%.1f C_u=%.1f "
               "eps_FD=%.3e kappa_tol=%.3f; budget <=%" ITGFORMAT " steps, <=%" ITGFORMAT
               " corrector iterations, <=%" ITGFORMAT " factorisations, <=%"
               ITGFORMAT " residual evaluations.%s",
               damage_ct_rhomin,damage_ct_clim,damage_ct_ulim,damage_ct_eps,
               damage_ct_kaptol,
               damage_ct_maxstep,damage_ct_maxcorr,damage_ct_maxfact,
               damage_ct_maxeval,"\n");
        fflush(stdout);
        if(damage_ct_selftest()!=0){
          printf("*ERROR: the continuation bordered-algebra self-test "
                 "FAILED.  Stopping rather than running a method whose "
                 "constraint row is wrong.%s","\n");
          fflush(stdout);FORTRAN(stop,());
        }
      }

      /* A requested rescue flag that silently fails to arm has already
         cost two whole runs: CCX_DAMAGE_REEQ_RESCUE3 and then
         CCX_DAMAGE_RESCUE_CORRIDOR were each parsed inside a block whose
         outer condition did not list them, so the run executed as the
         plain control and read like a negative result.  Shout instead of
         pretending. */
      if((getenv("CCX_DAMAGE_RESCUE_CORRIDOR")!=NULL)&&
         (damage_corr_mode==0)){
        printf("*ERROR: CCX_DAMAGE_RESCUE_CORRIDOR is set but the "
               "corridor did NOT arm; the run would silently be a plain "
               "control.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if((getenv("CCX_DAMAGE_REEQ_RESCUE3")!=NULL)&&
         (damage_reg_nlam==0)){
        printf("*ERROR: CCX_DAMAGE_REEQ_RESCUE3 is set but the "
               "regularized level did NOT arm.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if((getenv("CCX_DAMAGE_CONTINUATION")!=NULL)&&
         ((damage_ct_mode==0)||(damage_rescue_maxlevel!=4))){
        printf("*ERROR: CCX_DAMAGE_CONTINUATION is set but the continuation "
               "level did NOT arm; the run would silently be a plain "
               "Rescue2+dogleg control.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if((getenv("CCX_DAMAGE_TR_DOGLEG")!=NULL)&&
         ((damage_dl_mode==0)||(damage_rescue_maxlevel<3))){
        printf("*ERROR: CCX_DAMAGE_TR_DOGLEG is set but the dogleg level "
               "did NOT arm; the run would silently be a plain Rescue2 "
               "control.  Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }
      if(((getenv("CCX_DAMAGE_REEQ_RESCUE")!=NULL)||
          (getenv("CCX_DAMAGE_REEQ_RESCUE2")!=NULL))&&
         (damage_rescue_mode==0)){
        printf("*ERROR: a rescue flag is set but rescue did NOT arm.  "
               "Stopping.%s","\n");
        fflush(stdout);FORTRAN(stop,());
      }

      /* CCX_FRACTURE_DEADFACET is retired.  It made "a fully failed
         cohesive facet is not a load path" opt-in, and the one runner that
         armed the connectivity test - test/s3rad/run_s3rad.sh - explicitly
         unset it, which is how a 42807-element specimen came to be driven
         178 increments past its own severance.  A judgement that is only
         correct when somebody remembers to ask for it is not a judgement.
         loadpath.c now applies it unconditionally, and damstate_facet_dead
         still owns the "every integration point has failed" test. */

      /* CCX_DAMAGE_FACET_DELETE - the same failure flag as
         CCX_FRACTURE_DEADFACET, but acted on in the EQUATIONS instead of only
         in the termination test.

         The bulk phases create a displacement discontinuity when they fail;
         the interface never does, because g = max(gmin, 1-dvisc) floors the
         facet stiffness and a fully debonded facet keeps conducting gmin*Kn
         for ever.  Measured directly: 45.3% of the radial twin's backed
         facets sit at Dmin >= 0.999 with ZERO removed (E-113).  E-111 item 1
         named this as the remaining gap once the internal length (E-114) and
         the penalty stiffness (E-113) had both been tested and excluded, and
         E-119/E-120 localised the unresolved route to DEBONDING - the one
         route whose crack this prevents.

         Unlike CCX_FRACTURE_DEADFACET this DOES change the equations, so it
         carries the risk that sank damfloatface (E-57, E-72) and is gated the
         same way: default OFF, batch-capped, and travelling the unchanged
         transactional path so rollback restores it (E-64).  It must not be
         adopted without verify 65/65 and the full ladder. */
      if((damage_de13_env=getenv("CCX_DAMAGE_FACET_DELETE"))!=NULL){
        damage_facetdel=(strcmp(damage_de13_env,"0")==0)?0:1;
        if(damage_facetdel){
          printf("[DAMAGE FACET DELETE] a cohesive facet whose every "
                 "integration point has failed is REMOVED from the equations "
                 "(its compressive penalty is removed with it)\n");
        }
      }

      if(getenv("CCX_DISSIPATION_REPORT")!=NULL) damage_diss_report=1;
      damage_diss_env=getenv("CCX_DISSIPATION_TARGET");
      if(damage_diss_env!=NULL){
        damage_diss_target=atof(damage_diss_env);
        if(damage_diss_target>0.) damage_diss_report=1;
      }
      if(getenv("CCX_DISSIPATION_PROBE")!=NULL) damage_diss_probe=1;
      if(getenv("CCX_DAMAGE_BATCH_LIST")!=NULL) damage_batch_list=1;
      /* The stress is scaled by 1-Dvis whenever the viscosity is
         on, so the deletion trigger has to read the same variable.
         Measured on DHC1: elements were being removed while still
         carrying up to 41% of their effective stress (batch_Dvis
         down to 0.592), and correcting it moved the run from
         lambda=0.3793 to 0.4651.  With the viscosity off damvisc
         is never allocated and the trigger falls back to D, so
         this is a no-op there - verified on SP1. */
      if((damage_de13_env=getenv("CCX_DAMAGE_DELETE_VISC"))!=NULL){
        damage_delete_visc=(strcmp(damage_de13_env,"0")==0)?0:1;
      }
      if((damage_de13_env=getenv("CCX_STRUCT_FD_INC"))!=NULL)
        damage_fd_inc=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_STRUCT_FD_ITER"))!=NULL)
        damage_fd_it=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_STRUCT_FD_H"))!=NULL)
        damage_fd_h=atof(damage_de13_env);
      if(getenv("CCX_DAMAGE_FREE_PROBE")!=NULL) damage_free_probe=1;
      if((damage_de13_env=getenv("CCX_DAMAGE_NODE_DUMP"))!=NULL)
        damage_dump_node=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_DAMAGE_NODE_INC"))!=NULL)
        damage_dump_inc=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_DAMAGE_NULLVEC"))!=NULL)
        damage_null_inc=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_DAMAGE_NULLVEC_IT"))!=NULL)
        damage_null_nit=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_DAMAGE_PATH_DROP"))!=NULL)
        damage_path_drop=atof(damage_de13_env);
      if((damage_de13_env=getenv("CCX_DAMAGE_PATH_NSTEP"))!=NULL)
        damage_path_nstep=atoi(damage_de13_env);
      if((damage_de13_env=getenv("CCX_DAMAGE_PATH"))!=NULL){
        damage_path_on=atoi(damage_de13_env);
        if(damage_path_on<0) damage_path_on=0;
        if(damage_path_on>2) damage_path_on=2;
      }
      if(getenv("CCX_DAMAGE_STIFF_PROBE")!=NULL)
        damage_stiff_probe=1;
      if(getenv("CCX_DAMAGE_TANGENT_CENSUS")!=NULL)
        damage_unsym_census=1;
      /* resultsmech.f reads this name too and has no safe place to print
         from - it runs on several threads.  Reporting it here is what makes
         the difference between "the flag ran and changed nothing" and "the
         flag was never read", which is exactly the ambiguity that made the
         first CCX_DAMAGE_TANGENT_FULL A/B uninformative (J-10). */
      if(getenv("CCX_DAMAGE_TANGENT_FULL")!=NULL){
        damage_unsym_tanfull=1;
        printf("[DAMAGE TANGENT FULL] the consistent-tangent cut-off is "
               "taken on the SAME damage variable the stress uses instead "
               "of the stock D<0.999.  This CHANGES THE OPERATOR and "
               "therefore the answer\n");
        fflush(stdout);
      }
      /* Nonlocal damage: an internal length in the FORMULATION.
         Crack-band scales the dissipation and is correct, but it gives
         the localisation no width - three uniform bar meshes agree to
         0.94% before the peak and spread to 73-89% after it, and
         viscosity changes solvability without moving those curves at
         all (E-83).  Default off, so every stored baseline stands. */
      if((damage_de13_env=getenv("CCX_DAMAGE_NONLOCAL"))!=NULL){
        damage_nl_ell=atof(damage_de13_env);
        if(damage_nl_ell<0.) damage_nl_ell=0.;
        if(damage_nl_ell>0.){
          printf("[DAMAGE NONLOCAL] the damage driving variable is averaged over a neighbourhood of radius 2*ell with ell=%.4e; the stress update stays local\n",damage_nl_ell);
        }
      }
      FORTRAN(damnonlocalset,(&damage_nl_ell));
      /* Which regularisation backend computes that average.

         INTEGRAL is the default and stays the default: it is what E-84
         and E-85 measured, so `CCX_DAMAGE_NONLOCAL=<ell>` alone must keep
         meaning exactly what those records say it means.

         GRADIENT solves  ebar - div(ell^2 grad ebar) = e  instead.  Same
         internal length, no neighbour list: the integral form stores
         O(N*(2*ell/h)^3) neighbours, which is 6.5 GB at ell=0.15 and
         58 GB at ell=0.45 on demo_realistic_clusters3 against 31.6 GB of
         machine, because one hydride element's neighbourhood holds
         thousands of tiny elements.  The PDE form is O(nnz) on any mesh
         and is the only one that can be used on a graded mesh at all. */
      if((damage_de13_env=getenv("CCX_DAMAGE_NONLOCAL_MODE"))!=NULL){
        if((strcmp(damage_de13_env,"GRADIENT")==0)||
           (strcmp(damage_de13_env,"gradient")==0)){
          damage_nl_mode=1;
        }else if((strcmp(damage_de13_env,"INTEGRAL")==0)||
                 (strcmp(damage_de13_env,"integral")==0)){
          damage_nl_mode=0;
        }else if((strcmp(damage_de13_env,"FROZEN")==0)||
                 (strcmp(damage_de13_env,"frozen")==0)){
          /* Not a backend: the CONTROL that separates the internal length
             from the integration scheme.  With ell>0 the driving variable
             is refreshed only at the predictor and the commit, so damage
             stops being integrated implicitly inside Newton - a change of
             scheme, not a lag.  Every local-vs-nonlocal comparison in this
             tree therefore moves two factors at once.  FROZEN supplies the
             same staggered update with NO averaging and NO length. */
          damage_nl_mode=2;
        }else{
          printf("*ERROR: CCX_DAMAGE_NONLOCAL_MODE must be INTEGRAL, GRADIENT or FROZEN, not %s\n",damage_de13_env);
          FORTRAN(stop,());
        }
      }
      if((damage_nl_ell>0.)&&(damage_nl_mode==1)){
        printf("[DAMAGE NONLOCAL] backend GRADIENT: implicit gradient, ebar - div(ell^2 grad ebar) = e, lumped mass, Jacobi-CG; no neighbour list\n");
      }
      if((damage_nl_ell>0.)&&(damage_nl_mode==2)){
        printf("[DAMAGE NONLOCAL] backend FROZEN: CONTROL ONLY - the staggered damage update WITHOUT spatial averaging.  ell is read but NOT applied; this arm regularises nothing and exists to separate the internal length from the integration scheme\n");
      }
      FORTRAN(damnonlocalmode,(&damage_nl_mode));
      if((damage_de13_env=getenv("CCX_DAMAGE_QAM_FLOOR"))!=NULL){
        damage_qam_floor=atof(damage_de13_env);
        if(damage_qam_floor<0.) damage_qam_floor=0.;
        if(damage_qam_floor>1.) damage_qam_floor=1.;
        if(damage_qam_floor>0.){
          printf("[DAMAGE QAM FLOOR] DIAGNOSTIC: the reference force qam is held at or above %.3e of its running maximum, so the RELATIVE force criterion cannot collapse as the specimen unloads.  This CHANGES THE CONVERGENCE CRITERION and therefore the answer; going further with it is not by itself a success\n",damage_qam_floor);
        }
      }
      if(getenv("CCX_DAMAGE_AUTOSPC_NEG")!=NULL) damage_spc_neg=1;
      if((damage_de13_env=getenv("CCX_DAMAGE_AUTOSPC"))!=NULL){
        damage_spc_g=atof(damage_de13_env);
        if(damage_spc_g<0.) damage_spc_g=0.;
        if(damage_spc_g>1.e-1) damage_spc_g=1.e-1;
        if(damage_spc_g>0.){
          damage_stiff_probe=1;
          /* [DAMSTATE] prove the judgement before letting it decide
             anything, on every run, the discipline lsladder.c set. */
          if(damstate_selftest()!=0){
            printf("[DAMSTATE] *ERROR: the load-path judgement self test "
                   "failed; disabling AUTOSPC rather than judging nodes by a "
                   "rule that is not the one that was tested.\n");
            damage_spc_g=0.;
            damage_stiff_probe=0;
          }else{
            printf("[DAMAGE AUTOSPC] a node whose assembled diagonal has fallen below %.1e of its own intact value is excluded from the DISPLACEMENT convergence norm; the force residual is untouched and nothing is deleted\n",damage_spc_g);
          }
        }
      }
      if((damage_de13_env=getenv("CCX_DAMAGE_STIFF_MIN"))!=NULL){
        damage_stiff_min=atof(damage_de13_env);
        if(damage_stiff_min<0.) damage_stiff_min=0.;
        if(damage_stiff_min>0.1) damage_stiff_min=0.1;
        if(damage_stiff_min>0.) damage_stiff_probe=1;
      }
      /* Default 0 since E-22.  damdangle cannot distinguish a dangling
         sliver from the legitimate last element at a node, so it deletes
         healthy load-bearing material and on a small mesh removes every
         element, after which remastruct runs on an empty model and the
         process dies.  Measured cost: R2 and R3S PASS->FAIL,
         nc2_two_tet_wp and bk4_orphan_tet_wp access violation.  Measured
         benefit on DHC1: +0.0009, down from the +0.021 of E-04, which the
         E-05 deletion-on-Dtilde fix superseded.  Set to 1 only to
         reproduce the old behaviour. */
      if((damage_de13_env=getenv("CCX_DAMAGE_DANGLE"))!=NULL){
        damage_dangle_max=atoi(damage_de13_env);
        if(damage_dangle_max<0) damage_dangle_max=0;
        if(damage_dangle_max>4) damage_dangle_max=4;
      }
      damage_deadall_env=getenv("CCX_DAMAGE_DEADALL");
      if(damage_deadall_env!=NULL){
        damage_deadall_g=atof(damage_deadall_env);
        if(damage_deadall_g<0.) damage_deadall_g=0.;
        if(damage_deadall_g>0.5) damage_deadall_g=0.5;
        if(damage_deadall_g>0.){
          printf("[DAMAGE DEADALL] a node whose ENTIRE live support is dead "
                 "(g < %.1e) and which no cohesive facet holds has that "
                 "support deleted; only dead elements can be taken\n",
                 damage_deadall_g);
        }
      }
      damage_deadsole_env=getenv("CCX_DAMAGE_DEADSOLE");
      if(damage_deadsole_env!=NULL){
        damage_deadsole_g=atof(damage_deadsole_env);
        if(damage_deadsole_g<0.) damage_deadsole_g=0.;
        if(damage_deadsole_g>0.5) damage_deadsole_g=0.5;
        if(damage_deadsole_g>0.){
          printf("[DAMAGE DEADSOLE] a dead element (g < %.1e) that is the sole "
                 "support of a node is deleted\n",damage_deadsole_g);
        }
      }
      damage_stab_env=getenv("CCX_DAMAGE_STABILISE");
      if(damage_stab_env!=NULL){
        damage_stab_alpha=atof(damage_stab_env);
        if(damage_stab_alpha<0.) damage_stab_alpha=0.;
        if(damage_stab_alpha>1.e-1) damage_stab_alpha=1.e-1;
        if(damage_stab_alpha>0.){
          printf("[DAMAGE STABILISE] detached pieces held with alpha=%.3e "
                 "of the mean diagonal; nothing is deleted\n",
                 damage_stab_alpha);
        }
      }
      damage_diss_env=getenv("CCX_DISSIPATION_CONTROL");
      if((damage_diss_env!=NULL)&&(damage_diss_target>0.)){
        damage_diss_ctrl=(strcmp(damage_diss_env,"2")==0)?2:1;
        if(damage_diss_ctrl==2){
          NNEW(damage_diss_fhat,double,neq[1]);
          NNEW(damage_diss_uf,double,neq[1]);
        }
        damage_diss_report=1;
        printf("[DISSIPATION CONTROL] load factor solved from "
               "dG=%.6e per increment; theta is advanced only on "
               "acceptance\n",damage_diss_target);
      }

      /* Path following: give lambda an identity separate from theta.

         Everything in this file has so far identified the load factor WITH
         the step time: `damage_diss_lprev=theta`, and both dissipation
         schemes end by clamping `lamcur` back to `theta`
         (nonlingeo.c 5153 and 5635).  theta only ever advances, because
         dtime=dtheta*tper feeds the viscous update and must stay positive,
         so lambda could be pulled back INSIDE an increment and never end
         below the last converged value.  That is why the file's own comment
         says a genuine snap-back stays out of reach.

         The separation is the whole change:

             theta   monotone pseudo-time, still drives dtime
             lambda  the load factor, free to decrease

         and xbounact is built from lambda, which is legitimate because
         xbounold holds the value at the START of the step, so
         xbounold+(xboun-xbounold)*lambda is an absolute step fraction -
         the equivalence CCX_DAMAGE_PATH mode 1 was written to verify and
         did, at max|ramp-tempload| = 0.

         lambda still needs an EQUATION, and that is the bordered
         dissipation solve, so this requires CCX_DISSIPATION_CONTROL=2.
         Without it there is nothing to determine lambda and the flag is
         refused rather than silently doing something else.

         Measured motivation: soft50 with a released interface stops at 8.8%
         of peak on a force residual that fails at 81% NORMAL nodes, with a
         stable average force and a clean topology - a limit point, not an
         artefact (E-95).  m12_spc3 stops the same way at a healthy grip
         node (E-74). */
      /* The target does double duty: it sizes the step AND, through
         DAMAGE_DISS_ENGAGE, decides when the constraint takes over.  Those
         two pull in opposite directions - a target matched to the real
         dissipation rate engages while the response is still elastic, which
         is precisely the state the coupled solve cannot represent (dG is
         identically zero there).  Measured on the released-interface soft50:
         target 0.30 against a real post-peak rate of 0.75-2.2 per increment
         engaged at increment 9, lambda=0.0155, and died at increment 10.

         So give engagement its own control: stay in displacement control
         until the step time passes this value, whatever the dissipation is
         doing.  Set it past the load peak. */
      /* Separate the two halves of dissipation control.

         CCX_DISSIPATION_TARGET drives THREE things at once: it sizes the
         next increment, it gates engagement through DAMAGE_DISS_ENGAGE, and
         it is the right-hand side of the coupled constraint.  Measured on
         the released-interface soft50, a target matched to the real
         dissipation rate (0.75-2.2 per increment) makes the step controller
         grow the step by its 1.25 cap every increment and thrash into tmin
         at the load peak: step range 581x against the stock controller's
         64x, floor 2.83e-6 against 3.13e-5, and the run died with ZERO
         deletions before the constraint ever engaged (E-96).

         CCX_DISSIPATION_STEP=0 keeps the stock step controller and leaves
         the target to the constraint alone. */
      if((damage_de13_env=getenv("CCX_DISSIPATION_STEP"))!=NULL){
        damage_diss_step=(strcmp(damage_de13_env,"0")==0)?0:1;
        if(damage_diss_step==0){
          printf("[DISSIPATION] step sizing is OFF; the target drives the "
                 "constraint only, the stock controller sizes the step\n");
        }
      }
      if((damage_de13_env=getenv("CCX_DISSIPATION_ENGAGE_T"))!=NULL){
        damage_diss_engage_t=atof(damage_de13_env);
      }
      if((damage_de13_env=getenv("CCX_DAMAGE_ARCLENGTH"))!=NULL){
        damage_arc=(strcmp(damage_de13_env,"0")==0)?0:1;
      }
      if(damage_arc==1){
        if(damage_diss_ctrl!=2){
          printf("[PATH FOLLOWING] *WARNING: CCX_DAMAGE_ARCLENGTH needs "
                 "CCX_DISSIPATION_CONTROL=2 to have an equation for lambda; "
                 "it is ignored\n");
          damage_arc=0;
        }else{
          printf("[PATH FOLLOWING] lambda is decoupled from the step time "
                 "and may DECREASE; theta stays monotone for dtime\n");
        }
      }

      damage_delete_filter=getenv("CCX_DAMAGE_DELETE_MAT");

      damage_de13_env=getenv("CCX_DAMAGE_DELETE_D");
      if(damage_de13_env!=NULL){
        damage_de13_delete_d=atof(damage_de13_env);
        if((damage_de13_delete_d<0.5)||(damage_de13_delete_d>0.9999)){
          printf("[DAMAGE DE1.3.1] *WARNING: CCX_DAMAGE_DELETE_D=%s is "
                 "outside [0.5,0.9999]; keeping %.4f\n",
                 damage_de13_env,(double)DAMAGE_DE13_DELETE_D);
          damage_de13_delete_d=DAMAGE_DE13_DELETE_D;
        }
      }

      damage_topology_env=getenv("CCX_DAMAGE_TOPOLOGY");
      if((damage_topology_env!=NULL)&&
         ((strcmp(damage_topology_env,"DEFERRED")==0)||
          (strcmp(damage_topology_env,"deferred")==0)||
          (strcmp(damage_topology_env,"1")==0))){
        damage_topology_deferred_mode=1;
      }
    }
  }

  /* Path following is armed HERE, outside the *DAMAGE INITIATION block
     above.  It was inside it, which silently disabled it for any model
     without a damage material - including the cohesive benchmark, where
     the only nonlinearity is a UC6 interface and ndmat_ is zero.  The
     method has nothing to do with whether a bulk damage material exists. */

  /* ---- CCX_PATHFOLLOW ------------------------------------------
     A self-contained dissipation path following, independent of the
     CCX_DISSIPATION_CONTROL machinery above and of PARDISO.  It is
     armed only inside a domain where the bordered system is exactly
     the one pathfollow.c verifies, and it refuses to arm otherwise
     rather than degrade silently. */

  if(getenv("CCX_DAMAGE_BATCH_TRACE")!=NULL) td_trace=1;
  if(getenv("CCX_TOPODIAG")!=NULL) td_from=atoi(getenv("CCX_TOPODIAG"));
  if(getenv("CCX_DAMAGE_LS_PROBE")!=NULL) damage_ls_probe=1;
  if(getenv("CCX_DAMAGE_WALL_NULLVEC")!=NULL) damage_wall_null=1;

  /* [WALLDIAG] Arm the increment-numbered probes by LOAD FACTOR instead.

     The wall is reproducible in theta and NOT in the increment index: the
     same binary on the same deck lands on it at increment 348 with two
     threads and 353 with four, because MKL_CBWR fixes the instruction set
     and not the reduction order.  Naming an increment therefore arms the
     diagnostics at a state that the next run does not have.  theta is the
     physical coordinate the wall is recorded in, so the gate takes it. */

  if(getenv("CCX_DAMAGE_WALL_THETA")!=NULL){
    damage_wall_theta=atof(getenv("CCX_DAMAGE_WALL_THETA"));
    printf("[WALLDIAG] armed at theta >= %.9e: from the first attempt at or "
           "beyond that load factor the linearisation check, the topology "
           "diagnostic, the deflated null-vector probe and the line-search "
           "ladder probe run on the CURRENT increment, whatever its "
           "number.%s",damage_wall_theta,"\n");
    fflush(stdout);
  }
  /* [WALLDIAG] MASKSTEP.  The wall measurement says 99.96% of |p_N|^2 sits
     on the three dofs of ONE node whose assembled diagonal has collapsed,
     while the force imbalance is two nodes away.  This probe walks the same
     eps ladder along p_N with every component on an AUTOSPC-masked node
     zeroed, so the question "is the step unusable BECAUSE it is spent on
     collapsed-diagonal nodes" becomes a number instead of a story.  It reads
     the mask AUTOSPC already builds, changes no solution path, and is off
     unless asked for. */
  /* [DAMAGE AUTOSPC-FORCE] AUTOSPC already decides which nodes have lost
     their load path - assembled diagonal below CCX_DAMAGE_AUTOSPC times the
     node's OWN intact value - and removes them from the DISPLACEMENT
     convergence norm only, deliberately leaving the force residual alone so
     that "a node that still carries load still blocks convergence".

     At the second wall that reasoning inverts.  The peak force residual sits
     on a node held by ONE live bulk element whose six cohesive facets have
     all failed (g -> 0), it needs a displacement of order one to be
     equilibrated, and it moves 0.1% per Newton iteration.  It does not carry
     load; it blocks convergence anyway, and it blocks the very increment in
     which its own last element would damage and delete.  That is a deadlock:
     the solve cannot advance because of the fragment, and the fragment
     cannot resolve because the solve cannot advance.

     This lifts the veto for exactly the AUTOSPC set and for nothing else.
     The dof is still assembled, still solved and still moved - only its
     veto over ram[0] is removed - and every excluded residual is printed, so
     it can never hide a growing imbalance.  OFF unless asked for. */
  if(getenv("CCX_DAMAGE_AUTOSPC_FORCE")!=NULL){
    damage_spc_force=1;
    printf("[DAMAGE AUTOSPC-FORCE] armed: a node already masked by AUTOSPC "
           "is excluded from the FORCE residual ram[0] as well as from "
           "cam[0].  Its dof is still solved and still moved; the largest "
           "excluded residual is reported on every increment.%s","\n");
    fflush(stdout);
  }
  if(getenv("CCX_DAMAGE_WALL_MASKSTEP")!=NULL){
    damage_wall_maskstep=1;
    printf("[WALLDIAG] MASKSTEP armed: the linearisation check gains a third "
           "pass along p_N with the AUTOSPC-masked nodes' components zeroed. "
           " Diagnostic only - nothing on a solution path reads it.%s","\n");
    fflush(stdout);
  }
  if((td_trace!=0)||(td_from>0)){
    if(topodiag_selftest()!=0){
      printf("[TOPODIAG] *ERROR: self test failed; diagnostics disabled\n");
      td_trace=0;td_from=0;
    }else{
      td_armed=1;
      NNEW(td_comp,ITG,*nk);
      NNEW(td_sort,ITG,(*ne>0)?*ne:1);
      printf("[TOPODIAG] armed: batch trace %s, topology report from "
             "increment %" ITGFORMAT ".  Nothing here changes the "
             "calculation.\n",td_trace?"ON":"off",td_from);
    }
  }

  pf_env=getenv("CCX_PATHFOLLOW");
  if(pf_env!=NULL){
    pf_tauv=atof(pf_env);
    if(!(pf_tauv>0.)){
      printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW must be a positive "
             "dissipation increment; got \"%s\".  Not armed.\n",pf_env);
    }else if(damage_diss_ctrl>=1){
      printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW and "
             "CCX_DISSIPATION_CONTROL both drive the load factor; "
             "set only one.  Not armed.\n");
    }else if((*nmethod!=1)||(*ithermal>=2)||(*mortar>1)||(ncont!=0)||
             (*iexpl>1)||(*nboun<=0)||
             ((*isolver!=0)&&(*isolver!=7))){
      printf("[PATHFOLLOW] not armed: outside the verified domain "
             "(nmethod=%" ITGFORMAT " ithermal=%" ITGFORMAT " mortar=%"
             ITGFORMAT " ncont=%" ITGFORMAT " iexpl=%" ITGFORMAT
             " nboun=%" ITGFORMAT " isolver=%" ITGFORMAT
             "; needs static, ithermal<2, no contact, implicit, "
             "prescribed dofs, SPOOLES or PARDISO)\n",
             *nmethod,*ithermal,*mortar,ncont,*iexpl,*nboun,*isolver);
    }else if(pathfollow_selftest()!=0){
      printf("[PATHFOLLOW] *ERROR: the bordered-algebra self test "
             "failed; refusing to arm.\n");
    }else if(pathfollow_arm(pf_tauv,neq[1])==0){
      printf("[PATHFOLLOW] *ERROR: could not allocate; not armed.\n");
    }else{
      NNEW(pf_uf,double,neq[1]);
      NNEW(pf_uref,double,mt**nk);
      isiz=mt**nk;cpypardou(pf_uref,vold,&isiz,&num_cpus);
      pf_neqarm=neq[1];
      pf_taucur=pf_tauv;
      pf_on=1;
      if(getenv("CCX_PATHFOLLOW_CLIP")!=NULL)
        pf_clip=atof(getenv("CCX_PATHFOLLOW_CLIP"));
      if(!(pf_clip>0.)) pf_clip=0.05;
      if(getenv("CCX_PATHFOLLOW_LINCHECK")!=NULL)
        pf_lincheck=atoi(getenv("CCX_PATHFOLLOW_LINCHECK"));
      if(getenv("CCX_PATHFOLLOW_DTHETA")!=NULL)
        pf_dtheta_eng=atof(getenv("CCX_PATHFOLLOW_DTHETA"));
      if(!(pf_dtheta_eng>0.)) pf_dtheta_eng=1.e-3;
      printf("[PATHFOLLOW] armed: tau=%.6e per increment, "
             "|dlambda| clipped at %.3e per iteration.\n",
             pf_tauv,pf_clip);
      printf("[PATHFOLLOW] lambda is decoupled from the step time and "
             "MAY DECREASE; theta stays monotone so dtime>0.\n");
      printf("[PATHFOLLOW] ordinary control until the measured "
             "dissipation of an accepted increment reaches 0.2*tau, "
             "then the constraint takes over.\n");

      /* ---- crack-opening control ---------------------------------
         Build the control functional c once from the REFERENCE geometry
         of the UC6 facets: phi(u)=c^T u is the mean normal separation, a
         linear functional, so dg/du=c and dg/dlambda=0 are exact by
         construction.  See the block comment in pathfollow.c for why this
         replaces the dissipation constraint on a localised cohesive
         crack. */

      /* ---- mixed-mode crack control ------------------------------
         CCX_CRACK_CONTROL=<dphi> arms the generalisation of the above:
         the control coordinate is the effective separation the UC6 law
         itself advances along, deff^2 = max(dn,0)^2+beta*|ds|^2, frozen
         into an affine functional once per attempt.  The two are
         mutually exclusive; CCX_PATHFOLLOW_COD is kept unchanged so that
         the Mode-I result stays a regression test. */

      if((getenv("CCX_CRACK_CONTROL")!=NULL)&&
         (getenv("CCX_PATHFOLLOW_COD")!=NULL)){
        printf("[CRACKCTL] *ERROR: CCX_CRACK_CONTROL and "
               "CCX_PATHFOLLOW_COD both define the control coordinate; "
               "set only one.  Not armed.\n");
      }else if(getenv("CCX_CRACK_CONTROL")!=NULL){
        char *cce;
        ITG ce,ncoh=0;
        for(ce=0;ce<*ne;ce++){
          if((lakon[8*ce]!='U')||(lakon[8*ce+1]!='C')||
             (lakon[8*ce+2]!='6')) continue;
          ncoh++;
        }
        if(ncoh==0){
          printf("[CRACKCTL] *ERROR: CCX_CRACK_CONTROL needs UC6 "
                 "cohesive elements; none found.  Not armed.\n");
        }else if(crackcontrol_selftest()!=0){
          printf("[CRACKCTL] *ERROR: the kinematics self test failed; "
                 "refusing to arm.\n");
        }else{
          pf_dphi=atof(getenv("CCX_CRACK_CONTROL"));
          if(!(pf_dphi>0.)){
            printf("[CRACKCTL] *ERROR: CCX_CRACK_CONTROL must be a "
                   "positive control increment.  Not armed.\n");
          }else{
            /* Default DISS: the process zone restricted to where it is
               LOADING.  Measured on the target, the unrestricted zone
               mean runs backwards while the loading mean advances
               monotonically - see crackcontrol.c. */
            pf_ccmode=2;
            cce=getenv("CCX_CRACK_CONTROL_MODE");
            if(cce!=NULL){
              if((strcmp(cce,"MEAN")==0)||(strcmp(cce,"0")==0)) pf_ccmode=0;
              else if((strcmp(cce,"ZONE")==0)||(strcmp(cce,"1")==0)) pf_ccmode=1;
              else if((strcmp(cce,"DISS")==0)||(strcmp(cce,"2")==0)) pf_ccmode=2;
              else printf("[CRACKCTL] unknown CCX_CRACK_CONTROL_MODE "
                          "\"%s\"; keeping DISS\n",cce);
            }
            cce=getenv("CCX_CRACK_CONTROL_ENGAGE");
            if(cce!=NULL) pf_ccengage=atoi(cce);
            cce=getenv("CCX_CRACK_CONTROL_EPS");
            if(cce!=NULL) pf_eps=atof(cce);
            if(!(pf_eps>0.)) pf_eps=1.e-5;
            cce=getenv("CCX_CRACK_CONTROL_GROW");
            if(cce!=NULL) pf_ccgrow=atof(cce);
            if(!(pf_ccgrow>=1.)) pf_ccgrow=1.1;
            NNEW(pf_cvec,double,neq[1]);
            if(pathfollow_cod_arm(pf_cvec,neq[1])==1){
              pf_codmode=2;
              pf_ccarmed=1;
              pf_dphicur=pf_dphi;
              pf_engaged=0;      /* ordinary control until the crack has
                                    a process zone to control          */
              printf("[CRACKCTL] armed on %" ITGFORMAT " UC6 facet(s): "
                     "mode=%s, control increment dphi=%.6e per "
                     "increment.\n",ncoh,
                     (pf_ccmode==0)?"MEAN (all facets)":
                     ((pf_ccmode==1)?"ZONE (process zone)":
                      "DISS (cohesive dissipation)"),pf_dphi);
              printf("[CRACKCTL] the functional is refrozen from the "
                     "committed state at every attempt and the "
                     "constraint is incremental: g = c^T(u-u_n) - "
                     "dphi.\n");
              if(pf_ccengage>0)
                printf("[CRACKCTL] engagement deferred to increment %"
                       ITGFORMAT ".\n",pf_ccengage);
              else
                printf("[CRACKCTL] engages as soon as the process zone "
                       "is non-empty.\n");
            }else{
              printf("[CRACKCTL] *ERROR: could not arm.\n");
            }
          }
        }
      }

      if(getenv("CCX_PATHFOLLOW_COD")!=NULL){
        ITG ce,ci,ck,cip,cn,cnp,cdof,ncoh=0;
        double ca[3],cb[3],cnv[3],cnorm,csh[3],cw;
        for(ce=0;ce<*ne;ce++){
          if(ipkon[ce]<0) continue;
          if((lakon[8*ce]!='U')||(lakon[8*ce+1]!='C')||
             (lakon[8*ce+2]!='6')) continue;
          ncoh++;
        }
        if(ncoh==0){
          printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW_COD needs UC6 "
                 "cohesive elements; none found.  Not armed.\n");
        }else{
          NNEW(pf_cvec,double,neq[1]);
          cw=1./(3.*(double)ncoh);
          for(ce=0;ce<*ne;ce++){
            if(ipkon[ce]<0) continue;
            if((lakon[8*ce]!='U')||(lakon[8*ce+1]!='C')||
               (lakon[8*ce+2]!='6')) continue;
            for(ck=0;ck<3;ck++){
              ca[ck]=co[3*(kon[ipkon[ce]+1]-1)+ck]
                    -co[3*(kon[ipkon[ce]+0]-1)+ck];
              cb[ck]=co[3*(kon[ipkon[ce]+2]-1)+ck]
                    -co[3*(kon[ipkon[ce]+0]-1)+ck];
            }
            cnv[0]=ca[1]*cb[2]-ca[2]*cb[1];
            cnv[1]=ca[2]*cb[0]-ca[0]*cb[2];
            cnv[2]=ca[0]*cb[1]-ca[1]*cb[0];
            cnorm=sqrt(cnv[0]*cnv[0]+cnv[1]*cnv[1]+cnv[2]*cnv[2]);
            if(!(cnorm>1.e-30)) continue;
            for(ck=0;ck<3;ck++) cnv[ck]/=cnorm;
            for(cip=0;cip<3;cip++){
              for(ci=0;ci<3;ci++) csh[ci]=1./6.;
              csh[cip]=2./3.;
              for(ci=0;ci<3;ci++){
                cn =kon[ipkon[ce]+ci];      /* minus side */
                cnp=kon[ipkon[ce]+ci+3];    /* plus side  */
                for(ck=1;ck<mt;ck++){
                  cdof=nactdof[mt*(cnp-1)+ck];
                  if(cdof>0) pf_cvec[cdof-1]+=cw*csh[ci]*cnv[ck-1];
                  cdof=nactdof[mt*(cn-1)+ck];
                  if(cdof>0) pf_cvec[cdof-1]-=cw*csh[ci]*cnv[ck-1];
                }
              }
            }
          }
          pf_dphi=atof(getenv("CCX_PATHFOLLOW_COD"));
          if(!(pf_dphi>0.)) pf_dphi=5.e-5;
          if(pathfollow_cod_arm(pf_cvec,neq[1])==1){
            pf_codmode=1;
            pf_engaged=1;     /* no warm-up: phi is meaningful at once */
            printf("[PATHFOLLOW] crack-opening control armed on %"
                   ITGFORMAT " UC6 facet(s); mean normal separation "
                   "advances by %.6e per increment.\n",ncoh,pf_dphi);
            {ITG cnz=0; double cnn=0.;
             for(ce=0;ce<neq[1];ce++){
               cnn+=pf_cvec[ce]*pf_cvec[ce];
               if(pf_cvec[ce]!=0.) cnz++;}
             printf("[PATHFOLLOW] control functional: |c|=%.6e, %"
                    ITGFORMAT " non-zero of %" ITGFORMAT " equations\n",
                    sqrt(cnn),cnz,neq[1]);}
            printf("[PATHFOLLOW] the control is monotone through a "
                   "snap-back by construction; lambda is free.\n");
          }
        }
      }
    }
  }

  /* backup fields needed to roll back element deletion after
     an unsuccessful damage re-equilibration */

  if((*ndmat_>0)&&(*iexpl<=1)){
    NNEW(ipkondamageini,ITG,ne0);
    NNEW(damdamageini,double,mi[0]*ne0);
    NNEW(damde1prev,double,mi[0]*ne0);
    if(damage_de12_enabled){
      /* UNSYM stage 2 buffer: 6 effective-stress plus 6 dD/d(eps)
         components per integration point, filled by resultsmech.f and
         consumed by mafilldamas.f */

      if(damage_tangent_mode==2){
        NNEW(damage_damjac,double,12*mi[0]*ne0);
      }

      /* R4 viscous regularisation buffers.  damage_damvisc is the trial
         relaxed degradation rebuilt from damage_damviscini on every
         results() call; damage_damviscini is the committed value, taken
         over at the start of the next physical increment exactly like
         damdamageini.  A rejected increment simply discards the trial. */

      damage_visc_env=getenv("CCX_DAMAGE_VISCOSITY");
      if(damage_visc_env!=NULL) damage_visc_eta=atof(damage_visc_env);
      if(damage_visc_eta<0.) damage_visc_eta=0.;
      if(damage_visc_eta>0.){
        NNEW(damage_damvisc,double,mi[0]*ne0);
        NNEW(damage_damviscini,double,mi[0]*ne0);
      }
      NNEW(damage_de13_trigger_value,double,ne0);
      NNEW(damage_de13_trigger_ip,ITG,ne0);
      for(i=0;i<ne0;i++){
        damage_de13_trigger_value[i]=-1.;
        damage_de13_trigger_ip[i]=0;
      }
    }
    if(*nmethod!=4) NNEW(veolddamageini,double,mt**nk);

    /* committed hard-deletion history; one file per CalculiX job */

    strcpy2(damagefilename,jobnamec,132);
    strcat(damagefilename,".damage");

    if(*istep==1){
      fdamage=fopen(damagefilename,"w");
    }else{
      fdamage=fopen(damagefilename,"a");
    }

    if(fdamage==NULL){
      printf(" *ERROR in nonlingeo: cannot open %s for writing\n",
             damagefilename);
      FORTRAN(stop,());
    }

    printf("[DAMAGE PATCH A3] adaptive-fast transactional damage enabled "
           "(batch_max=%d floor=%.2f direct_alpha=%.2f)\n",
           DAMAGE_FAST_BATCH_MAX,DAMAGE_FAST_MIN_FRACTION,
           DAMAGE_FAST_DIRECT_ALPHA);
    if(damage_de12_enabled){
      printf("[DAMAGE DE1.2] Newton-integrated trial damage enabled "
             "(%" ITGFORMAT " material(s)); BK1 one-pass local update / "
             "one global Newton solve; exact cell VTK enabled\n",
             damage_de12_matcount);
      if(damage_dm20_matcount>0){
        printf("[DAMAGE DM2.0] tabulated ductile initiation enabled "
               "(%" ITGFORMAT " material(s)); eps_f(eta) linear interpolation\n",
               damage_dm20_matcount);
      }
      printf("[DAMAGE DE1.3.1] terminal failure enabled "
             "(Ddelete=%.4f batch_max=%d); transactional A3 topology "
             "backend\n",damage_de13_delete_d,DAMAGE_DE13_BATCH_MAX);
      printf("[DAMAGE SOLVER NC1] adaptive slow-Newton controller enabled "
             "(stock_ic=%" ITGFORMAT " cap=%d max_extra=%d); "
             "final stock convergence criteria unchanged\n",
              (ITG)icref,DAMAGE_SLOW_NEWTON_MAX_ITERS,
              DAMAGE_SLOW_NEWTON_MAX_EXTRA);
      /* crack-band check: L*sigma_y/(u_f*E) must stay below 1 or the
         element softens faster than it can elastically unload */

      FORTRAN(damsnapcheck,(co,kon,ipkon,lakon,&ne0,ielmat,mi,ndmcon,
                            dmcon,ndmat_,ntmat_,elcon,ncmat_,plicon,
                            nplicon,npmat_,&damage_snap_ratio,
                            &damage_snap_elem,&damage_snap_bad));
      printf("[DAMAGE CRACK-BAND] max L*sigma_y/(u_f*E) = %.4f in "
             "element %" ITGFORMAT "\n",damage_snap_ratio,
             damage_snap_elem);
      if(damage_snap_bad>0){
        printf("[DAMAGE CRACK-BAND] *WARNING: %" ITGFORMAT " element(s) "
               "have r>1, i.e. the element-level softening branch is "
               "steeper than its own elastic unloading.  This is an "
               "indicator, not a proof: the structure may still find "
               "equilibrium by redistributing into neighbours.  It "
               "does mean convergence with a consistent tangent can "
               "be poor; refining the mesh or raising u_f removes "
               "it.\n",damage_snap_bad);
      }

      if(damage_tangent_mode==1){
        printf("[DAMAGE TANGENT BK2] FD_SYM enabled; finite-difference "
               "dD/deps with symmetric damage Jacobian correction\n");
      }else if(damage_tangent_mode==2){
        printf("[DAMAGE TANGENT UNSYM] stage 1: asymmetric assembly and "
               "solve with the unchanged symmetric g(D)*Cep tangent; "
               "results must match the symmetric path exactly\n");
      }else{
        printf("[DAMAGE TANGENT BK2] secant g(D)*Cep baseline enabled; "
               "set CCX_DAMAGE_TANGENT=FD_SYM for tangent trial\n");
      }
      if(damage_reeq_scale_mode==1){
        printf("[DAMAGE SOLVER NC2] terminal same-load correction scale "
               "uses the converged physical-increment displacement norm; "
               "stock residual/correction tolerances retained\n");
      }else{
        printf("[DAMAGE SOLVER NC2] stock zero-load re-equilibration scale; "
               "set CCX_DAMAGE_REEQ_SCALE=PHYSICAL for NC2 trial\n");
      }
      if(damage_reeq_uam_floor>0.){
        printf("[DAMAGE REEQ FLOOR] DIAGNOSTIC: the re-equilibration "
               "displacement reference is held at or above %.3e of its "
               "running maximum, so the RELATIVE correction criterion "
               "cannot collapse as the step is cut.  The perturbation of a "
               "same-load topology solve does not scale with dt, so cutting "
               "the step cannot converge it (J-09).  This CHANGES THE "
               "CONVERGENCE CRITERION and therefore the answer; going "
               "further with it is not by itself a success\n",
               damage_reeq_uam_floor);
      }
      if(damage_linesearch_mode==1){
        printf("[DAMAGE SOLVER BK3] adaptive line search enabled "
               "(growth=%.2f lambda=[%.2f,%.2f] trials=%d); full Newton retained "
               "while residual contracts\n",DAMAGE_LINESEARCH_GROWTH,
               DAMAGE_LINESEARCH_MIN,DAMAGE_LINESEARCH_MAX,
               DAMAGE_LINESEARCH_MAX_TRIALS);
      /* The damage line search is clamped to a compiled-in floor and trial
         cap.  On m12_field_soft50 both bind SIMULTANEOUSLY on every iteration
         of the failing increment - lambda pinned at exactly 0.100 with
         trials=3, seven iterations running, residual oscillating in
         1.70e-3..2.30e-3 with no contraction and a stable active set (E-81).
         A search that asks for a shorter step on every iteration and is
         refused on every iteration is a search whose floor is the binding
         constraint, not a search that has converged.

         The structural FD probe says the tangent there is NOT the problem:
         2.52e-3 median at the wall against a 1.33e-4 elastic noise floor on
         the same deck and SP1's healthy 1.92e-3 (E-91).  So make the two
         clamps measurable instead of assumed.  Defaults are the compiled-in
         values, so an unset environment is bit-identical. */
      if((damage_de13_env=getenv("CCX_DAMAGE_LS_MIN"))!=NULL){
        damage_ls_min=atof(damage_de13_env);
        if(damage_ls_min<1.e-6) damage_ls_min=1.e-6;
        if(damage_ls_min>DAMAGE_LINESEARCH_MAX) damage_ls_min=DAMAGE_LINESEARCH_MAX;
      }
      if((damage_de13_env=getenv("CCX_DAMAGE_LS_TRIALS"))!=NULL){
        damage_ls_trials=atoi(damage_de13_env);
        if(damage_ls_trials<1) damage_ls_trials=1;
        if(damage_ls_trials>32) damage_ls_trials=32;
      }

      /* ---- THE BACKTRACKING LADDER ---------------------------------
         Measured on the target with CCX_DAMAGE_LS_PROBE, over the 245
         line-search activations of a run to the recorded wall:

           in 114 of them (46.5%) the best step length lies BELOW the
           floor of 0.1, i.e. outside what the ladder {1.0, 0.5, 0.1}
           can reach;

           at the wall increment 353, iterations 12 to 16, the floor step
           RAISES |R|inf - 2.852e-03, 3.053e-03, 3.353e-03, 3.649e-03,
           3.822e-03 - while alpha=0.03 lowers it at every one of them.
           The iteration was converging until iteration 11 and the ladder
           turned it round.

         So the ladder is deepened (floor 1e-3, up to 8 halvings) and the
         search now falls back on the BEST alpha it measured rather than
         the last.  Neither changes any convergence criterion: what an
         increment must satisfy to be accepted is untouched, this only
         chooses how far along a direction to step.

         CCX_DAMAGE_LS_LEGACY=1 restores the old ladder exactly, which is
         what makes a controlled A/B possible from ONE binary. */

      /* The ladder is the thing that was wrong, so prove it is right
         before using it, on every run, the way pathfollow and
         crackcontrol do.  A failure here is not something to carry on
         through: the search would silently go back to handing back a
         step it had measured to be worse. */

      if(lsladder_selftest()!=0){
        printf("[DAMAGE LINESEARCH] *ERROR: the backtracking ladder self "
               "test failed; disabling the adaptive line search rather "
               "than running with a rule that is not the one that was "
               "tested.\n");
        damage_linesearch_mode=0;
      }

      if(getenv("CCX_DAMAGE_LS_LEGACY")!=NULL){
        damage_ls_legacy=1;
        damage_ls_min=DAMAGE_LINESEARCH_MIN;
        damage_ls_trials=DAMAGE_LINESEARCH_MAX_TRIALS;
        printf("[DAMAGE LINESEARCH] LEGACY ladder: floor %.2f, %d trials, "
               "and the last trial is taken whether or not it contracts.\n",
               DAMAGE_LINESEARCH_MIN,DAMAGE_LINESEARCH_MAX_TRIALS);
      }else{
        if(getenv("CCX_DAMAGE_LS_MIN")==NULL) damage_ls_min=1.e-3;
        if(getenv("CCX_DAMAGE_LS_TRIALS")==NULL) damage_ls_trials=8;
        printf("[DAMAGE LINESEARCH] ladder: floor %.4e, up to %"
               ITGFORMAT " trials, and the BEST alpha measured is taken "
               "when none contracts.\n",damage_ls_min,damage_ls_trials);
      }
      if((damage_ls_min!=DAMAGE_LINESEARCH_MIN)||
         (damage_ls_trials!=DAMAGE_LINESEARCH_MAX_TRIALS)){
        printf("[DAMAGE LINESEARCH] floor %.4e (default %.2f), "
               "max trials %" ITGFORMAT " (default %d)\n",
               damage_ls_min,DAMAGE_LINESEARCH_MIN,
               damage_ls_trials,DAMAGE_LINESEARCH_MAX_TRIALS);
      }

      }else{
        printf("[DAMAGE SOLVER BK3] disabled; set "
               "CCX_DAMAGE_LINESEARCH=ADAPTIVE for damage globalization\n");
      }
      if(damage_topology_deferred_mode==1){
        printf("[DAMAGE TOPOLOGY BK4] deferred sparse-structure compaction "
               "enabled; terminal elements leave the assembly immediately, "
               "remastruct is deferred until an active node becomes orphaned\n");
      }else{
        printf("[DAMAGE TOPOLOGY BK4] immediate remastruct baseline; set "
               "CCX_DAMAGE_TOPOLOGY=DEFERRED for deferred compaction\n");
      }
    }else{
      printf("[DAMAGE DE1.1] fixed-point controller available for legacy "
             "DE1 states; exact cell VTK enabled\n");
    }
    fflush(stdout);

    if(*istep==1){
      fprintf(fdamage,"# CalculiX committed damage deletion history v2\n");
      fprintf(fdamage,"# element step increment step_time total_time "
              "material damage critical_ip batch\n");
      fprintf(fdamage,"# DE1.3.1: damage is the captured terminal-trigger D "
              "at first topology marking\n");
      fflush(fdamage);
    }
  }
  
  /*********************************************************************/
  
  /* calculating of the acceleration due to force discontinuities
     (external - internal force) at the start of a step */
  
  /*********************************************************************/
  
  if((*nmethod==4)&&(*ithermal!=2)&&(icfd==0)){
    bet=(1.-alpha[0])*(1.-alpha[0])/4.;
    gam=0.5-alpha[0];
      
    /* calculating the stiffness and mass matrix */
      
    reltime=0.;
    time=0.;
    dtime=0.;
      
    FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
		      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
		      xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,nam,ampli,
		      &time,&reltime,ttime,&dtime,ithermal,nmethod,xbounold,
		      xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,nodeforc,
		      ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,
		      ilboun,nelemload,sideload,mi,ntrans,trab,inotr,veold,
		      integerglob,doubleglob,tieset,istartset,iendset,ialset,
		      ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,ipobody,
		      iponoeln,inoeln,ipkon,kon,ielprop,prop,ielmat,shcon,nshcon,
		      rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));
      
    time=0.;
    dtime=1.;
    
    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's)
	if contact arises the number of MPC's can also change */
      
    cam[0]=0.;cam[1]=0.;cam[2]=0.;
      
    if(icascade==2){
      memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
      RENEW(nodempc,ITG,3*memmpcref_);
      for(k=0;k<3*memmpcref_;k++){
	nodempc[k]=nodempcref[k];}
      RENEW(coefmpc,double,memmpcref_);
      for(k=0;k<memmpcref_;k++){
	coefmpc[k]=coefmpcref[k];}
    }

    newstep=0;
    FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
		       nmpc,ikboun,ilboun,nboun,xbounold,aux,iaux,
		       &maxlenmpc,ikmpc,ilmpc,&icascade,
		       kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
		       &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi,&kchdep));
    if(icascade==2){
      for(k=0;k<3*memmpc_;k++){
	nodempcref[k]=nodempc[k];}
      for(k=0;k<memmpc_;k++){
	coefmpcref[k]=coefmpc[k];}
    }

    /* recalculating the matrix structure */
    
    if(icascade>0){
      remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		 &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		 labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		 kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		 neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		 &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar,
		 typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
		 nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus,
		 ielmat,matname);
    }

    /* invert nactdof */

    NNEW(nactdofinv,ITG,1);
      
    iout=-1;
    ielas=1;
      
    MNEW(fn,double,mt**nk);
    NNEW(stx,double,6*mi[0]**ne);
      
    if((*iexpl<=1)||(*mortar==-1)){intscheme=1;}
      
    if(ne1d2d==1)NNEW(inum,ITG,*nk);
    results(co,nk,kon,ipkon,lakon,ne,vold,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1old,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xbounold,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,&bet,
	    &gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
	    ncmat_,nstate_,sti,vini,ikboun,ilboun,ener,enern,emeini,xstaten,
	    eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
	    ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,
	    &ne0,thicke,shcon,nshcon,
	    sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
	    inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	    itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	    islavquadel,aut,irowt,jqt,&mortartrafoflag,
	    &intscheme,physcon,dam,damn,iponoel);
      
    SFREE(fn);SFREE(stx);if(ne1d2d==1)SFREE(inum);
      
    if(*mortar<2){
      iout=0;
      ielas=0;
	  
      reltime=0.;
      time=0.;
      dtime=0.;
    }
      
    if(*iexpl>1){

      mscalmethod=0;
      nloadrhs=*nload;nbodyrhs=*nbody;
	  
      /* Explicit: Calculation of stable time increment according to
	 Courant's Law  Carlo Monjaraz Tec (CMT) and Selctive Mass Scaling CC*/

      /*Mass Scaling
	mscalmethod < 0: no explicit dynamics
	mscalmethod = 0: no mass scaling
	mscalmethod = 1: selective mass scaling for nonlinearity after 
	Olovsson et. al 2005

        mscalmethod=2 and mscalmethod=3 correspond to 0 and 1, 
        respectively with in addition contact scaling active; contact
        scaling is activated if the user time increment cannot be satisfied */

      dtset=*tmin*(*tper);
      NNEW(smscale,double,*ne);
	  
      FORTRAN(calcstabletimeincvol,(&ne0,elcon,nelcon,rhcon,nrhcon,alcon,
				    nalcon,orab,ntmat_,ithermal,alzero,plicon,
				    nplicon,plkcon,nplkcon,npmat_,mi,&dtime,
				    xstiff,ncmat_,vold,ielmat,t0,t1,matname,
				    lakon,wavespeed,nmat,ipkon,co,kon,&dtvol,
				    alpha,smscale,&dtset,&mscalmethod,mortar,
				    jobnamef,iperturb));

      printf(" Explicit time integration: Volumetric COURANT initial stable time increment:%e\n\n",dtvol);

      if(dtvol>(*tmax*(*tper))){
	*tinc=*tmax*(*tper);}
      else if(dtvol<dtset){
	*tinc=dtset;}
      else{
	*tinc=dtvol;
      }
	  
      dtheta=(*tinc)/(*tper);
      dthetaref=dtheta;
      printf(" SELECTED time increment (not considering penalty contact):%e\n\n",*tinc);
    }
      
    /* in mafillsm the stiffness and mass matrix are computed;
       The primary aim is to calculate the mass matrix (not 
       lumped for an implicit dynamic calculation, lumped for an
       explicit dynamic calculation). However:
       - for an implicit calculation the mass matrix is "doped" with
       a small amount of stiffness matrix, therefore the calculation
       of the stiffness matrix is needed.
       - for an explicit calculation the stiffness matrix is not 
       needed at all. Since the calculation of the mass matrix alone
       is not possible in mafillsm, the determination of the stiffness
       matrix is taken as unavoidable "ballast". */
      
    NNEW(ad,double,neq[1]);
    NNEW(au,double,nzs[1]);

    mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
		 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		 nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		 nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		 nmethod,ikmpc,ilmpc,ikboun,ilboun,
		 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		 ielmat,ielorien,norien,orab,ntmat_,
		 t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		 xstiff,npmat_,&dtime,matname,mi,
		 ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		 physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		 &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		 xstateini,xstate,thicke,integerglob,doubleglob,
		 tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		 pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
		 iponoeln,inoeln,network,ntrans,inotr,trab,smscale,&mscalmethod,
		 set,nset,islavquadel,aut,irowt,jqt,&mortartrafoflag);
    
    if(*nmethod==0){
	  
      /* error occurred in mafill: storing the geometry in frd format */
	  
      ++*kode;
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	MNEW(ipneigh,ITG,*nk);
      }
	  
      ptime=*ttime+time;
      frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	  kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	  ielprop,prop,sti,damn,&errn);
	  
      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}      
#ifdef COMPANY
      FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif	  
      if(nmethodold==0){FORTRAN(stopwithout201,());}else{FORTRAN(stop,());}
	  
    }

    /* massless contact: setting up the system matrices based on
       the stiffness and mass matrix; these matrices are not
       assumed to change during the step;
       factorization of the LHS matrix */
    
    if(*mortar==-1){
      if(ncont!=0){
	nmasts=nmastnode[*ntie];

	NNEW(kslav,ITG,3**nslavs);
	NNEW(lslav,ITG,3**nslavs);
	NNEW(ktot,ITG,3**nslavs+3*nmasts);
	NNEW(ltot,ITG,3**nslavs+3*nmasts);
	NNEW(fric,double,*nslavs);

	/*  Create set of slave and slave+master contact DOFS (sorted);
	    assign a friction coefficient to each slave node */
      
	FORTRAN(create_contactdofs,(kslav,lslav,ktot,ltot,nslavs,islavnode,
				    &nmasts,imastnode,nactdof,mi,&neqtot,
				    nslavnode,fric,tieset,tietol,ntie,elcon,
				    ncmat_,ntmat_));
      }else{
	neqtot=0;
      }
      
      /*   RENEW(kslav,ITG,3**nslavs);
      RENEW(lslav,ITG,3**nslavs);
      RENEW(ktot,ITG,neqtot);
      RENEW(ltot,ITG,neqtot);*/

      /* create RHS of system:  M/dt - (aM + bK)/2 and store in adc,auc */
      
      NNEW(adc,double,neq[0]);
      for(k=0;k<neq[0];k++){
	adc[k]=adb[k]/(*tinc)-(alpham*adb[k]+betam*ad[k])/2.0;
      }

      NNEW(auc,double,nzs[0]);
      for(k=0;k<nzs[0];k++){
	auc[k]=aub[k]/(*tinc)-(alpham*aub[k]+betam*au[k])/2.0;
      }

      /* create LHS of system:  M/dt + (aM + bK)/2  and store in adb,aub */

      for(k=0;k<neq[0];k++){
	adb[k]=adb[k]/(*tinc)+(alpham*adb[k]+betam*ad[k])/2.0;
      }

      for(k=0;k<nzs[0];k++){
	aub[k]=aub[k]/(*tinc)+(alpham*aub[k]+betam*au[k])/2.0;
      }

      /* reduce LHS and RHS by removing contact dofs (diagonal terms
         are set to 1, off-diagonal terms to 0 */

      if(ncont!=0){
	FORTRAN(reducematrix,(aub,adb,jq,irow,neq,&neqtot,ktot));
	FORTRAN(reducematrix,(auc,adc,jq,irow,neq,&neqtot,ktot));
      }

      /* factorize the LHS */

      if(*isolver==0){
#ifdef SPOOLES

	spooles_factor(adb,aub,adb,aub,&sigma,icol,irow,
		       &neq[0],&nzs[0],&symmetryflag,&inputformat,&nzs[0]);

#else
	printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],token);
#else
	printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_factor(adb,&aub,adb,aub,&sigma,icol,&irow,&neq[0],&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
		       &symmetryflag,&inputformat,jq,&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_factor_main(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			   &symmetryflag,&inputformat,jq,&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }

      // Storing contact force vector initial solution

      if(ncont!=0){
	NNEW(aloc,double,3**nslavs);
	NNEW(alglob,double,neqtot);
      }

      /* no nlgeom and no nonlinear material for massless explicit dynamics */
      
      if((iperturb[0]<3)&&(iperturb[1]==0)) masslesslinear=1;

      /* check whether the output consists of displacements only */

      FORTRAN(checkdispoutonly,(prlab,nprint,nlabel,filab,&idispfrdonly));

      if(idispfrdonly==1){
	NNEW(inumcp,ITG,*nk);
	strcpy1(&cflag[0],&filab[4],1);
	FORTRAN(createinum,(ipkon,inumcp,kon,lakon,nk,ne,&cflag[0],nelemload,
			    nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,
			    ielmat,ielprop,prop));
      }

      
    } //endif massless
      
    /* mass x acceleration = f(external)-f(internal) 
       only for the mechanical loading*/
      
    /* not needed for massless contact */
    
    if(*mortar!=-1){
      for(k=0;k<neq[0];++k){b[k]=fext[k]-f[k];}
    }
      
    if(*iexpl<=1){
	  
      /* a small amount of stiffness is added to the mass matrix
	 otherwise the system leads to huge accelerations in 
	 case of discontinuous load changes at the start of the step */
	  
      dtime=*tinc/10.;
      scal1=bet*dtime*dtime*(1.+alpha[0]);
      for(k=0;k<neq[0];++k){
	ad[k]=adb[k]+scal1*ad[k];
      }
      for(k=0;k<nzs[0];++k){
	au[k]=aub[k]+scal1*au[k];
      }
      if(*isolver==0){
#ifdef SPOOLES
	spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		&symmetryflag,&inputformat,&nzs[2]);
#else
	printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
#endif
      }
      else if((*isolver==2)||(*isolver==3)){
	preiter(ad,&au,b,&icol,&irow,&neq[0],&nzs[0],isolver,iperturb);
      }
      else if(*isolver==4){
#ifdef SGI
	token=1;
	sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
#else
	printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[0],&nzs[0]);
#else
	printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		     &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
	printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
	pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		    &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
	printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	FORTRAN(stop,());
#endif
      }
    }
      
    else{

      /* explicit dynamics; no selective mass scaling
         (at most spring scaling) */
      
      /* if massless contact: no acceleration needed */
      
      if((mscalmethod==0)||(mscalmethod==2)){
        if(*mortar!=-1){
	  for(k=0;k<neq[0];++k){b[k]=(fext[k]-f[k])/adb[k];}
	}
      }
	  
      else{

	/* explicit dynamics with selective mass scaling */

	inputformat=0;
	if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,&nzs[2]);
	  spooles_solve(b,&neq[0]);
#else
	  printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	  token=1;
	  sgi_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],token);
	  sgi_solve(b,token);
#else
	  printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	  tau_factor(adb,&aub,adb,aub,&sigma,icol,&irow,&neq[0],&nzs[0]);
	  tau_solve(b,&neq[0]);
#else
	  printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,jq,&nzs[0]);

	  pardiso_solve(b,&neq[0],&symmetryflag,&inputformat,&nrhs);
#else
	  printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
	  pastix_factor_main(adb,aub,adb,aub,&sigma,icol,irow,&neq[0],&nzs[0],
			&symmetryflag,&inputformat,jq,&nzs[0]);

	  pastix_solve(b,&neq[0],&symmetryflag,&nrhs);
#else
	  printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
      }
    }
      
    /* for thermal loading the acceleration is set to zero */
      
    for(k=neq[0];k<neq[1];++k){
      b[k]=0.;
    }
      
    /* calculating the displacements, stresses and forces */
      
    if(*mortar!=-1){
      NNEW(v,double,mt**nk);
      isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
      
      NNEW(stx,double,6*mi[0]**ne);
      MNEW(fn,double,mt**nk);
      
      /* setting a "special" time consisting of the first primes;
	 used to recognize the initial acceleration procedure
	 in file resultsini.f */

      if(ne1d2d==1)NNEW(inum,ITG,*nk);
      dtime=1.235711130e-20;
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	      islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
	      inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	      itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	      islavquadel,aut,irowt,jqt,&mortartrafoflag,
	      &intscheme,physcon,dam,damn,iponoel);
      if(ne1d2d==1)SFREE(inum);
      dtime=0.;

      isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
      if(*ithermal!=2){
	isiz=6*mi[0]*ne0;	    
	cpypardou(sti,stx,&isiz,&num_cpus);
      }

      SFREE(v);SFREE(stx);SFREE(fn);
    }
    SFREE(ad);SFREE(au);
      
    /* the mass matrix is kept for subsequent calculations, therefore,
       no new mass calculation is necessary for the remaining iterations
       in the present step */
      
    mass[0]=0;intscheme=0;
    energyref=energy[0]+energy[1]+energy[2]+energy[3];

    if(*iexpl<=1){
	  
      NNEW(tmp,double,neq[1]);
      NNEW(adblump,double,neq[1]);
      for(k=0;k<neq[1];k++){
	tmp[k] = 1;
      }
      if(nasym==0){
	opmain(&neq[1],tmp,adblump,adb,aub,jq,irow); 
      }else{
	FORTRAN(opas,(&neq[1],tmp,adblump,adb,aub,jq,irow,nzs)); 
      }
      SFREE(tmp);
    }
  }

  /* warning: for C3D8R-elements the stiffness is needed in subroutine
              hgforce, therefore, in explicit dynamic steps
              with C3D8R-elements icmd should not be set to 3 */
  
  if(*iexpl>1) icmd=3;

  
  /**************************************************************/
  /* starting the loop over the increments                      */
  /**************************************************************/
  
  newstep=1;
	  
  //    MPADD start
  if((*nmethod==4)&&(*ithermal<2)&&(*iexpl<=1)){
    neini=*ne;
    for(k=0;k<4;k++){
      energystartstep[k]=energy[k];
    }
    emax=0.1*energyref;
    // Anti-stick at the beginning of simulation
  } 
  //    MPADD end

  /* saving the distributed loads (volume heating will be
     added because of friction heating) */

  if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
    nloadref=*nload;
    NNEW(nelemloadref,ITG,2**nload);
    if(*nam>0) NNEW(iamloadref,ITG,2**nload);
    NNEW(sideloadref,char,20**nload);
      
    isiz=2**nload;cpyparitg(nelemloadref,nelemload,&isiz,&num_cpus);
    if(*nam>0){
      isiz=2**nload;cpyparitg(iamloadref,iamload,&isiz,&num_cpus);
    }
    memcpy(&sideloadref[0],&sideload[0],sizeof(char)*20**nload);
  }
  
  while(((1.-theta>1.e-6)||(negpres==1))&&
        (damage_fracture_complete==0)){
      
    if((icutb==0)&&(idamagereeq==0)){

      /* ---- [STALLSTATE] one pass per accepted increment -------------
         Ask whether the run is still making progress, and if it has stopped,
         arm the wall diagnostics by setting the load factor they were
         waiting to be told.  Re-using damage_wall_theta rather than adding a
         second arming path means the probes, their ordering and the reason
         the gate is on theta instead of an increment index all stay exactly
         as they were - the only thing that changes is who supplies the
         number.  Diagnostics only: nothing here is on a solution path. */
      if(damage_stall_ok&&(iinc>0)){
        ITG sv=stallstate_judge(&damage_stall,dtheta);
        if(sv==3){
          printf("[STALLSTATE] inc=%" ITGFORMAT " theta=%.7f dtheta is %.3e of "
                 "this deck's median pace, %" ITGFORMAT " of %" ITGFORMAT
                 " consecutive; not a stall until confirmed\n",
                 iinc,theta,stallstate_ratio(&damage_stall,dtheta),
                 damage_stall.ndisc,damage_stall.nconfirm);
          fflush(stdout);
        }else if(sv==1){
          printf("\n[STALLSTATE] STALLED at inc=%" ITGFORMAT " theta=%.7f: "
                 "dtheta has been below %.3g of this deck's own median pace "
                 "for %" ITGFORMAT " consecutive accepted increments "
                 "(ratio %.3e).\n",
                 iinc,theta,damage_stall.ratio,damage_stall.nconfirm,
                 damage_stall.rseen);
          if(damage_wall_theta<0.){
            damage_wall_theta=theta;
            printf("                   arming the wall diagnostics here, at "
                   "theta=%.9e, without having been told it in advance.\n",
                   damage_wall_theta);
          }else{
            printf("                   the wall diagnostics were already "
                   "armed at theta=%.9e.\n",damage_wall_theta);
          }
          printf("\n");
          fflush(stdout);
        }
        stallstate_note(&damage_stall,dtheta);
      }

      /* ---- [DAMAGE CORR] one pass per converged increment ---------- */
      if(damage_corr_mode==1){
        if((damage_rec_used_in_inc==0)&&(damage_reg_on==0)){
          /* The crawl breaker needs a reference from HEALTHY operation.
             Taking it from the last clean increment before entry is
             wrong: that increment sits at the wall, where the stock
             controller has already cut dtheta hard.  Measured on
             run_c3_s3rad the reference came out 6.592e-06 while the run
             had been healthy at 5e-5..1.2e-4, so the mean inside the
             corridor (1.211e-05) was ABOVE it and the breaker could
             never fire.  Keep a ring of the last clean increments and
             use their MEDIAN instead. */
          damage_dtheta_healthy=dtheta;
          damage_dth_ring[damage_dth_i]=dtheta;
          damage_dth_i=(damage_dth_i+1)%20;
          if(damage_dth_n<20) damage_dth_n++;
        }
        if(damage_corr_on==1){
          ITG cbrk=0;
          double cmean;
          damage_corr_ninc++;
          damage_corr_nsince++;
          damage_corr_nwallstab=0;
          damage_corr_dtsum+=dtheta;damage_corr_dtn++;
          damage_corr_nfact+=damage_reg_napply;
          /* A probe stays a PROBE until it has survived stableneed
             converged increments.  Clearing the flag after the FIRST
             success - the earlier draft did - makes a failure on the
             second or third increment look like a wall on the HELD
             lambda, so the wall handler would ESCALATE away from the
             proven value instead of falling back to it.  The flag and
             the promotion of lamstable must therefore happen together,
             and only once nsince has been reached. */
          if((damage_corr_trial==1)&&
             (damage_corr_nsince>=damage_corr_stableneed)){
            damage_corr_lamstable=damage_corr_lam;
            damage_corr_trial=0;
            printf("[DAMAGE CORR] lambda=%.3e has held %" ITGFORMAT
                   " increments and is now the PROVEN value%s",
                   damage_corr_lam,damage_corr_nsince,"\n");
          }
          if(damage_corr_lam<=0.) damage_corr_clean++;
          else damage_corr_clean=0;
          cmean=(damage_corr_dtn>0)?damage_corr_dtsum/damage_corr_dtn:0.;
          if(damage_corr_clean>=damage_corr_exit){
            printf("[DAMAGE CORR] EXIT: %" ITGFORMAT " consecutive "
                   "increments converged with NO help.  Corridor: %"
                   ITGFORMAT " increments, %" ITGFORMAT " intervention(s), "
                   "%" ITGFORMAT " regularized factorisation(s), step time "
                   "gained %.6e, mean dtime inside %.6e against reference "
                   "%.6e.  Stock Newton has the run back.%s",
                   damage_corr_clean,damage_corr_ninc,damage_corr_nint,
                   damage_corr_nfact,(theta-damage_corr_theta0)**tper,
                   cmean**tper,damage_corr_dtref**tper,"\n");
            fflush(stdout);
            damage_corr_on=0;damage_reg_on=0;damage_corr_lam=0.;
          }else{
            if(damage_corr_ninc>damage_corr_maxinc) cbrk=1;
            if((cbrk==0)&&(damage_corr_ninc>damage_corr_grace)&&
               (damage_corr_dtref>0.)&&
               (cmean<damage_corr_minfrac*damage_corr_dtref)) cbrk=2;
            if(cbrk>0){
              printf("[DAMAGE CORR] CIRCUIT BREAKER (%s): %" ITGFORMAT
                     " increments, %" ITGFORMAT " intervention(s), %"
                     ITGFORMAT " regularized factorisation(s), step time "
                     "gained %.6e, mean dtime %.6e against reference %.6e."
                     "  Carried, not recovering: the corridor closes, the "
                     "mechanism DISARMS, and the next wall goes to the "
                     "original stock stop.%s",
                     (cbrk==1)?"length":"crawl",damage_corr_ninc,
                     damage_corr_nint,damage_corr_nfact,
                     (theta-damage_corr_theta0)**tper,cmean**tper,
                     damage_corr_dtref**tper,"\n");
              fflush(stdout);
              damage_corr_on=0;damage_reg_on=0;damage_corr_lam=0.;
              damage_rec_disarmed=1;
            }else{
              damage_corr_try--;
              if((damage_corr_try<=0)&&(damage_corr_lam>0.)&&
                 (damage_corr_trial==0)&&
                 (damage_corr_nsince>=damage_corr_stableneed)){
                damage_corr_try=damage_corr_tryevery;
                damage_corr_lam=damage_corr_lamstable*0.25;
                if(damage_corr_lam<1.e-4) damage_corr_lam=0.;
                damage_corr_trial=1;
                damage_corr_nsince=0;
                printf("[DAMAGE CORR] probing weaker help: lambda %.3e -> "
                       "%.3e (proven value kept for fallback)%s",
                       damage_corr_lamstable,damage_corr_lam,"\n");
                fflush(stdout);
              }
              damage_reg_lambda=damage_corr_lam;
              damage_reg_on=(damage_corr_lam>0.)?1:0;
            }
          }
          damage_reg_napply=0;
        }
      }

      /* [DAMAGE RESCUE] the previous increment converged.  If it was the
         rescue attempt, it is now committed together with any re-equilibration
         it needed: switch backtracking off at once and re-arm for the next
         INDEPENDENT wall.  theta has advanced, so this cannot loop. */

      /* [DAMAGE CT] acceptance needs BOTH the stock residual criteria and
         the constraint.  checkconvergence has already accepted on R; the
         constraint is verified here, before anything is committed. */
      if((damage_ct_on==1)&&(fabs(damage_ct_cprev)>damage_ct_tolc)){
        printf("[DAMAGE CT] increment accepted on R but the CONSTRAINT is "
               "not satisfied: |c|=%.6e > tol_c=%.6e.  Requesting the clean "
               "PARTIAL exit rather than committing a state that does not "
               "lie on the constraint.%s",
               fabs(damage_ct_cprev),damage_ct_tolc,"\n");
        fflush(stdout);
        damage_ct_partial=1;damage_ct_on=0;
      }
      /* [DAMAGE CT] lambda >= 1 is NOT a completion in this MVP. */
      if((damage_ct_on==1)&&(damage_ct_lam>=1.)){
        printf("[DAMAGE CT] lambda has reached %.12e >= 1.  The MVP has no "
               "terminal landing, so this is a PARTIAL stop, not a "
               "completed step.%s",damage_ct_lam,"\n");
        fflush(stdout);
        damage_ct_partial=1;damage_ct_on=0;
      }
      if(damage_rec_used_in_inc==0){
        damage_rec_healthy++;
        /* [DAMAGE TR] five consecutive clean increments, derived from
           INCREMENT NUMBERS rather than from a running flag.  The previous
           counter was incremented in one place and cleared in two, and its
           message carried no inc= token, so nothing reading the log could
           attribute it - both defects are fixed here.  damage_dl_lasthelp is
           the increment of the last Rescue or trust-region firing of ANY
           level; iinc is the increment just committed, so the difference is
           the number of increments since then, and every one of them
           committed without help or lasthelp would have moved. */
        if((damage_dl_narm>0)&&(damage_dl_lasthelp>0)&&
           (iinc>damage_dl_lasthelp)){
          damage_dl_selfrec=iinc-damage_dl_lasthelp;
          if((damage_dl_selfrec>=5)&&
             (damage_dl_recdone!=damage_dl_lasthelp)){
            damage_dl_recdone=damage_dl_lasthelp;
            printf("[DAMAGE TR] SELF-RECOVERY inc=%" ITGFORMAT
                   ": increments %" ITGFORMAT "..%" ITGFORMAT
                   " - five consecutive - committed with PLAIN NEWTON, no "
                   "Rescue and no trust region, since the help at inc=%"
                   ITGFORMAT "; step time is now %.12e.  This, and not "
                   "t_end, is the criterion the method was built to meet.%s",
                   iinc,damage_dl_lasthelp+1,iinc,damage_dl_lasthelp,
                   theta**tper,"\n");
            fflush(stdout);
          }
        }
        if((damage_rec_healthy==damage_rec_window)&&
           (damage_rec_unrec>0)){
          printf("[DAMAGE RESCUE] RECOVERED: %" ITGFORMAT " consecutive "
                 "increments converged with no intervention; the "
                 "un-recovered counter is cleared%s",
                 damage_rec_window,"\n");
          fflush(stdout);
          damage_rec_unrec=0;
        }
      }
      damage_rec_used_in_inc=0;

      if(damage_rescue_bt_on==1){
        damage_rescue_bt_on=0;
        damage_rescue_used=0;
        damage_evt_on=0;
        if(damage_dl_on==1){
          printf("[DAMAGE TR] inc=%" ITGFORMAT " the TRUST-REGION attempt "
                 "CONVERGED at step time %.12e, judged by checkconvergence "
                 "on the UNMODIFIED residual.  The trust region switches OFF and the run "
                 "continues on plain Newton.  Steps kept in this attempt and "
                 "before it: %" ITGFORMAT " Newton, %" ITGFORMAT " Cauchy, %"
                 ITGFORMAT " dogleg; %" ITGFORMAT " iteration(s) accepted "
                 "nothing and kept the full Newton step; %" ITGFORMAT
                 " rejected trial(s); %" ITGFORMAT " residual evaluations; %"
                 ITGFORMAT " armed factorisation(s).  Whether the method "
                 "WORKED is decided by the next five increments, not by "
                 "this line.%s",
                 iinc,theta**tper,damage_dl_nnewt,damage_dl_ncau,
                 damage_dl_ndog,
                 damage_dl_nfail,damage_dl_nrej,damage_dl_neval,
                 damage_dl_nfact,"\n");
          fflush(stdout);
        }
        damage_dl_on=0;
        damage_dl_delta=0.;
        if((damage_reg_on==1)&&(damage_corr_mode==1)&&
           (damage_corr_on==0)){
          damage_corr_on=1;
          damage_corr_lam=damage_reg_lambda;
          damage_corr_lamstable=damage_reg_lambda;
          damage_corr_ninc=0;damage_corr_nint=1;
          damage_corr_nfact=damage_reg_napply;
          damage_corr_clean=0;damage_corr_try=damage_corr_tryevery;
          damage_corr_dtsum=0.;damage_corr_dtn=0;
          damage_corr_trial=0;damage_corr_nsince=0;
          damage_corr_nwallstab=0;damage_corr_nwalltot=0;
          damage_corr_theta0=theta;
          {ITG mi1,mj1;double mtmp,msort[20];
           for(mi1=0;mi1<damage_dth_n;mi1++) msort[mi1]=damage_dth_ring[mi1];
           for(mi1=1;mi1<damage_dth_n;mi1++){
             mtmp=msort[mi1];
             for(mj1=mi1;(mj1>0)&&(msort[mj1-1]>mtmp);mj1--)
               msort[mj1]=msort[mj1-1];
             msort[mj1]=mtmp;
           }
           damage_corr_dtref=(damage_dth_n>0)?msort[damage_dth_n/2]:
                                              damage_dtheta_healthy;
          }
          printf("[DAMAGE CORR] ENTERED at step time %.12e with "
                 "lambda=%.3e after %" ITGFORMAT " factorisation(s).  "
                 "Reference dtime: MEDIAN over the last %" ITGFORMAT
                 " clean increments = %.6e (the last clean increment "
                 "alone was %.6e and sits at the wall, so it is not a "
                 "healthy reference)%s",
                 theta**tper,damage_corr_lam,damage_reg_napply,
                 damage_dth_n,damage_corr_dtref**tper,
                 damage_dtheta_healthy**tper,"\n");
          fflush(stdout);
        }else if(damage_reg_on==1){
          printf("[DAMAGE REG] the regularized attempt CONVERGED after %"
                 ITGFORMAT " factorisation(s) at lambda=%.3e; the shift "
                 "is switched OFF and the run continues stock%s",
                 damage_reg_napply,damage_reg_lambda,"\n");
          fflush(stdout);
        }
        if(damage_corr_on==0) damage_reg_on=0;
        damage_rescue_nok++;
        printf("[DAMAGE RESCUE] ACCEPTED: the rescue increment converged at "
               "step time %.12e; backtracking switched OFF, rescue re-armed "
               "(%" ITGFORMAT " accepted of %" ITGFORMAT " fired)%s",
               theta**tper,damage_rescue_nok,damage_rescue_nfired,"\n");
        fflush(stdout);
      }
	  
      /* ---- [DAMAGE CT] ring, filled passively at the ONE commit point ---
         Six committed endpoints give five intervals.  Endpoint 0 is taken
         the first time this block runs; the next five commits complete the
         set.  Taken BEFORE vini<-vold and xstateini<-xstate, so vold/vini
         still bracket the increment that has just been committed and
         xstate still holds its converged values.  sti is used rather than
         stx for the compression flag: it is the committed stress copy and
         is alive here, whereas stx belongs to the iteration.
         Commit-only, so a cutback cannot corrupt it. */

      if(damage_ct_mode==1){
        ITG cti,ctj,ctn;
        double ctd;
        if(damage_ct_alloc==0){
          NNEW(damage_ct_ring,double,18*mi[0]*ne0);
          NNEW(damage_ct_fl,ITG,6*mi[0]*ne0);
          damage_ct_alloc=1;damage_ct_head=0;damage_ct_nring=0;
          printf("[DAMAGE CT] ring allocated: 6 committed endpoints x %"
                 ITGFORMAT " integration points (%.1f MB).  It fills from "
                 "this commit onwards and changes nothing until arming.%s",
                 mi[0]*ne0,
                 (double)(18*mi[0]*ne0*8+6*mi[0]*ne0*4)/1048576.,"\n");
          fflush(stdout);
        }
        ctn=mi[0]*ne0;
        damage_ct_head=(damage_ct_nring==0)?0:((damage_ct_head+1)%6);
        damage_ct_snap(co,kon,ipkon,lakon,vold,sti,xstate,ne0,mi[0],
                       *nstate_,mt,
                       &damage_ct_ring[3*ctn*damage_ct_head],
                       &damage_ct_fl[ctn*damage_ct_head]);
        if(damage_ct_nring>0){
          ctd=0.;
          for(cti=0;cti<*nk;cti++)
            for(ctj=1;ctj<mt;ctj++)
              if(fabs(vold[mt*cti+ctj]-vini[mt*cti+ctj])>ctd)
                ctd=fabs(vold[mt*cti+ctj]-vini[mt*cti+ctj]);
          damage_ct_duinf[damage_ct_head]=ctd;
          damage_ct_dt[damage_ct_head]=dtime;
          damage_ct_dlam[damage_ct_head]=dtheta;
        }
        if(damage_ct_nring<6) damage_ct_nring++;

        /* ---- [DAMAGE CT] continuation commit lifecycle ----------------
           A committed increment while armed IS one accepted continuation
           step: advance lambda_c and delta_c, count the step, re-establish
           dtime = kappa*ds for the next one, and measure physical progress
           P1/P2/P3.  delta_c itself is re-anchored by the corrector at the
           first iteration of the next step (damage_ct_newstep), so only the
           bookkeeping is done here. */
        if(damage_ct_on==1){
          ITG cp1=0,cp2=0,cp3=0,cq;
          double cpd0;
          damage_ct_step++;damage_ct_ncommit++;
          damage_ct_lamc=damage_ct_lam;
          damage_ct_newstep=1;
          for(cti=0;cti<ne0;cti++){
            if(ipkon[cti]<0){cp2++;continue;}
            if(lakon[8*cti]!='U') continue;
            if(ielprop[cti]<0) continue;
            cpd0=prop[ielprop[cti]+1];
            if(cpd0<=0.) continue;
            cpd0=cpd0/prop[ielprop[cti]];       /* d0 = Tn0/Kn */
            for(ctj=0;ctj<3;ctj++){
              if(ctj>=mi[0]) break;
              cq=mi[0]*cti+ctj;
              if(xstate[*nstate_*cq+3]>=0.5) cp1++;
              if(xstate[*nstate_*cq]>cpd0) cp3++;
            }
          }
          if((cp1>damage_ct_p1r)||(cp2>damage_ct_p2r)||
             (cp3>=damage_ct_p3r+5)) damage_ct_nprog=0;
          else damage_ct_nprog++;
          printf("[DAMAGE CT] COMMIT step=%" ITGFORMAT " inc=%" ITGFORMAT
                 " lambda=%.12e theta=%.12e ds=%.6e dtime=%.6e ; "
                 "P1 failed facets %" ITGFORMAT " (was %" ITGFORMAT
                 "), P2 deleted elements %" ITGFORMAT " (was %" ITGFORMAT
                 "), P3 ip with dmax>d0 %" ITGFORMAT " (was %" ITGFORMAT
                 ") ; quiet windows %" ITGFORMAT "/10 ; evals %" ITGFORMAT
                 " fact %" ITGFORMAT "%s",
                 damage_ct_step,iinc,damage_ct_lam,theta,damage_ct_ds,dtime,
                 cp1,damage_ct_p1r,cp2,damage_ct_p2r,cp3,damage_ct_p3r,
                 damage_ct_nprog,damage_ct_neval,damage_ct_nfact,"\n");
          fflush(stdout);
          damage_ct_p1r=cp1;damage_ct_p2r=cp2;damage_ct_p3r=cp3;
          if((damage_ct_nprog>=10)||
             (damage_ct_step>=damage_ct_maxstep)||
             (damage_ct_neval>=damage_ct_maxeval)||
             (damage_ct_nfact>=damage_ct_maxfact)){
            printf("[DAMAGE CT] segment ENDS: %s.  Requesting the clean "
                   "PARTIAL exit.%s",
                   (damage_ct_nprog>=10)?
                   "no physical progress (P1=P2=P3 unchanged) over 10 "
                   "committed continuation steps":
                   ((damage_ct_step>=damage_ct_maxstep)?"step budget":
                    ((damage_ct_neval>=damage_ct_maxeval)?
                     "residual-evaluation budget":"factorisation budget")),
                   "\n");
            fflush(stdout);
            damage_ct_partial=1;damage_ct_on=0;
          }
        }
      }

      /* previous increment converged: update the initial values */

      /* [LOADPATH] Ask the owner whether that increment was still a
         specimen.

         The site matters as much as the judgement.  The pre-existing
         connectivity test sat inside `if(damage_tent_count>0)`, the bulk
         deletion transaction, so it could only fire in an increment that
         terminally deleted a C3D4 - which means an interface-dominated
         fracture, where facets fail and no bulk erodes, could never reach
         it however carefully it was configured.  Here it is asked once per
         converged increment, unconditionally.

         Cost: one O(nkon+nk) graph walk against a PARDISO solve. */

      if(damage_lp_ready&&(iinc>0)){
        loadpath_census(&damage_lp,ipkon,kon,lakon,ne,xstate,*nstate_,mi[0],
                        stx,co,vold,mt);
        loadpath_note(&damage_lp,iinc,theta**tper);

        /* Policy, applied by the consumer, not by the owner: once the
           specimen has come apart there is nothing left to solve, so the
           run ends.  This reuses damage_fracture_complete, the existing
           and only loop exit, rather than adding a second one.

           CCX_FRACTURE_PAST_SEVERANCE=1 keeps going, and is not a way of
           ignoring the result: every further increment is stamped
           [LOADPATH PHANTOM] and the run's summary says how many there
           were.  test/pathfollow/close.inp needs it - that benchmark
           drives a facet past failure on purpose and then closes it, to
           measure the compressive branch of a crack face, so it is
           severed in tension by design and the interesting half of the
           test comes afterwards. */
        if(loadpath_severed(&damage_lp)&&damage_lp_stop){
          damage_fracture_complete=1;
          /* `continue` goes straight to the loop condition, which
             damage_fracture_complete has just made false, so the run ends
             ON the increment that severed rather than one past it.  The
             last converged state is already in vold/sti and is what the
             post-loop output writes. */
          continue;
        }
      }

      iinc++;
      jprint++;
      damage_active_pass=0;
      damage_soft_reeq=0;
      damage_fast_retry=0;
      damage_fast_used=0;
      damage_fast_recover=0;
      damage_reeq_uam_ref[0]=0.;
      damage_reeq_uam_ref[1]=0.;

      /* A new physical increment starts with no pending DE1.3 terminal
         triggers.  Entries are populated only at the exact topology-change
         event and are preserved across same-load active-set extensions. */
      if(damage_de13_trigger_value!=NULL){
        for(i=0;i<ne0;i++){
          damage_de13_trigger_value[i]=-1.;
          damage_de13_trigger_ip[i]=0;
        }
      }

      /* Establish the end point of the current normal CCX increment.
         If this increment later needs a cutback, theta_goal is kept
         fixed while reduced physical substeps traverse the same interval.
         Do not move the goal while local substepping is already active. */

      if((ilocalsubstep==0)&&(*ndmat_>0)&&(*iexpl<=1)&&
         (*nmethod!=4)&&(*idrct==0)){
        theta_local_start=theta;
        theta_goal=theta+dtheta;
        if(theta_goal>1.) theta_goal=1.;
        dtheta_restore=dtheta;
      }

      /* store number of elements (important for implicit dynamic
	 contact */

      neini=*ne;
	  
      /* ---- CCX_PATHFOLLOW: commit the increment just accepted --------
         This runs exactly once per accepted physical increment, and it
         must run BEFORE vini is overwritten with vold: the increment's
         displacement change is vold-vini, and after the copy that
         difference is identically zero.  A rejected attempt never reaches
         this block, which is what makes the constraint transactional. */

      if((pf_on==1)&&(pf_pending==1)&&(pathfollow_have()==1)){
        pf_pdu=pf_project(pathfollow_fhat(),vold,vini,nactdof,*nk,mt);
        if(pf_codmode==1){
          pathfollow_cod_settarget(pathfollow_cod_target()+pf_dphi);
          pf_dgc=pf_cu;
        }else if(pf_codmode==2){

          /* CONSTRAINT RESIDUAL AT THE ACCEPTED STATE.  Everything else
             the driver prints is measured mid-iteration, before the last
             correction was applied, so it cannot answer "did the
             constraint converge".  This does: vold is the accepted
             state, vini is still the committed one, and the projection
             uses the same c the corrector used. */

          pf_phiacc=pf_project(pathfollow_cod_c(),vold,vini,nactdof,*nk,mt);
          pf_gacc=pf_phiacc-pf_dphicur;
        }
        pathfollow_commit(pf_lam,pf_pdu,&pf_dgc);
        if(pf_engaged==1){
          pf_taucur*=1.4;
          if(pf_taucur>pf_tauv) pf_taucur=pf_tauv;
          pathfollow_settau(pf_taucur);
          if(pf_codmode==2){

            /* 1.4 was measured to be too greedy on the target: dphi grew
               past what the increment could take, attempt 1 failed, dphi
               was halved, attempt 2 converged - one wasted attempt per
               increment, every increment.  Grow gently instead. */

            pf_dphicur*=pf_ccgrow;
            if(pf_dphicur>pf_dphi) pf_dphicur=pf_dphi;
          }
        }
        pf_dlampred=pf_lam-pf_lamprev;
        pf_lamprev=pf_lam;
        if(!(fabs(pf_dlampred)>1.e-30)) pf_dlampred=dtheta;
        if((pf_engaged==0)&&(pf_dgc>0.2*pf_tauv)){
          pf_engaged=1;
          printf("[PATHFOLLOW] engaged at inc=%" ITGFORMAT " lambda=%.8f: "
                 "measured dG=%.6e has reached 0.2*tau=%.6e\n",
                 iinc,pf_lam,pf_dgc,0.2*pf_tauv);
        }
        if(pf_codmode==2){

          /* Census at the ACCEPTED state (vold, xstate), which is what
             "the front advanced" has to be measured on.  pf_cvec is
             refrozen from the committed state at the top of the next
             attempt, so using it as the scratch vector here costs
             nothing. */

          crackcontrol_build(pf_cvec,neq[1],pf_ccmode,co,kon,ipkon,lakon,
                             *ne,ielprop,prop,xstate,*nstate_,mi,vold,
                             nactdof,*nk,mt,&pf_cs);
          printf("[CRACKCTL] inc=%" ITGFORMAT " ACCEPTED lambda=%.8f "
                 "dphi=%.6e achieved=%.6e g=%.3e |R|=%.3e du=%.3e "
                 "zone=%" ITGFORMAT " load=%" ITGFORMAT " init=%"
                 ITGFORMAT " fail=%" ITGFORMAT " dead=%" ITGFORMAT
                 " deffmax=%.6e shear=%.4f w=%.6e fhat=%.6f/%.4f\n",
                 iinc,pf_lam,pf_dphicur,pf_phiacc,pf_gacc,ram[0],ram[1],
                 pf_cs.nzone,pf_cs.nload,pf_cs.ninit,pf_cs.nfail,
                 pf_cs.ndead,pf_cs.deffmax,pf_cs.shearfrac,pf_cs.weight,
                 pf_fhcos,pf_fhrat);
        }else if(pf_codmode==1){
          printf("[PATHFOLLOW] inc=%" ITGFORMAT " ACCEPTED lambda=%.8f "
                 "phi=%.6e target=%.6e\n",iinc,pf_lam,pf_cu,
                 pathfollow_cod_target());
        }else{
          printf("[PATHFOLLOW] inc=%" ITGFORMAT " ACCEPTED lambda=%.8f "
                 "dG=%.6e P=%.6e ff=%.6e engaged=%" ITGFORMAT " refusals=%"
                 ITGFORMAT "\n",iinc,pf_lam,pf_dgc,pathfollow_Pn(),
                 pathfollow_ff(),pf_engaged,pathfollow_refusals());
        }
        fflush(stdout);
      }
      pf_pending=pf_on;

      /* vold is copied into vini */
	  
      isiz=mt**nk;cpypardou(vini,vold,&isiz,&num_cpus);
	  
      isiz=*nboun;cpypardou(xbounini,xbounact,&isiz,&num_cpus);


      if((*ithermal==1)||(*ithermal>=3)){
	isiz=*nk;cpypardou(t1ini,t1act,&isiz,&num_cpus);
      }
      isiz=neq[1];cpypardou(fini,f,&isiz,&num_cpus);
      if(*nmethod==4){
	if(*iexpl<=1){
	  isiz=mt**nk;
	  cpypardou(veini,veold,&isiz,&num_cpus);
	  cpypardou(accini,accold,&isiz,&num_cpus);
	}
	isiz=mt**nk;cpypardou(fnextini,fnext,&isiz,&num_cpus);

	isiz=neq[1];
	cpypardou(fextini,fext,&isiz,&num_cpus);
	cpypardou(cvini,cv,&isiz,&num_cpus);
	      
	if(*ithermal<2){
	  allwkini=allwk;
	  // MPADD start
	  if(idamping==1)dampwkini = dampwk;
	  for(k=0;k<4;k++){
	    energyini[k]=energy[k];
	  }
	  // MPADD end
	}
      }
      if(*ithermal!=2){
	isiz=6*mi[0]*ne0;	    
	cpypardou(stiini,sti,&isiz,&num_cpus);
	cpypardou(emeini,eme,&isiz,&num_cpus);
      }

      /* the contact friction energy is stored at the slave nodes for
         mortar not equal to 1 */
      
      if(*nener==1){
	  isiz=2*mi[0]*ne0;
	cpypardou(enerini,ener,&isiz,&num_cpus);
      }
	      

      if(*mortar!=1){
	if(*nstate_!=0){
	  isiz=*nstate_*mi[0]*(ne0+*nslavs);
	  cpypardou(xstateini,xstate,&isiz,&num_cpus);
	}
      }

      /* store topology and damage at the start of the physical increment */

      if((*ndmat_>0)&&(*iexpl<=1)){
	isiz=ne0;
	cpyparitg(ipkondamageini,ipkon,&isiz,&num_cpus);
	isiz=mi[0]*ne0;
	cpypardou(damdamageini,dam,&isiz,&num_cpus);
	if(damage_damviscini!=NULL){
	  cpypardou(damage_damviscini,damage_damvisc,&isiz,&num_cpus);
	}
	if(*nmethod!=4){
	  isiz=mt**nk;
	  cpypardou(veolddamageini,veold,&isiz,&num_cpus);
	}
      }
	
      /* From this point on every results() call in the physical increment
         reconstructs DE1.2 trial damage from the immutable committed
         baseline.  The public results() ABI remains unchanged. */
      results_set_de12_context(damage_de12_enabled,damage_tangent_mode,
                               ndmat_,ndmcon,dmcon,damdamageini,
                               damage_damjac,damage_damvisc,
                               damage_damviscini,damage_visc_eta);

      if(*mortar>1){
	for (i=0;i<*ntie;i++){
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    islavactini[j]=islavact[j];
	    bpini[j]=bp[j];
	    for(k=0;k<mt;k++){
	      cstressini[mt*j+k]=cstress[mt*j+k];
	    }		      
	  }    
	}
      }
    }
      
    /* check for max. # of increments */
      
    if(iinc>jmax[0]){
      printf(" *ERROR in nonlingeo: max. # of increments reached\n\n");
      FORTRAN(stop,());
    }

    if(*iexpl<=1){
      printf(" increment %" ITGFORMAT " attempt %" ITGFORMAT " \n",iinc,icutb+1);
      printf(" increment size= %e\n",dtheta**tper);
      printf(" sum of previous increments=%e\n",theta**tper);
      printf(" actual step time=%e\n",(theta+dtheta)**tper);
      printf(" actual total time=%e\n\n",*ttime+(theta+dtheta)**tper);
      
      printf(" iteration 1\n\n");
    }
      
    qamold[0]=qam[0];
    qamold[1]=qam[1];

    icntrl=0;

    /* restoring the distributed loading before adding the
       friction heating */

    if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
      *nload=nloadref;
      isiz=2**nload;cpyparitg(nelemload,nelemloadref,&isiz,&num_cpus);
      if(*nam>0){
	isiz=2**nload;cpyparitg(iamload,iamloadref,&isiz,&num_cpus);
      }
      memcpy(&sideload[0],&sideloadref[0],sizeof(char)*20**nload);
    }
      
    /* store the load level in case damage requires re-equilibration */

    if(idamagereeq==0){
      thetadamage=theta;
      dthetadamage=dtheta;
      dthetarefdamage=dthetaref;
    }

    /* determining the actual loads at the end of the new increment*/

    /* Once the constraint drives lambda, theta is only a monotone counter
       that keeps dtime positive and paces the output.  The stock
       controller still sizes it, but it is capped here: theta reaching 1
       ends the step, and an uncapped theta would end the run in the middle
       of the snap-back while lambda still has most of the branch to go. */

    /* CAP, not set.  The original code SET dtheta to a constant once the
       constraint took over, on the argument that theta is only a counter
       when lambda is decoupled from it, so shrinking it cannot make a
       retry easier.  That argument holds only for a rate-INDEPENDENT
       model.

       The s3rad target is not one: CCX_DAMAGE_VISCOSITY=1e-4 and the UC6
       law both relax with alpha = dt/(mu+dt), so dtime is a constitutive
       parameter and a theta cutback really does change the problem.
       Measured: the stock run needs dtime around 1e-6 at increment 250
       (alpha=0.0095) and recovers there by cutting it; pinning dtheta at
       1e-4 holds alpha at 0.5, i.e. removes almost all of the viscous
       regularisation, and the first engaged attempt failed through five
       cutbacks that could not touch dtime.

       Capping it and leaving the stock controller to size it below the
       cap was tried and MEASURED, and it does not work either, for a
       different reason: once engaged, an attempt that fails is fixed by
       halving dphi, not by halving dtheta, but the stock controller
       halves dtheta anyway and then refuses to grow it back because the
       accepted attempt took 29 iterations.  dtheta therefore ratchets
       down one notch per increment and the run dies on "increment size
       smaller than minimum" - measured, 14 accepted engaged increments
       and then that exact stop.

       So dtheta is SET again, and it is a CHOICE of time step: with
       lambda decoupled, theta is the clock the rate-dependent laws run
       on and nothing else, and holding it fixed makes the viscous
       relaxation per increment constant.  The continuation step size is
       dphi and that alone is what a cutback shrinks. */

    if((pf_on==1)&&(pf_engaged==1)){
      dtheta=pf_dtheta_eng;
      if(theta+dtheta>1.) dtheta=1.-theta;
    }

    reltime=theta+dtheta;

    /* Content hash of the COMMITTED state at the top of this attempt.
       vini/xstateini/damdamageini/ipkondamageini are exactly what a
       rollback restores, so two attempts that print the same hash really
       did start from the same state - which is the half of the
       "repeating event" claim that cannot be read off the existing log. */

    if((td_armed!=0)&&(td_trace!=0)){
      td_state=topodiag_hash_seed();
      td_state=topodiag_hash_d(vini,mt**nk,td_state);
      if(xstateini!=NULL)
        td_state=topodiag_hash_d(xstateini,*nstate_*mi[0]**ne,td_state);
      if(damdamageini!=NULL)
        td_state=topodiag_hash_d(damdamageini,mi[0]*ne0,td_state);
      if(ipkondamageini!=NULL)
        td_state=topodiag_hash_i(ipkondamageini,ne0,td_state);
      td_state=topodiag_hash_i(ipkon,*ne,td_state);
      printf("[BATCHTRACE] attempt inc=%" ITGFORMAT " icutb=%" ITGFORMAT
             " theta=%.12e dtheta=%.12e committed_state=%016llx\n",
             iinc,icutb,theta,dtheta,td_state);
      fflush(stdout);
    }

    /* [WALLDIAG] theta is the increment's starting load factor here, before
       dtheta is added, so the gate opens on the first attempt whose base
       state is at or beyond the recorded wall. */

    if((damage_wall_theta>=0.)&&(theta>=damage_wall_theta)){
      if(damage_wall_armed==0){
        damage_wall_armed=1;
        printf("[WALLDIAG] GATE OPEN at inc=%" ITGFORMAT " icutb=%" ITGFORMAT
               " theta=%.12e dtheta=%.12e%s",iinc,icutb,theta,dtheta,"\n");
        fflush(stdout);
      }
      damage_dl_lincheck=iinc;
      damage_ls_probe=1;
      if(td_from<=0) td_from=iinc;

      /* NOT the null-vector probe.  CCX_DAMAGE_NULLVEC ends with
         stopwithout201() by design - it is a terminal measurement, taken
         once, on a run that is not meant to continue.  Arming it from a
         load-factor gate killed the run at the gate instead of at the wall
         - measured, R0-early stopped at increment 92 with the ladder
         unmeasured.  It stays opt-in and its terminal nature is stated. */
      if((damage_wall_null!=0)&&(damage_null_inc<=0)){
        damage_null_inc=iinc;
        printf("[WALLDIAG] arming the deflated null-vector probe on this "
               "increment.  IT IS TERMINAL: the run stops after it "
               "prints.%s","\n");
        fflush(stdout);
      }
    }

    FORTRAN(uc6setinc,(&iinc));

    /* trial load factor for dissipation control; theta itself is
       left alone until the increment is accepted */

    if(damage_diss_ctrl>=1){
      damage_diss_lamcur=(damage_arc==1)?(damage_arc_lam+dtheta):(theta+dtheta);
      damage_diss_have=0;
    }
    time=reltime**tper;
    dtime=dtheta**tper;
      
    FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
		      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
		      xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,nam,ampli,
		      &time,&reltime,ttime,&dtime,ithermal,nmethod,xbounold,
		      xboun,xbounact,iamboun,nboun,nodeboun,ndirboun,nodeforc,
		      ndirforc,istep,&iinc,co,vold,itg,&ntg,amname,ikboun,
		      ilboun,nelemload,sideload,mi,ntrans,trab,inotr,veold,
		      integerglob,doubleglob,tieset,istartset,iendset,ialset,
		      ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc,ipobody,
		      iponoeln,inoeln,ipkon,kon,ielprop,prop,ielmat,shcon,nshcon,
		      rhcon,nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));

    /* tempload has just built xbounact from reltime=theta+dtheta.  With
       lambda decoupled that is no longer the load factor, so rebuild the
       prescribed values from lambda itself.  xbounold is the value at the
       START of the step, so this is an absolute step fraction - exactly the
       ramp CCX_DAMAGE_PATH mode 1 verified against tempload at max|diff|=0. */

    if(damage_arc==1){
      for(k=0;k<*nboun;k++){
        xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*damage_diss_lamcur;
      }
    }

    /* ---- CCX_PATHFOLLOW: predictor for this attempt -----------------
       pathfollow_incstart discards whatever a rejected attempt left, so
       a cutback re-enters here from the committed state.  Before the
       constraint engages lambda simply follows the step time, which is
       the ordinary control the stock solver would have run. */

    if(pf_on==1){
      if(neq[1]!=pf_neqarm){
        if(pathfollow_resize(neq[1])==0){
          printf("[PATHFOLLOW] disarmed: could not follow the equation "
                 "count from %" ITGFORMAT " to %" ITGFORMAT "\n",
                 pf_neqarm,neq[1]);
          pf_on=0;
        }else{
          RENEW(pf_uf,double,neq[1]);
          printf("[PATHFOLLOW] equation count %" ITGFORMAT " -> %" ITGFORMAT
                 "; constraint origin restarted at the current state\n",
                 pf_neqarm,neq[1]);
          pf_neqarm=neq[1];

          /* The control functional lives in EQUATION space.  After a
             remastruct its entries refer to dofs that no longer exist,
             so a vector kept from before the change is not merely stale,
             it is read out of bounds.  The mixed-mode functional is
             rebuilt from the model a few lines below, every attempt; the
             absolute Mode-I functional cannot be rebuilt without redoing
             the geometry pass, so it is disarmed loudly rather than used
             wrong. */

          if(pf_codmode==2){
            RENEW(pf_cvec,double,neq[1]);
          }else if(pf_codmode==1){
            printf("[PATHFOLLOW] *ERROR: CCX_PATHFOLLOW_COD builds its "
                   "control functional once, from a fixed equation "
                   "numbering.  The equation count has changed, so the "
                   "vector no longer refers to the same dofs.  "
                   "Disarming rather than continuing with a stale "
                   "constraint; use CCX_CRACK_CONTROL, which refreezes "
                   "the functional every attempt.\n");
            pf_on=0;pf_engaged=0;pf_codmode=0;
          }
        }
      }
    }
    if(pf_on==1){

      /* The dissipation increment IS the step size here.  A stock cutback
         shrinks theta, but lambda no longer follows theta, so the retry
         would be the identical problem.  Shrink tau instead - that is the
         path-following cutback, and without it an engaged increment can
         only fail down to the minimum step time. */

      if((pf_engaged==1)&&(icutb>pf_icutbprev)){
        pf_taucur*=0.5;
        if(pf_taucur<1.e-6*pf_tauv) pf_taucur=1.e-6*pf_tauv;
        pf_ncut++;
        pathfollow_settau(pf_taucur);

        /* For crack control the CONTROL INCREMENT is the step size, for
           exactly the same reason: lambda no longer follows theta, so a
           theta cutback retries the identical problem. */

        if(pf_codmode==2){
          pf_dphicur*=0.5;
          if(pf_dphicur<1.e-8*pf_dphi){
            pf_dphicur=1.e-8*pf_dphi;

            /* TERMINATION.  With dtheta held fixed the stock stop
               "increment size smaller than minimum" can never fire, so
               without this the run has no termination criterion at all
               and a wall becomes an infinite loop instead of a reported
               result.  Measured on the target: the deletion batch at
               increment 348 leaves an orphan node, its same-load
               re-equilibration does not converge, the batch is rolled
               back, and the identical event repeats forever.

               dphi at its floor and still failing means the obstruction
               is not the size of the continuation step. */

            pf_ncutfloor++;
            if(pf_ncutfloor>=3){
              printf("[CRACKCTL] *WALL: the control increment is at its "
                     "floor (%.6e, 1e-8 of the requested %.6e) and the "
                     "attempt still fails.  The obstruction is not the "
                     "size of the continuation step.  Stopping.\n",
                     pf_dphicur,pf_dphi);
              fflush(stdout);
              FORTRAN(stop,());
            }
          }else{
            pf_ncutfloor=0;
          }
          printf("[CRACKCTL] cutback: dphi -> %.6e\n",pf_dphicur);
        }

        /* the predictor is proportional to tau, so it rescales itself */

        printf("[PATHFOLLOW] cutback %" ITGFORMAT ": tau -> %.6e\n",
               pf_ncut,pf_taucur);
        fflush(stdout);
      }
      pf_icutbprev=icutb;
      pathfollow_incstart();

      /* ---- mixed-mode crack control: refreeze the functional --------
         Built from the COMMITTED state (vini, xstateini), which is what
         a rollback restores, so a retry rebuilds the identical vector.
         The constraint is incremental, so the target is simply the
         control increment currently in force. */

      if(pf_codmode==2){
        pf_ccnw=crackcontrol_build(pf_cvec,neq[1],pf_ccmode,co,kon,ipkon,
                                   lakon,*ne,ielprop,prop,xstateini,
                                   *nstate_,mi,vini,nactdof,*nk,mt,&pf_cs);

        /* ZONE and DISS are supported on the process zone.  Before the
           first initiation, and again if every initiated point has
           failed, that set is empty and the functional carries no
           weight.  Fall back to the mean over all live facets for this
           attempt rather than handing the bordered row a zero vector. */

        if((pf_ccnw==0)&&(pf_ccmode!=0)){
          pf_ccnw=crackcontrol_build(pf_cvec,neq[1],0,co,kon,ipkon,
                                     lakon,*ne,ielprop,prop,xstateini,
                                     *nstate_,mi,vini,nactdof,*nk,mt,
                                     &pf_cs);
        }
        pathfollow_cod_arm(pf_cvec,neq[1]);
        pathfollow_cod_settarget(pf_dphicur);

        /* With the extrapolation suppressed, b/dlamjump at the first
           iteration IS -dR/dlambda, so there is no reason to keep a
           vector captured many increments and one remastruct ago.  The
           freeze existed to stop the dissipation constraint rescaling
           its own tau; the crack-control constraint does not depend on
           the scale of f_hat, only on its being the true derivative. */

        if(pf_engaged==1) pathfollow_unfreeze();

        if((pf_engaged==0)&&(pf_ccnw>0)){
          if(((pf_ccengage>0)&&(iinc>=pf_ccengage))||
             ((pf_ccengage<=0)&&(pf_cs.nzone>0))){
            pf_engaged=1;

            /* The cap above was applied before this decision was taken,
               so the attempt that engages would otherwise run with the
               uncapped stock dtheta - which is exactly the attempt where
               the constitutive relaxation must not jump.  Apply it here
               too.  Only the time-like quantities are rebuilt: xbounact
               is overwritten from pf_lam a few lines below anyway. */

            dtheta=pf_dtheta_eng;
            if(theta+dtheta>1.) dtheta=1.-theta;
            reltime=theta+dtheta;
            time=reltime**tper;
            dtime=dtheta**tper;
            printf("[CRACKCTL] engaged at inc=%" ITGFORMAT " lambda=%.8f: "
                   "process zone %" ITGFORMAT " ip(s) of which %"
                   ITGFORMAT " loading, initiated %" ITGFORMAT
                   ", failed %" ITGFORMAT ", shear fraction %.4f, "
                   "weight %.6e\n",iinc,pathfollow_lamn(),pf_cs.nzone,
                   pf_cs.nload,pf_cs.ninit,pf_cs.nfail,pf_cs.shearfrac,
                   pf_cs.weight);
            fflush(stdout);
          }
        }
      }

      if(pf_engaged==0){
        pf_lam=theta+dtheta;
      }else{

        /* Tangent predictor, equation (7) in pathfollow.c.  Its SIGN is
           not chosen here: it comes out of lambda_n*ff-P_n, which changes
           sign at the limit point.  That is what turns the branch. */

        if(pf_codmode==0){
          if(pathfollow_predictor(&pf_dlampred)==0){
            printf("[PATHFOLLOW] predictor degenerate at inc=%" ITGFORMAT
                   "; keeping dlambda=%+.4e\n",iinc,pf_dlampred);
          }
        }
        /* For crack-opening control the previous accepted dlambda is the
           predictor: it carries the sign of the branch, and the constraint
           corrects the length in one step because phi is linear in u.

           It must not be zero on the first engaged increment: f_hat is
           captured by dividing the first residual by this jump, so a zero
           predictor leaves f_hat undefined, the corrector never runs and
           lambda sits still while the step controller burns its cutbacks.
           Measured as "too many cutbacks" with no [PATHFOLLOW] it= line at
           all. */
        /* PROBE JUMP, not predictor jump, once the crack constraint is
           driving.

           f_hat is obtained as b/dlamjump at the first iteration.  With
           the displacement extrapolation suppressed that quotient is a
           SECANT of -dR/dlambda over the jump, and it is the tangent only
           in the limit of a small jump.  Using the previous accepted
           dlambda as that jump makes it a secant over the whole physical
           step: at increment 60 of the target the step is 2e-03, over
           which the boundary-adjacent material yields and damages, and
           an honest finite difference of the residual put
           |f_hat|/|dR/dlambda| at 1.43.

           A scale error there is not harmless.  The corrector adds
           dlambda*K^-1 f_hat to u and dlambda to lambda; if f_hat = s*q
           the DISPLACEMENT half is right and the LOAD half is a factor
           1/s short, so the prescribed dofs end up where the free dofs
           did not assume.  Measured on the target: the constraint row
           converged to 1e-19 every iteration while the equilibrium
           residual stalled at 0.17 - and it stalled on nodes next to the
           loaded face, which is exactly where -dR/dlambda lives.  With
           dlambda clamped to 1e-9 the same increment converged to 8e-04
           and was accepted, so the machinery is sound and the reference
           vector is what is wrong.

           So take the jump as a small PROBE of size pf_eps.  It costs
           nothing - no extra residual evaluation, no state to restore -
           because the first iteration of the increment has to evaluate a
           residual anyway.  The constraint then supplies the real
           predictor: its own dlambda closes the whole step in one
           corrector, since phi is linear in u. */

        if(pf_codmode==2) pf_dlampred=pf_eps;
        if(!(fabs(pf_dlampred)>1.e-12)) pf_dlampred=dtheta;
        pf_lam=pathfollow_lamn()+pf_dlampred;
      }
      pf_dlamjump=pf_lam-pathfollow_lamn();
      for(k=0;k<*nboun;k++){
        xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*pf_lam;
      }

      /* Bisection probe: is the model state at the START of an attempt
         really the committed one?  vold-vini must be identically zero
         here; anything else means the rollback did not restore what the
         constraint assumes it restored. */

      if(getenv("CCX_PATHFOLLOW_ACCUMCHECK")!=NULL){
        double pfd=0.,pfm=0.;
        ITG pfi,pfj,pfk;
        for(pfi=0;pfi<*nk;pfi++){
          for(pfj=1;pfj<mt;pfj++){
            pfk=nactdof[mt*pfi+pfj];
            if(pfk>0){
              pfd+=(vold[mt*pfi+pfj]-vini[mt*pfi+pfj])*
                   (vold[mt*pfi+pfj]-vini[mt*pfi+pfj]);
              if(fabs(vold[mt*pfi+pfj]-vini[mt*pfi+pfj])>pfm)
                pfm=fabs(vold[mt*pfi+pfj]-vini[mt*pfi+pfj]);
            }
          }
        }
        printf("[PF-START] inc=%" ITGFORMAT " attempt icutb=%" ITGFORMAT
               " |vold-vini|=%.6e max=%.6e lambda=%.8f\n",
               iinc,icutb,sqrt(pfd),pfm,pf_lam);
        fflush(stdout);
      }
    }


    /* ================= [DAMAGE CT] ARMING =========================
       Runs once, at the top of the attempt that follows the level-4
       firing.  The standard cutback rollback has already restored the
       increment-start state, so every quantity below is read from the last
       COMMITTED state.  ANY refusal leaves that state untouched, sets
       damage_ct_refused and hands the wall to the ORIGINAL stock stop. */

    if(damage_ct_arm==1){
      ITG ai,aj,ak,an,anp,anm,abest=-1,aip=-1,acnt=0,aok,ansel=0;
      ITG afl0,afl1,asupp[18],anod[18],adir[18];
      double aw[18],adl[3],armat[9],ash[3],ad0[3],ad1[3];
      double am[3],adeff,abeta,atn0,ats0,agc,atau,adf;
      double aq[5],akp[5],atmp,abestv=-1.,aclam,ag[3];
      ITG ah,ahp,as0,as1,ne_uc6=0;

      damage_ct_arm=0;
      damage_ct_refused=1;               /* pessimistic until every guard passes */

      /* ---- admissibility domain (SPEC FREEZE v1) ---- */
      aok=1;
      if(*nam!=0) aok=0;  if(*nforc!=0) aok=0;  if(*nload!=0) aok=0;
      if(*nbody!=0) aok=0; if(*nmpc!=0) aok=0;  if(ncont!=0) aok=0;
      if(*mortar>1) aok=0; if(*ithermal>=2) aok=0; if(*nmethod!=1) aok=0;
      if(*idrct!=0) aok=0; if(*iprestr!=0) aok=0; if(*nstate_<4) aok=0;
      if(*isolver!=7) aok=0;
      {
        char *ae=getenv("CCX_PARDISO_REUSE_SYMBOLIC");
        if((ae==NULL)||((strcmp(ae,"1")!=0)&&(strcmp(ae,"ON")!=0)&&
                        (strcmp(ae,"on")!=0)&&(strcmp(ae,"YES")!=0)&&
                        (strcmp(ae,"yes")!=0))) aok=0;
      }
      for(ai=0;ai<ne0;ai++)
        if((ipkon[ai]>=0)&&(lakon[8*ai]=='U')) ne_uc6++;
      if(ne_uc6==0) aok=0;
      if(damage_ct_nring<6) aok=0;
      if(aok==0){
        printf("[DAMAGE CT] REFUSING TO ARM: admissibility domain or ring not "
               "satisfied (nam=%" ITGFORMAT " nforc=%" ITGFORMAT " nload=%"
               ITGFORMAT " nbody=%" ITGFORMAT " nmpc=%" ITGFORMAT " ncont=%"
               ITGFORMAT " mortar=%" ITGFORMAT " ithermal=%" ITGFORMAT
               " nmethod=%" ITGFORMAT " idrct=%" ITGFORMAT " iprestr=%"
               ITGFORMAT " nstate_=%" ITGFORMAT " isolver=%" ITGFORMAT
               " uc6_elements=%" ITGFORMAT " ring=%" ITGFORMAT "/6, "
               "PARDISO symbolic reuse required).  The last committed state "
               "is untouched; this wall goes to the ORIGINAL stock stop.%s",
               *nam,*nforc,*nload,*nbody,*nmpc,ncont,*mortar,*ithermal,
               *nmethod,*idrct,*iprestr,*nstate_,*isolver,ne_uc6,
               damage_ct_nring,"\n");
        fflush(stdout);
      }else{

      /* ---- candidate scan over the five committed intervals ---- */
      ah=damage_ct_head;
      an=mi[0]*ne0;
      printf("[DAMAGE CT] arming at inc=%" ITGFORMAT " theta=%.12e: scanning "
             "%" ITGFORMAT " UC6 elements over 6 committed endpoints.%s",
             iinc,theta,ne_uc6,"\n");
      for(ai=0;ai<ne0;ai++){
        if(ipkon[ai]<0) continue;
        if(lakon[8*ai]!='U') continue;
        if(ielprop[ai]<0) continue;
        atn0=prop[ielprop[ai]+1];
        ats0=prop[ielprop[ai]+2];
        agc =prop[ielprop[ai]+3];
        if((atn0<=0.)||(ats0<=0.)||(agc<=0.)) continue;
        abeta=(ats0/atn0)*(ats0/atn0);
        adf=2.*agc/atn0;
        atau=1.e-9*adf;
        for(aj=0;aj<3;aj++){
          ak=an*ah+ (mi[0]*ai+aj);
          afl0=damage_ct_fl[an*ah+mi[0]*ai+aj];
          if(afl0<0) continue;                   /* not a live UC6 ip      */
          if(afl0&1) continue;                   /* already failed         */
          /* frozen mixed-mode direction m at the newest committed endpoint */
          ad1[0]=damage_ct_ring[3*(an*ah+mi[0]*ai+aj)];
          ad1[1]=damage_ct_ring[3*(an*ah+mi[0]*ai+aj)+1];
          ad1[2]=damage_ct_ring[3*(an*ah+mi[0]*ai+aj)+2];
          adeff=(ad1[0]>0.?ad1[0]*ad1[0]:0.)
                +abeta*(ad1[1]*ad1[1]+ad1[2]*ad1[2]);
          adeff=(adeff>0.)?sqrt(adeff):0.;
          if(!(adeff>0.)) continue;
          if(fabs(ad1[0])<0.05*adeff) continue;  /* too close to the kink  */
          am[0]=(ad1[0]>0.?ad1[0]:0.)/adeff;
          am[1]=abeta*ad1[1]/adeff;
          am[2]=abeta*ad1[2]/adeff;
          /* five historical intervals, ALL must be admissible */
          aok=1;
          for(ak=0;ak<5;ak++){
            as1=(ah-ak+6)%6;  as0=(ah-ak-1+6)%6;
            afl1=damage_ct_fl[an*as1+mi[0]*ai+aj];
            afl0=damage_ct_fl[an*as0+mi[0]*ai+aj];
            if((afl1<0)||(afl0<0)){aok=0;break;}
            if((afl1&1)||(afl0&1)){aok=0;break;}        /* alive both ends */
            if((afl1&2)!=(afl0&2)){aok=0;break;}        /* ucomp unchanged */
            for(anp=0;anp<3;anp++){
              ad1[anp]=damage_ct_ring[3*(an*as1+mi[0]*ai+aj)+anp];
              ad0[anp]=damage_ct_ring[3*(an*as0+mi[0]*ai+aj)+anp];
            }
            if((ad1[0]>0.)!=(ad0[0]>0.)){aok=0;break;}  /* sign(dn) kept   */
            atmp=am[0]*(ad1[0]-ad0[0])+am[1]*(ad1[1]-ad0[1])
                +am[2]*(ad1[2]-ad0[2]);
            if(!(atmp>atau)){aok=0;break;}
            aq[ak]=atmp;
            if(!(damage_ct_dt[as1]>0.)){aok=0;break;}
            akp[ak]=damage_ct_dt[as1]/atmp;
            if(!(akp[ak]>0.)){aok=0;break;}
          }
          if(aok==0) continue;
          acnt++;
          if(aq[0]>abestv){
            abestv=aq[0];abest=ai;aip=aj;
            for(anp=0;anp<3;anp++) damage_ct_m[anp]=am[anp];
            for(anp=0;anp<5;anp++){
              damage_ct_dt[anp]=damage_ct_dt[anp];      /* untouched      */
            }
            /* keep this candidate's five ratios and advances */
            for(anp=0;anp<5;anp++){damage_ct_g[0]=0.;}
            for(anp=0;anp<5;anp++){aq[anp]=aq[anp];akp[anp]=akp[anp];}
            /* median of five: insertion sort on local copies */
            {
              double sq[5],sk[5];
              ITG bi,bj2;
              for(bi=0;bi<5;bi++){sq[bi]=aq[bi];sk[bi]=akp[bi];}
              for(bi=1;bi<5;bi++){
                atmp=sq[bi];for(bj2=bi;(bj2>0)&&(sq[bj2-1]>atmp);bj2--)
                  sq[bj2]=sq[bj2-1]; sq[bj2]=atmp;
                atmp=sk[bi];for(bj2=bi;(bj2>0)&&(sk[bj2-1]>atmp);bj2--)
                  sk[bj2]=sk[bj2-1]; sk[bj2]=atmp;
              }
              damage_ct_ds0=sq[2];
              damage_ct_kappa=sk[2];
              damage_ct_tau=(sk[0]>0.)?sk[4]/sk[0]:1.e30;
            }
          }
        }
      }
      printf("[DAMAGE CT] candidates with five admissible intervals: %"
             ITGFORMAT "%s",acnt,"\n");
      if((acnt==0)||(abest<0)){
        printf("[DAMAGE CT] REFUSING TO ARM: no UC6 integration point has "
               "five consecutive admissible committed intervals at this "
               "state.  The last committed state is untouched; this wall "
               "goes to the ORIGINAL stock stop.%s","\n");
        fflush(stdout);
      }else if(!(damage_ct_tau<=damage_ct_kaptol)){
        printf("[DAMAGE CT] REFUSING TO ARM: kappa is not stable over the "
               "five intervals of the best candidate (max/min=%.4f > 1.5, "
               "median kappa=%.6e).  Refusing rather than freezing an "
               "unmeasured control-to-time coupling.%s",
               damage_ct_tau,damage_ct_kappa,"\n");
        fflush(stdout);
      }else{
        /* ---- freeze the control point, m, the stencil and c_lambda ---- */
        damage_ct_elem=abest;damage_ct_ip=aip;
        damage_ct_kin(co,kon,ipkon[abest],vold,mt,aip,adl,armat,ash);
        for(ak=0;ak<3;ak++){
          ag[ak]=damage_ct_m[0]*armat[ak]+damage_ct_m[1]*armat[3+ak]
                +damage_ct_m[2]*armat[6+ak];
          damage_ct_g[ak]=ag[ak];
          /* delta_c is taken from the SAME committed snapshot m came from
             - the newest ring endpoint - not recomputed from vold.  Mixing
             the two sources is what made the first constraint value differ
             from -ds by 83x on the short deck: the ring endpoint and vold
             are not the same instant when the wall is reached in a later
             increment than the last commit. */
          damage_ct_dc[ak]=damage_ct_ring[3*(mi[0]*ne0*ah+mi[0]*abest+aip)+ak];
        }
        printf("[DAMAGE CT] delta_c source check: ring=(%.9e,%.9e,%.9e) "
               "vold=(%.9e,%.9e,%.9e) m.(vold-ring)=%.6e%s",
               damage_ct_dc[0],damage_ct_dc[1],damage_ct_dc[2],
               adl[0],adl[1],adl[2],
               damage_ct_m[0]*(adl[0]-damage_ct_dc[0])
              +damage_ct_m[1]*(adl[1]-damage_ct_dc[1])
              +damage_ct_m[2]*(adl[2]-damage_ct_dc[2]),"\n");
        damage_ct_nsupp=0;aclam=0.;aok=1;
        for(ai=0;ai<3;ai++){
          anm=kon[ipkon[abest]+ai]-1;
          anp=kon[ipkon[abest]+ai+3]-1;
          for(ak=0;ak<3;ak++){
            anod[damage_ct_nsupp]=anp;adir[damage_ct_nsupp]=ak+1;
            aw[damage_ct_nsupp]= ash[ai]*ag[ak];damage_ct_nsupp++;
            anod[damage_ct_nsupp]=anm;adir[damage_ct_nsupp]=ak+1;
            aw[damage_ct_nsupp]=-ash[ai]*ag[ak];damage_ct_nsupp++;
          }
        }
        if(damage_ct_w==NULL) NNEW(damage_ct_w,double,18);
        for(ai=0;ai<damage_ct_nsupp;ai++){
          ak=nactdof[mt*anod[ai]+adir[ai]];
          damage_ct_w[ai]=aw[ai];
          damage_ct_node[ai]=anod[ai];damage_ct_dir[ai]=adir[ai];
          if(ak>0){
            asupp[ai]=ak-1;
          }else{
            asupp[ai]=-1;
            if(ak<0){
              if(ak!=2*(ak/2)) aok=0;          /* odd negative = MPC       */
            }
            /* SPC: find its boundary entry and add w * d(u_c)/d(lambda) */
            for(aj=0;aj<*nboun;aj++){
              if((nodeboun[aj]-1==anod[ai])&&(ndirboun[aj]==adir[ai])){
                aclam+=aw[ai]*(xboun[aj]-xbounold[aj]);
                break;
              }
            }
          }
          damage_ct_supp[ai]=asupp[ai];
        }
        damage_ct_clam=aclam;
        if(aok==0){
          printf("[DAMAGE CT] REFUSING TO ARM: an MPC-dependent degree of "
                 "freedom is in the control stencil; the MPC expansion is "
                 "not implemented.  ORIGINAL stock stop.%s","\n");
          fflush(stdout);
        }else{
          damage_ct_lamc=theta;
          damage_ct_lam=theta;
          damage_ct_ds=damage_ct_ds0;
          damage_ct_dsmin=1.e-3*damage_ct_ds0;
          damage_ct_dsmax=20.*damage_ct_ds0;
          damage_ct_tolc=1.e-3*damage_ct_ds0;
          if(damage_ct_tolc<1.e-9) damage_ct_tolc=1.e-9;
          /* healthy references for the absolute predictor bounds */
          {
            double sl[5],su[5];ITG bi,bj2;
            for(bi=0;bi<5;bi++){
              as1=(ah-bi+6)%6;
              sl[bi]=damage_ct_dlam[as1];su[bi]=damage_ct_duinf[as1];
            }
            for(bi=1;bi<5;bi++){
              atmp=sl[bi];for(bj2=bi;(bj2>0)&&(sl[bj2-1]>atmp);bj2--)
                sl[bj2]=sl[bj2-1]; sl[bj2]=atmp;
              atmp=su[bi];for(bj2=bi;(bj2>0)&&(su[bj2-1]>atmp);bj2--)
                su[bj2]=su[bj2-1]; su[bj2]=atmp;
            }
            damage_ct_lamref=sl[2];damage_ct_duref=su[2];
          }
          damage_ct_on=1;damage_ct_refused=0;damage_ct_step=0;
          damage_ct_it=0;damage_ct_ncommit=0;
          dtheta=damage_ct_kappa*damage_ct_ds/(*tper);
          if(dtheta<=0.) dtheta=damage_ct_dt[ah];
          dthetaref=dtheta;
          printf("[DAMAGE CT] ARMED.  control point: element %" ITGFORMAT
                 " ip %" ITGFORMAT "; m=(%.6f,%.6f,%.6f); |m.delta_c|=%.6e; "
                 "c_lambda=%.6e; support dofs %" ITGFORMAT " of 18 free; "
                 "kappa=%.6e (max/min=%.4f over 5 intervals); ds0=%.6e; "
                 "ds_min=%.6e ds_max=%.6e; tol_c=%.6e; lambda_c=%.12e; "
                 "dtheta=kappa*ds/tper=%.6e; healthy refs dlambda=%.6e "
                 "|du|inf=%.6e.  lambda now OWNS the boundary for the rest "
                 "of the step; every ending is PARTIAL.%s",
                 abest+1,aip+1,damage_ct_m[0],damage_ct_m[1],damage_ct_m[2],
                 damage_ct_m[0]*adl[0]+damage_ct_m[1]*adl[1]
                 +damage_ct_m[2]*adl[2],
                 damage_ct_clam,damage_ct_nsupp,damage_ct_kappa,
                 damage_ct_tau,damage_ct_ds0,damage_ct_dsmin,
                 damage_ct_dsmax,damage_ct_tolc,damage_ct_lamc,dtheta,
                 damage_ct_lamref,damage_ct_duref,"\n");
          fflush(stdout);
        }
      }
      }
      if(damage_ct_on==0) damage_ct_refused=1;
    }

    /* [DAMAGE CT] once armed, lambda owns the boundary for the rest of the
       step.  xbounold is the value at the START of the step, so this is an
       absolute step fraction, and it is the same ramp CCX_DAMAGE_PATH mode 1
       verified against tempload at max|diff|=0.  theta stays monotone
       pseudo-time and carries dtime only. */
    if(damage_ct_on==1){
      for(k=0;k<*nboun;k++){
        xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*damage_ct_lam;
      }
      /* dtime = kappa*ds is re-established for every armed attempt, so a
         reduced ds reduces the viscous time by the same factor and the
         total viscous time across a given control advance is invariant. */
      if(damage_ct_kappa>0.){
        dtheta=damage_ct_kappa*damage_ct_ds/(*tper);
        if(dtheta>0.) dthetaref=dtheta;
      }
    }

    /* Path control: give the load factor an identity of its own.

       tempload builds xbounact from reltime=theta+dtheta, and theta only
       ever advances, so lambda cannot be reduced.  That is the one
       capability missing when a deletion releases energy and the
       neighbouring equilibrium sits at a LOWER load.  xbounold holds the
       value at the start of the step (it is refreshed only after the
       increment loop), so xbounold+(xboun-xbounold)*lambda is an
       absolute step fraction and lambda is free to move either way.

       Mode 1 only measures: it rebuilds the ramp and reports the largest
       deviation from what tempload produced, which has to be zero before
       anything is allowed to rely on the equivalence.  Amplitudes would
       break it, and this reports that rather than assuming it. */

    if(damage_path_on>0){
      damage_path_lam=theta+dtheta;

      /* Diagnostic descent.  After a deletion transaction has been
         rolled back, the question is whether an equilibrium exists
         at all for the post-deletion topology, or only at a LOWER
         load.  theta cannot answer it: it only advances.  Here the
         retry is placed below the committed load factor, keeping
         theta as a monotone counter so dtime stays positive for the
         viscous update.  If Newton converges there, the stall is a
         limit point and path control is the answer; if it converges
         nowhere, it is not. */

      /* Arm on ANY repeated attempt, not only a rolled-back deletion.
         The first version armed on the deletion rollback alone and
         fired exactly once in a 533-increment run, never reaching the
         stall, which fails as slow convergence with no deletion in
         that increment at all. */
      /* Controlled descent, one accepted increment at a time.

         The first version dropped lambda by up to 80% inside a single
         retry.  That does not test whether an equilibrium exists lower
         down; it asks Newton to jump from a converged field at one load
         to a distant field at another, and it destroyed the state -
         displacement increment 229 in a specimen of size 4, residual
         3.4e8, damage saturated in 9195 elements.  A descent has to be
         walked, with each step accepted, exactly like a loading path. */

      if((damage_path_on>=2)&&(damage_path_desc>0)){
        damage_path_lam=damage_path_lamcom*(1.-damage_path_drop);
        if(damage_path_lam<1.e-6) damage_path_lam=1.e-6;

        /* The descent is armed when dtheta has already collapsed to
           the floor, so without this it fails the minimum-size test
           immediately and never gets to walk - three attempts and
           out.  With lambda decoupled, dtheta only carries time for
           the viscous update; the load step is the lambda drop.
           Give the descent a workable step and a fresh set of
           attempts, bounded so a descent that never converges cannot
           spin. */
        damage_path_att++;
        if(damage_path_att>4*damage_path_nstep){
          damage_path_desc=0;
          printf("[DAMAGE PATH] descent abandoned after %" ITGFORMAT
                 " attempts without an accepted step\n",
                 damage_path_att);
          fflush(stdout);
        }else{
          if(dtheta<0.5*dthetaref) dtheta=0.5*dthetaref;
          icutb=0;
        }
        printf("[DAMAGE PATH] inc=%" ITGFORMAT " descent step %"
               ITGFORMAT " left: lambda %.6f -> %.6f\n",
               iinc,damage_path_desc,damage_path_lamcom,damage_path_lam);
        fflush(stdout);
      }
      damage_path_dev=0.;
      for(k=0;k<*nboun;k++){
        damage_path_ref=xbounold[k]+(xboun[k]-xbounold[k])*damage_path_lam;
        if(fabs(damage_path_ref-xbounact[k])>damage_path_dev)
          damage_path_dev=fabs(damage_path_ref-xbounact[k]);
        if(damage_path_on>=2) xbounact[k]=damage_path_ref;
      }
      if(damage_path_dev>damage_path_devmax){
        damage_path_devmax=damage_path_dev;
        printf("[DAMAGE PATH] inc=%" ITGFORMAT " lambda=%.6f "
               "max|ramp-tempload|=%.6e (new maximum)\n",
               iinc,damage_path_lam,damage_path_dev);
        fflush(stdout);
      }
    }

    for(i=0;i<3;i++){
      cam[i]=0.;}
    for(i=3;i<5;i++){
      cam[i]=0.5;}
    if(*ithermal>1){
      radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,
		  ac,bc,nload,sideload,nelemload,xloadact,lakon,ipiv,ntmat_,
		  vold,
		  shcon,nshcon,ipkon,kon,co,
		  kontri,&ntri,nloadtr,tarea,tenv,physcon,erad,&adview,&auview,
		  nflow,ikboun,xbounact,nboun,ithermal,&iinc,&iit,
		  cs,mcs,inocs,&ntrit,nk,fenv,istep,&dtime,ttime,&time,ilboun,
		  ikforc,ilforc,xforcact,nforc,cam,ielmat,&nteq,prop,ielprop,
		  nactdog,nacteq,nodeboun,ndirboun,network,
		  rhcon,nrhcon,ipobody,ibody,xbodyact,nbody,iviewfile,jobnamef,
		  ctrl,xloadold,&reltime,nmethod,set,mi,istartset,iendset,
		  ialset,nset,
		  ineighe,nmpc,nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,
		  iamload,
		  jqrad,irowrad,&nzsrad,icolrad,ne,iaxial,qa,cocon,ncocon,
		  iponoeln,
		  inoeln,nprop,amname,namta,amta,iexpl);
             
      /* check whether network iterations converged */

      if(qa[2]>0){
	checkdivergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod, 
			kode,filab,een,t1act,&time,epn,ielmat,matname,enern, 
			xstaten,nstate_,istep,&iinc,iperturb,ener,mi,output,
			ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
			ielorien,norien,description,sti,&icutb,&iit,&dtime,qa,
			vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,&icntrl,
			&theta,&dtheta,veold,vini,idrct,tper,&istab,tmax, 
			nactdof,b,tmin,ctrl,amta,namta,itpamp,inext,&dthetaref,
			&itp,&jprint,jout,&uncoupled,t1,&iitterm,nelemload,
			nload,nodeboun,nboun,itg,ndirboun,&deltmx,&iflagact,
			set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,
			mortar,nmat,ielprop,prop,&ialeatoric,&kscale,
			energy, &allwk,&energyref,&emax,&r_abs,&enetoll,
			energyini,
			&allwkini,&temax,&sizemaxinc,&ne0,&neini,&dampwk,
			&dampwkini,energystartstep);

	/* the divergence is flagged by icntrl!=0
	   icutb is reset to zero in order to generate
	   regular contact elements etc.. */

	icutb--;
      }
    }
      
    if(icfd==2){
      compfluidfem(&cof,&nkf,&ipkonf,&konf,&lakonf,nef,&sideface,
		   ifreestream,&nfreestream,isolidsurf,neighsolidsurf,
		   &nsolidsurf,iponoelf,inoelf,nshcon,shcon,nrhcon,rhcon,
		   &voldf,ntmat_,nodebounf,ndirbounf,&nbounf,&ipompcf,
		   &nodempcf,&nmpcf,&ikmpcf,&ilmpcf,ithermal,ikbounf,ilbounf,
		   &iturbulent,isolver,iexpl,ttime,&time,&dtime,nodeforc,
		   ndirforc,xforc,nforc,nelemloadf,sideloadf,
		   xloadf,&nloadf,xbody,ipobodyf,nbody,&ielmatf,matname,mi,
		   ncmat_,physcon,istep,&iinc,ibody,xloadold,xbounf,&coefmpcf,
		   nmethod,xforcold,xforcact,iamforc,iamloadf,xbodyold,xbodyact,
		   t1old,t1,t1act,iamt1,amta,namta,nam,ampli,xbounold,xbounact,
		   iambounf,itg,&ntg,amname,t0,&nelemface,&nface,cocon,ncocon,
		   xloadact,tper,jmax,jout,set,nset,istartset,iendset,ialset,
		   prset,prlab,nprint,trab,inotr,ntrans,filab,&labmpc,sti,
		   norien,orab,jobnamef,tieset,ntie,mcs,ics,cs,nkon,&mpcfree,
		   &memmpc_,&fmpc,nef,&inomat,qfx,kode,ipface,ielprop,prop,
		   orname,tincf,&ifreesurface,&nkftot,ielorienf,nelold,nkold,
		   nknew,nelnew);

      for(i=0;i<nkftot;i++){
	for(j=0;j<mt;j++){
	  vold[mt*(nkold[i]-1)+j]=voldf[mt*i+j];
	}
      }
    }
      
    if(icascade==2){
      memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
      RENEW(nodempc,ITG,3*memmpcref_);
      isiz=3*memmpcref_;cpyparitg(nodempc,nodempcref,&isiz,&num_cpus);
      RENEW(coefmpc,double,memmpcref_);
      isiz=memmpcref_;cpypardou(coefmpc,coefmpcref,&isiz,&num_cpus);
    }

    /* generating contact elements */
      
    if((ncont!=0)&&(*mortar<=1)&&

       /*       for purely thermal calculations: determine contact integration
		points only at the start of a step */

       ((*ithermal!=2)||(iit==-1))){

      *ne=ne0;*nkon=nkon0;

      /* at start of new increment: 
	 - copy state variables (node-to-face)
	 - determine slave integration points (face-to-face)
	 - interpolate state variables (face-to-face) */

      if(icutb==0){
	if(*mortar==1){

	  if(*nstate_!=0){
	    if(maxprevcontel!=0){
	      if(iit!=-1){
		NNEW(islavsurfold,ITG,2**ifacecount+2);
		NNEW(pslavsurfold,double,3**nintpoint);
		isiz=2**ifacecount+2;
		cpyparitg(islavsurfold,islavsurf,&isiz,&num_cpus);
		isiz=3**nintpoint;
		cpypardou(pslavsurfold,pslavsurf,&isiz,&num_cpus);
	      }
	    }
	  }

	  *nintpoint=0;

	  /* determine the location of the slave integration
	     points */

	  precontact(&ncont,ntie,tieset,nset,set,istartset,
                     iendset,ialset,itietri,lakon,ipkon,kon,koncont,ne,
                     cg,straight,co,vold,istep,&iinc,&iit,itiefac,
                     islavsurf,islavnode,imastnode,nslavnode,nmastnode,
                     imastop,mi,ipe,ime,tietol,
		     nintpoint,&pslavsurf,xmastnor,cs,mcs,ics,clearini,
                     nslavs);
		  
	  /* changing the dimension of element-related fields */
		  
	  RENEW(kon,ITG,*nkon+22**nintpoint);
	  RENEW(springarea,double,2**nintpoint);
	  RENEW(pmastsurf,double,6**nintpoint);
		  
	  if(*nener==1){
	    RENEW(ener,double,mi[0]*(*ne+*nintpoint)*2);

	    /* setting the entries for the contact energy to zero */

	    DOUMEMSET(ener,2*mi[0]**ne,2*mi[0]*(*ne+*nintpoint),0.);

	  }
	  RENEW(ipkon,ITG,*ne+*nintpoint);
	  RENEW(lakon,char,8*(*ne+*nintpoint));
		  
	  if(*norien>0){
	    RENEW(ielorien,ITG,mi[2]*(*ne+*nintpoint));
	    ITGMEMSET(ielorien,mi[2]**ne,mi[2]*(*ne+*nintpoint),0);
	  }
	  RENEW(ielmat,ITG,mi[2]*(*ne+*nintpoint));
	  isiz=mi[2]**nintpoint;
	  ITGMEMSET(ielmat,mi[2]**ne,mi[2]*(*ne+*nintpoint),1);

	  /* interpolating the state variables */

	  if(*nstate_!=0){
	    if(maxprevcontel!=0){
	      RENEW(xstateini,double,
		    *nstate_*mi[0]*(ne0+maxprevcontel));
	      isiz=*nstate_*mi[0]*maxprevcontel;
	      cpypardou(&xstateini[*nstate_*mi[0]*ne0],
			&xstate[*nstate_*mi[0]*ne0],&isiz,&num_cpus);
	    }
		      
	    RENEW(xstate,double,*nstate_*mi[0]*(ne0+*nintpoint));
	    isiz=*nstate_*mi[0]**nintpoint;
	    DOUMEMSET(xstate,*nstate_*mi[0]*ne0,
		      *nstate_*mi[0]*(ne0+*nintpoint),0.);
		      
	    if((*nintpoint>0)&&(maxprevcontel>0)){
			  
	      /* interpolation of xstate */
			  
	      interpolatestatemain(ne,ipkon,kon,lakon,
				   &ne0,mi,xstate,pslavsurf,nstate_,
				   xstateini,islavsurf,islavsurfold,
				   pslavsurfold,tieset,ntie,itiefac);
			  
	    }

	    if(maxprevcontel!=0){
	      SFREE(islavsurfold);SFREE(pslavsurfold);
	    }

	    maxprevcontel=*nintpoint;

	    RENEW(xstateini,double,*nstate_*mi[0]*(ne0+*nintpoint));
	    isiz=*nstate_*mi[0]*(ne0+*nintpoint);
	    cpypardou(xstateini,xstate,&isiz,&num_cpus);
	  }

	}

	/* set the contact spring energy to zero at the start of
           an increment. The friction energy is summed in energy[3]
           based on energyini[3] */
	
	if(*nener==1){
	  if((*mortar!=1)||(ncont==0)){
	    DOUMEMSET(enerini,2*mi[0]**ne,2*mi[0]*(*ne+*nslavs),0.);
	  }else{
	    RENEW(enerini,double,2*mi[0]*(*ne+*nintpoint));
	    DOUMEMSET(enerini,2*mi[0]**ne,2*mi[0]*(*ne+*nintpoint),0.);
	  }
	}
	
      }

      /* massless contact: calculate matrix Wb */
      
      if(*mortar==-1){

	if((masslesslinear==0)||(iinc==1)){
	  nzsw=5*9**nslavs;
	  
	  /* 5 = 1 slave + maximal 4 master, so 5 terms in the equation times
	     3 dofs = 15 terms; for each slave node 3 dofs, so 3 equations */ 
	  
	  NNEW(auw,double,nzsw);
	  NNEW(jqw,ITG,3**nslavs+1);
	  NNEW(iroww,ITG,nzsw);
	}
	if((masslesslinear>0)&&(iinc==1)){
	  NNEW(fullgmatrix,double,9**nslavs**nslavs);
	  NNEW(fullr,double,3**nslavs);
	}
      }

      if((*mortar!=-1)||(masslesslinear==0)||(iinc==1)){
	contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,straight,nkon,
		co,vold,ielmat,cs,elcon,istep,&iinc,&iit,ncmat_,ntmat_,
		&ne0,nmethod,
		iperturb,ikboun,nboun,mi,imastop,nslavnode,islavnode,
		islavsurf,
		itiefac,areaslav,iponoels,inoels,springarea,tietol,&reltime,
		imastnode,nmastnode,xmastnor,filab,mcs,ics,&nasym,
		xnoels,mortar,pslavsurf,pmastsurf,clearini,&theta,
		xstateini,xstate,nstate_,&icutb,&ialeatoric,jobnamef,
		&alea,auw,jqw,iroww,&nzsw);
      }
   
      /* check whether, for a dynamic calculation, contact damping 
	 is involved */

      if(*nmethod==4){
	if(*iexpl<=1){
	  if(idampingwithoutcontact==0){
	    for(i=0;i<*ne;i++){
	      if(ipkon[i]<0) continue;
	      if(*ncmat_>=5){
		if(strcmp1(&lakon[i*8],"ES")==0){
		  if(strcmp1(&lakon[i*8+6],"C")==0){
		    imat=ielmat[i*mi[2]];
		    if(elcon[(*ncmat_+1)**ntmat_*(imat-1)+4]>0.){
		      idamping=1;break;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
	  
      if(*iexpl<=1) printf(" Number of contact spring elements=%"
			   ITGFORMAT "\n\n",*ne-ne0);
            
      /* carlo start */

      /* dynamic time step estimation for explicit dynamics under penalty 
	 contact(CMT) start */

      if((*iexpl>1)&&(*mortar!=-1)){
	      
	if((*ne-ne0)<ncontacts){

	  /* number of contact elements has decreased */
	  
	  ncontacts=*ne-ne0;
	  inccontact=0;
	}  
	else if((*ne-ne0)>ncontacts)  {

	  /* number of contact elements has increased */
	  
	  RENEW(smscale,double,*ne);
		  
	  FORTRAN(calcstabletimeinccont,(ne,lakon,kon,ipkon,mi,ielmat,elcon,
					 mortar,adb,alpha,nactdof,springarea,
					 &ne0,ntmat_,ncmat_,&dtcont,smscale,
					 &dtset,&mscalmethod));
	  if(dtcont<dtvol){
	    dtmin=dtcont;
	  }else{
	    dtmin=dtvol;
	  }
	  
	  if(dtmin>(*tmax*(*tper))){
	    dtime=*tmax*(*tper);}
	  else if(dtmin<dtset){
	    dtime=dtset;}
	  else {
	    dtime=dtmin;
	  }
	  
	  dtheta=(dtime)/(*tper);
	  reltime=theta+dtheta;
	  time=reltime**tper;
	  dthetaref=dtheta;
	  printf(" SELECTED time increment (based on contact):%e\n\n",dtime);

	  ncontacts=*ne-ne0; 
	  inccontact=0;
	}else if((inccontact==500)&&(ncontacts==0)){

          /* no contact elements (either no contact or massless contact) */

	  if(dtvol>(*tmax*(*tper))){
	    dtime=*tmax*(*tper);
	  }else if(dtvol<dtset){
	    dtime=dtset;
	  }else{
	    dtime=dtvol;
	  }
	  dtheta=(dtime)/(*tper);
	  reltime=theta+dtheta;
	  time=reltime**tper;
	  dthetaref=dtheta;
	  printf(" SELECTED time increment (based on contact):%e\n\n",*tinc);

	  dtcont=1.e30;
	}
	inccontact++;
      }
	  
      /* CMT end */

    }
      
    /*  updating the nonlinear mpc's (also affects the boundary
	conditions through the nonhomogeneous part of the mpc's) */
      
    FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
		       nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
		       &maxlenmpc,ikmpc,ilmpc,&icascade,
		       kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,
		       &iit,&idiscon,&ncont,trab,ntrans,ithermal,mi,&kchdep));
      
    if(icascade==2){
      isiz=3*memmpc_;cpyparitg(nodempcref,nodempc,&isiz,&num_cpus);
      isiz=memmpc_;cpypardou(coefmpcref,coefmpc,&isiz,&num_cpus);
    }

    /* recalculating the matrix structure; only needed if:
       1) MPC's are cascaded
       2) contact occurs in an implicit calculation 
          (in penalty explicit no matrices are needed, in massless
           explicit no contact elements are generated) */
      
    if((icascade>0)||((ncont!=0)&&(*iexpl<=1)))
      remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		 &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,ilboun,
		 labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		 kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		 neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		 &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar,
		 typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
		 nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus,
		 ielmat,matname);

    /* invert nactdof (not for dynamic explicit calculations) */

    if(*iexpl<=1){
      SFREE(nactdofinv);
      NNEW(nactdofinv,ITG,mt**nk);
      MNEW(nodorig,ITG,*nk);
      FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			     ipkon,lakon,kon,ne));
      SFREE(nodorig);
    }
      
    /* check whether the forced displacements changed; if so, and
       if the procedure is static, the first iteration has to be
       purely linear elastic, in order to get an equilibrium
       displacement field; otherwise huge (maybe nonelastic)
       stresses may occur, jeopardizing convergence */
      
    iforbou=0;
      
    /* only for iinc=1 a linearized calculation is performed, since
       for iinc>1 a reasonable displacement field is predicted by using the
       initial velocity field at the end of the last increment */
      
    if((iinc==1)&&(*ithermal<2)&&((*nmethod!=4)||(*mortar==-1))){
      dev=0.;
      for(k=0;k<*nboun;++k){
	err=fabs(xbounact[k]-xbounini[k]);
	if(err>dev){dev=err;}
      }
      if(dev>1.e-5) iforbou=1;
    }
    if((*mortar==-1)&&(iforbou==1)){
      printf(" *ERROR in nonlingeo: nonzero boundary conditions are not allowed\n");
      printf("        in combination with massless contact\n\n");
      FORTRAN(stop,());
    }
      
    /* prediction of the kinematic vectors  */
      
    NNEW(v,double,mt**nk);
    
    /* for massless contact there is no need for prediction,
       since scheme is on velocity level */
      
    /* ---- CCX_PATHFOLLOW: the stock predictor is the wrong predictor ---
       prediction() extrapolates v = vold + dtime*veold, i.e. it repeats
       the previous increment's displacement change.  That is right when
       the step time IS the load parameter.  It is wrong here for two
       independent reasons, both measured on the target at increment 60:

       1. It moves the CONTROL COORDINATE before the corrector runs.  The
          extrapolation advanced phi to 1.65e-04 against a target of
          1.0e-05, so every corrector opened by having to undo 94% of the
          predictor, and a cutback made that worse rather than better
          because it shrinks dphi faster than it shrinks the predictor.

       2. It destroys the capture of f_hat.  f_hat is obtained as
          b/dlamjump at the first iteration, which is -dR/dlambda only if
          the ONLY thing that moved since the committed state is the
          prescribed pattern.  With an extrapolated u, b also contains
          K*du_pred - which is designed to cancel the load term - so the
          quotient is a difference of two nearly equal quantities.
          Measured against an honest finite difference of the residual at
          fixed u, the frozen f_hat came out at cos=+0.983 but
          |f_hat|/|q| = 1.43.  A scale error s there leaves the
          DISPLACEMENT correction right and the LOAD FACTOR step a factor
          1/s short, so the prescribed dofs end up where the free dofs did
          not assume: the constraint row converged to 1e-19 while the
          equilibrium row stalled with lambda marching down 6.7e-05 per
          iteration.  That is exactly the trace that was measured.

       prediction() skips the extrapolation when idiscon != 0 and resets
       the flag itself, so suppressing it is one assignment.  The method
       keeps its OWN predictor: the previous accepted dlambda, applied to
       lambda alone, which is the standard arc-length predictor. */

    if((pf_on==1)&&(pf_engaged==1)) idiscon=1;

    if (*mortar==-1){
      memcpy(&v[0],&vold[0],sizeof(double)*mt **nk);
    }else{
      prediction(uam,nmethod,&bet,&gam,&dtime,ithermal,nk,veold,accold,v,
		 &iinc,&idiscon,vold,nactdof,mi,&num_cpus);
    }
      
    MNEW(fn,double,mt**nk);
    NNEW(stx,double,6*mi[0]**ne);
      
    /* determining the internal forces at the start of the increment
	 
       for a static calculation with increased forced displacements
       the linear strains are calculated corresponding to
	 
       the displacements at the end of the previous increment, extrapolated
       if appropriate (for nondispersive media) +
       the forced displacements at the end of the present increment +
       the temperatures at the end of the present increment (this sum is
       v) -
       the displacements at the end of the previous increment (this is vold)
	 
       these linear strains are converted in stresses by multiplication
       with the tangent element stiffness matrix and converted into nodal
       forces. 
	 
       this boils down to the fact that the effect of forced displacements
       should be handled in a purely linear way at the
       start of a new increment, in order to speed up the convergence and
       (for dissipative media) guarantee smooth loading within the increment.
	 
       for all other cases the nodal force calculation is based on
       the true stresses derived from the appropriate strain tensor taking
       into account the extrapolated displacements at the end of the 
       previous increment + the forced displacements and the temperatures
       at the end of the present increment */
      
    iout=-1;
    if(istrainfree==1) iout=-2;
    iperturb_sav[0]=iperturb[0];
    iperturb_sav[1]=iperturb[1];
      
    /* first iteration in first increment: elastic tangent */
      
    if((*nmethod!=4)&&(iforbou==1)){
	  
      ielas=1;
	  
      iperturb[0]=-1;
      iperturb[1]=0;
	  
      isiz=neq[1];cpypardou(b,f,&isiz,&num_cpus);
      if(ne1d2d==1)NNEW(inum,ITG,*nk);
      results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	      elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	      ielorien,norien,orab,ntmat_,t1ini,t1act,ithermal,
	      prestr,iprestr,filab,eme,emn,een,iperturb,
	      f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	      ndirboun,xbounact,nboun,ipompc,
	      nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
	      &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	      xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	      &icmd, ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	      emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	      iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	      fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	      &reltime,&ne0,thicke,shcon,nshcon,
	      sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	      mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	      islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
	      inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	      itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	      islavquadel,aut,irowt,jqt,&mortartrafoflag,
	      &intscheme,physcon,dam,damn,iponoel);
      iperturb[0]=0;if(ne1d2d==1)SFREE(inum);
	  
      /* check whether any displacements or temperatures are changed
	 in the new increment */
	  
      for(k=0;k<neq[1];++k){
	f[k]=f[k]+b[k];}
	  
    }
    else{

      if(*mortar!=-1){
	if(ne1d2d==1)NNEW(inum,ITG,*nk);
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		ndirboun,xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		&reltime,&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
		inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavquadel,aut,irowt,jqt,&mortartrafoflag,
		&intscheme,physcon,dam,damn,iponoel);
	if(ne1d2d==1)SFREE(inum);
	  
	isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
	  
	if(*ithermal!=2){
	  isiz=6*mi[0]*ne0;	    
	  cpypardou(sti,stx,&isiz,&num_cpus);
	}
      }  
    }
      
    ielas=0;
    iout=0;

    /* ---- [DAMAGE RELEASE] report ------------------------------------
       f now holds f_int on the NEW topology, at a displacement state that is
       bit-identical to the converged u*: the topology site set idiscon=1
       before looping back, and prediction() with idiscon!=0 copies vold into
       v with no extrapolation (prediction.c:99).  The difference below is
       therefore the released internal force and nothing else.

       Self-checks printed with the numbers, because a probe that is trusted
       without them is worse than none:
         - action=reuse-sparse-graph means nactdof did not change, so
           removed MUST be exactly zero;
         - anom counts DOF active AFTER but not BEFORE, which a deletion
           cannot produce - nonzero means the mapping is wrong;
         - iforbou=1 would mean f had a boundary term added to it
           (f[k]+=b[k] on the other branch above) and the reading is void. */
    if((damage_release_probe)&&(damage_release_armed)&&
       (damage_frel!=NULL)&&(damage_ract!=NULL)){
      ITG ri,rj,rk,rns=0,rds=0,rnr=0,rdr=0,ranom=0;
      double dfv,dfa,smax=0.,sl1=0.,sl2=0.,rmax=0.,rl1=0.,rl2=0.;
      for(ri=0;ri<*nk;ri++){
        for(rj=0;rj<mt;rj++){
          rk=nactdof[mt*ri+rj];
          dfa=(rk>0)?f[rk-1]:0.;
          if(damage_ract[mt*ri+rj]){
            dfv=damage_frel[mt*ri+rj]-dfa;
            if(rk>0){
              sl1+=fabs(dfv); sl2+=dfv*dfv;
              if(fabs(dfv)>smax){smax=fabs(dfv);rns=ri+1;rds=rj;}
            }else{
              rl1+=fabs(dfv); rl2+=dfv*dfv;
              if(fabs(dfv)>rmax){rmax=fabs(dfv);rnr=ri+1;rdr=rj;}
            }
          }else if(rk>0){
            ranom++;
          }
        }
      }
      sl2=sqrt(sl2); rl2=sqrt(rl2);
      printf("[DAMAGE RELEASE] inc=%" ITGFORMAT " pass=%" ITGFORMAT
             " time=%.12e action=%s dt=%.6e%s"
             "   surv:    dF_max=%.6e node=%" ITGFORMAT " dof=%" ITGFORMAT
             " dF_l1=%.6e dF_l2=%.6e%s"
             "   removed: dF_max=%.6e node=%" ITGFORMAT " dof=%" ITGFORMAT
             " dF_l1=%.6e dF_l2=%.6e%s"
             "   qa=%.6e qam=%.6e  dF_surv_max/qam=%.6e  dF_max/qam=%.6e%s"
             "   deleted: terminal=%" ITGFORMAT " deadall+deadsole=%" ITGFORMAT
             " islands=%" ITGFORMAT " cohfacets=%" ITGFORMAT "%s"
             "   selfcheck: anom=%" ITGFORMAT " iforbou=%" ITGFORMAT
             " removed_must_be_zero=%s%s",
             iinc,damage_release_pass,theta**tper,
             damage_release_rebuild?"remastruct":"reuse-sparse-graph",
             damage_release_dt,"\n",
             smax,rns,rds,sl1,sl2,"\n",
             rmax,rnr,rdr,rl1,rl2,"\n",
             damage_release_qa,damage_release_qam,
             (damage_release_qam>0.)?smax/damage_release_qam:-1.,
             (damage_release_qam>0.)?
               ((smax>rmax?smax:rmax)/damage_release_qam):-1.,"\n",
             damage_release_nterm,damage_release_nother,
             damage_release_nisl,damage_release_ncoh,"\n",
             ranom,damage_release_iforbou,
             damage_release_rebuild?"no":"YES","\n");
      fflush(stdout);
      damage_release_armed=0;
    }

      
    SFREE(fn);SFREE(v);
    if((*ithermal!=3)||(ncont==0)||(*mortar!=1)||(*ncmat_<11)) SFREE(stx);
      
    /***************************************************************/
    /* iteration counter and start of the loop over the iterations */
    /***************************************************************/

    if(*mortar>1){	  
      NNEW(bhat,double,neq[1]);
      NNEW(islavactdof,ITG,neq[1]);
    } 
      
    iit=1;

    /* change due to previous checkdivergence routine */

    if(icntrl!=0) icutb++;

    ctrl[0]=i0ref;ctrl[1]=irref;ctrl[3]=icref;
    damage_slow_extended=0;
    damage_slow_camprev1=1.e300;
    damage_slow_camprev2=1.e300;
    damage_slow_maxiters=DAMAGE_SLOW_NEWTON_MAX_ITERS;
    if((ITG)icref>damage_slow_maxiters)
      damage_slow_maxiters=(ITG)icref;
    if(*nmethod!=4)NNEW(resold,double,neq[1]);
    if(uncoupled){
      *ithermal=2;
      NNEW(iruc,ITG,nzs[1]-nzs[0]);
      for(k=0;k<nzs[1]-nzs[0];k++){
	iruc[k]=irow[k+nzs[0]]-neq[0];}
    }

    /* Second half of the rollback probe: PF-START reports |vold-vini| at
       the top of the attempt, this reports it at the last statement before
       the Newton loop.  Anything nonzero here is state that moved between
       the two, which the constraint would otherwise attribute to the
       increment. */

    if((pf_on==1)&&(getenv("CCX_PATHFOLLOW_ACCUMCHECK")!=NULL)){
      double pfd=0.,pft;
      ITG pfi,pfj,pfk;
      for(pfi=0;pfi<*nk;pfi++){
        for(pfj=1;pfj<mt;pfj++){
          pfk=nactdof[mt*pfi+pfj];
          if(pfk>0){
            pft=vold[mt*pfi+pfj]-vini[mt*pfi+pfj];
            pfd+=pft*pft;
          }
        }
      }
      printf("[PF-PRELOOP] inc=%" ITGFORMAT " icutb=%" ITGFORMAT
             " |vold-vini|=%.6e\n",iinc,icutb,sqrt(pfd));
      fflush(stdout);
    }

    while(icntrl==0){

#ifdef COMPANY
      FORTRAN(uiter,(&iit));
#endif	  

      /*  updating the nonlinear mpc's (also affects the boundary
	  conditions through the nonhomogeneous part of the mpc's) */

      if((iit!=1)||((uncoupled)&&(*ithermal==1))){

	printf(" iteration %" ITGFORMAT "\n\n",iit);

	/* restoring the distributed loading before adding the
	   friction heating */
	  
	if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
	  *nload=nloadref;
	  isiz=2**nload;cpyparitg(nelemload,nelemloadref,&isiz,&num_cpus);
	  if(*nam>0){
	    isiz=2**nload;cpyparitg(iamload,iamloadref,&isiz,&num_cpus);
	  }
	  memcpy(&sideload[0],&sideloadref[0],sizeof(char)*20**nload);
	}
	  
	FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
			  xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,
			  xbodyact,t1old,t1,t1act,iamt1,nk,amta,namta,nam,
			  ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
			  xbounold,xboun,xbounact,iamboun,nboun,nodeboun,
			  ndirboun,nodeforc,ndirforc,istep,&iinc,co,vold,itg,
			  &ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
			  ntrans,trab,inotr,veold,integerglob,doubleglob,
			  tieset,istartset,iendset,ialset,ntie,nmpc,ipompc,
			  ikmpc,ilmpc,nodempc,coefmpc,ipobody,iponoeln,inoeln,
			  ipkon,kon,ielprop,prop,ielmat,shcon,nshcon,rhcon,
			  nrhcon,cocon,ncocon,ntmat_,lakon,set,nset));

	for(i=0;i<3;i++){
	  cam[i]=0.;}
	for(i=3;i<5;i++){
	  cam[i]=0.5;}
	if(*ithermal>1){
	  radflowload(itg,ieg,&ntg,&ntr,adrad,aurad,bcr,ipivr,ac,bc,nload,
		      sideload,nelemload,xloadact,lakon,ipiv,ntmat_,vold,shcon,
		      nshcon,ipkon,kon,co,kontri,&ntri,nloadtr,tarea,tenv,
		      physcon,erad,&adview,&auview,nflow,ikboun,xbounact,nboun,
		      ithermal,&iinc,&iit,cs,mcs,inocs,&ntrit,nk,fenv,istep,
		      &dtime,ttime,&time,ilboun,ikforc,ilforc,xforcact,nforc,
		      cam,ielmat,&nteq,prop,ielprop,nactdog,nacteq,nodeboun,
		      ndirboun,network,rhcon,nrhcon,ipobody,ibody,xbodyact,
		      nbody,iviewfile,jobnamef,ctrl,xloadold,&reltime,nmethod,
		      set,mi,istartset,iendset,ialset,nset,ineighe,nmpc,
		      nodempc,ipompc,coefmpc,labmpc,&iemchange,nam,iamload,
		      jqrad,irowrad,&nzsrad,icolrad,ne,iaxial,qa,cocon,ncocon,
		      iponoeln,inoeln,nprop,amname,namta,amta,iexpl);
             
	  /* check whether network iterations converged */

	  if(qa[2]>0){
	    checkdivergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod,kode,filab,
			    een,t1act,&time,epn,ielmat,matname,enern,xstaten,
			    nstate_,istep,&iinc,iperturb,ener,mi,output,
			    ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
			    ielorien,norien,description,sti,&icutb,&iit,&dtime,
			    qa,vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,
			    &icntrl,&theta,&dtheta,veold,vini,idrct,tper,
			    &istab,tmax,nactdof,b,tmin,ctrl,amta,namta,itpamp,
			    inext,&dthetaref,&itp,&jprint,jout,&uncoupled,t1,
			    &iitterm,nelemload,nload,nodeboun,nboun,itg,
			    ndirboun,&deltmx,&iflagact,set,nset,istartset,
			    iendset,ialset,emn,thicke,jobnamec,mortar,nmat,
			    ielprop,prop,&ialeatoric,&kscale,energy,&allwk,
			    &energyref,&emax,&r_abs,&enetoll,energyini,
			    &allwkini,&temax,&sizemaxinc,&ne0,&neini,&dampwk,
			    &dampwkini,energystartstep);
	    continue;
	  }
	}

	if(icascade==2){
	  memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
	  RENEW(nodempc,ITG,3*memmpcref_);
	  isiz=3*memmpcref_;cpyparitg(nodempc,nodempcref,&isiz,&num_cpus);
	  RENEW(coefmpc,double,memmpcref_);
	  isiz=memmpcref_;cpypardou(coefmpc,coefmpcref,&isiz,&num_cpus);
	}

	if((ncont!=0)&&(*mortar<=1)&&(ismallsliding==0)&&
	   /*           for node-to-face contact: freeze contact elements for
			iterations 8 and higher */
	   ((iit<=8)||(*mortar==1))&&
	   /*           for purely thermal calculations: freeze contact elements
			during complete step */
	   ((*ithermal!=2)||(iit==-1))){

	  neold=*ne;
	  *ne=ne0;*nkon=nkon0;
	  contact(&ncont,ntie,tieset,nset,set,istartset,iendset,
		  ialset,itietri,lakon,ipkon,kon,koncont,ne,cg,
		  straight,nkon,co,vold,ielmat,cs,elcon,istep,
		  &iinc,&iit,ncmat_,ntmat_,&ne0,
		  nmethod,iperturb,
		  ikboun,nboun,mi,imastop,nslavnode,islavnode,islavsurf,
		  itiefac,areaslav,iponoels,inoels,springarea,tietol,
		  &reltime,imastnode,nmastnode,xmastnor,
		  filab,mcs,ics,&nasym,xnoels,mortar,pslavsurf,pmastsurf,
		  clearini,&theta,xstateini,xstate,nstate_,&icutb,
		  &ialeatoric,jobnamef,&alea,auw,jqw,iroww,&nzsw);

	  /* check whether, for a dynamic calculation, contact damping is 
	     involved */
	      
	  if(*nmethod==4){
	    if(*iexpl<=1){
	      if(idampingwithoutcontact==0){
		for(i=0;i<*ne;i++){
		  if(ipkon[i]<0) continue;
		  if(*ncmat_>=5){
		    if(strcmp1(&lakon[i*8],"ES")==0){
		      if(strcmp1(&lakon[i*8+6],"C")==0){
			imat=ielmat[i*mi[2]];
			if(elcon[(*ncmat_+1)**ntmat_*(imat-1)+4]>0.){
			  idamping=1;break;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	      
	  if(*mortar==0){
	    if(*ne!=neold){iflagact=1;}
	  }else if(*mortar==1){
	    if(((*ne-ne0)<(neold-ne0)*(1.-delcon))||
	       ((*ne-ne0)>(neold-ne0)*(1.+delcon))){iflagact=1;}
	  }

	  printf(" Number of contact spring elements=%" ITGFORMAT "\n\n",
		 *ne-ne0);

	}
	  
	if(*ithermal==3){
	  for(k=0;k<*nk;++k){
	    t1act[k]=vold[mt*k];}
	}

	FORTRAN(nonlinmpc,(co,vold,ipompc,nodempc,coefmpc,labmpc,
			   nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,
			   &maxlenmpc,ikmpc,ilmpc,&icascade,
			   kon,ipkon,lakon,ne,&reltime,&newstep,xboun,fmpc,&iit,
			   &idiscon,&ncont,trab,ntrans,ithermal,mi,&kchdep));

	if(icascade==2){
	  isiz=3*memmpc_;cpyparitg(nodempcref,nodempc,&isiz,&num_cpus);
	  isiz=memmpc_;cpypardou(coefmpcref,coefmpc,&isiz,&num_cpus);
	}

	/* recalculating the matrix structure */

	/* for face-to-face contact (mortar=1) this is only done if
	   the dependent term in nonlinear MPC's changed */
	
	if((icascade>0)||(ncont!=0)){
	  if((*mortar!=1)||(kchdep==1)){
	    remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		       &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,ikboun,
		       ilboun,labmpc,nk,&memmpc_,&icascade,&maxlenmpc,
		       kon,ipkon,lakon,ne,nactdof,icol,jq,&irow,isolver,
		       neq,nzs,nmethod,&f,&fext,&b,&aux2,&fini,&fextini,
		       &adb,&aub,ithermal,iperturb,mass,mi,iexpl,mortar,
		       typeboun,&cv,&cvini,&iit,network,itiefac,&ne0,&nkon0,
		       nintpoint,islavsurf,pmastsurf,tieset,ntie,&num_cpus,
		       ielmat,matname);
	  }

	  /* invert nactdof */
	      
	  SFREE(nactdofinv);
	  NNEW(nactdofinv,ITG,mt**nk);
	  MNEW(nodorig,ITG,*nk);
	  FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
				 ipkon,lakon,kon,ne));
	  SFREE(nodorig);
	      
	  MNEW(v,double,mt**nk);
	  NNEW(stx,double,6*mi[0]**ne);
	  MNEW(fn,double,mt**nk);
      
	  isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
	  iout=-1;
	      
	  if(ne1d2d==1)NNEW(inum,ITG,*nk);
	  results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		  prestr,iprestr,filab,eme,emn,een,iperturb,
		  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		  ndirboun,xbounact,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		  ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		  xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		  nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		  &reltime,&ne0,thicke,shcon,nshcon,
		  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		  mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		  islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
		  inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
		  itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		  islavquadel,aut,irowt,jqt,&mortartrafoflag,
		  &intscheme,physcon,dam,damn,iponoel);
	  
	  isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
	      
	  if(*ithermal!=2){
	    isiz=6*mi[0]*ne0;	    
	    cpypardou(sti,stx,&isiz,&num_cpus);
	  }
	      
	  SFREE(v);SFREE(fn);if(ne1d2d==1)SFREE(inum);
	  if((*ithermal!=3)||(ncont==0)||(*mortar!=1)||(*ncmat_<11)) SFREE(stx);
	  iout=0;
	}    
      }
	  
      /* add friction heating  */
      
      if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
	nload_=*nload+2*(*ne-ne0);

	RENEW(nelemload,ITG,2*nload_);
	ITGMEMSET(nelemload,2**nload,2*nload_,0);
	if(*nam>0){
	  RENEW(iamload,ITG,2*nload_);
	  ITGMEMSET(iamload,2**nload,2*nload_,0);
	}
	RENEW(xloadact,double,2*nload_);
	DOUMEMSET(xloadact,2**nload,2*nload_,0.);
	RENEW(sideload,char,20*nload_);
	DMEMSET(sideload,20**nload,20*nload_,'\0');

	MNEW(idefload,ITG,nload_);
	ITGMEMSET(idefload,0,nload_,1);
	FORTRAN(frictionheating,(&ne0,ne,ipkon,lakon,ielmat,mi,elcon,ncmat_,
				 ntmat_,kon,islavsurf,pmastsurf,springarea,co,
				 vold,veold,pslavsurf,xloadact,nload,&nload_,
				 nelemload,iamload,idefload,sideload,stx,nam,
				 &time,ttime,matname,istep,&iinc));
	SFREE(idefload);SFREE(stx);
      }

      /* calculate the stiffness matrix for:
         - implicit calculations
         - linear massless explicit calculations in the first increment
         - nonlinear massless explicit calculations */
      
      if((*iexpl<=1)||((*mortar==-1)&&((masslesslinear==0)||(iinc==1)))){

	/* calculating the local stiffness matrix and external loading */

	NNEW(ad,double,neq[1]);
	NNEW(au,double,nzs[1]);

	if(*nmethod==4){
	  DOUMEMSET(fnext,0,mt**nk,0.);
	}

	/* UNSYM stage 1: route the bulk assembly through the existing
	   asymmetric storage while the constitutive tangent is still the
	   stock symmetric g(D)*C_ep.

	   mafillsmmain writes the symmetric part into the lower triangle
	   exactly as before; mafillsmasmain then mirrors it into the upper
	   half at offset nzs[2] and adds asymmetric element contributions.
	   With no contact elements present it only performs the mirror, so
	   the assembled operator is numerically identical to the symmetric
	   one and PARDISO mtype=1 must reproduce the symmetric answer to
	   the last digit.  That is the acceptance test for this stage: it
	   validates the plumbing before any damage term is added to it.

	   Gates: implicit mechanical static path, no contact (contact owns
	   nasym itself), progressive damage material present. */

	damage_unsym_active=0;
	if((damage_tangent_mode==2)&&(damage_de12_enabled)&&
	   (ncont==0)&&(*iexpl<=1)&&(*nmethod!=4)&&(*nmethod!=5)&&
	   (*ithermal<2)&&(*mortar!=-1)){
	  damage_unsym_active=1;
	  if(nasym==0){
	    nasym=1;
	    printf("[DAMAGE TANGENT UNSYM] bulk assembly switched to the "
		   "asymmetric path (nasym=1); PARDISO symbolic reuse is "
		   "available here since 2026-08-28 (mtype=1, structurally "
		   "symmetric) - set CCX_PARDISO_REUSE_SYMBOLIC=1\n");
	    fflush(stdout);
	  }
	}

	mafillsmmain(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
		     ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		     nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
		     nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
		     nmethod,ikmpc,ilmpc,ikboun,ilboun,
		     elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		     ielmat,ielorien,norien,orab,ntmat_,
		     t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
		     nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
		     xstiff,npmat_,&dtime,matname,mi,
		     ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
		     physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,
		     &coriolis,ibody,xloadold,&reltime,veold,springarea,nstate_,
		     xstateini,xstate,thicke,integerglob,doubleglob,
		     tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
		     pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,&kscale,
		     iponoeln,inoeln,network,ntrans,inotr,trab,smscale,
		     &mscalmethod,set,nset,islavquadel,aut,irowt,jqt,
		     &mortartrafoflag);
	//		     &nslavquadel);

	if(nasym==1){
	  RENEW(au,double,2*nzs[1]);
	  if(*nmethod==4){
	    RENEW(aub,double,2*nzs[1]);}
	  symmetryflag=2;
	  inputformat=1;

	  mafillsmasmain(co,nk,kon,ipkon,lakon,ne,nodeboun,
			 ndirboun,xbounact,nboun,
			 ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
			 nforc,nelemload,sideload,xloadact,nload,xbodyact,
			 ipobody,
			 nbody,cgr,ad,au,fext,nactdof,icol,jq,irow,neq,nzl,
			 nmethod,ikmpc,ilmpc,ikboun,ilboun,
			 elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
			 ielmat,ielorien,norien,orab,ntmat_,
			 t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
			 nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			 xstiff,npmat_,&dtime,matname,mi,
			 ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,
			 physcon,shcon,nshcon,cocon,ncocon,ttime,&time,istep,
			 &iinc,
			 &coriolis,ibody,xloadold,&reltime,veold,springarea,
			 nstate_,
			 xstateini,xstate,thicke,
			 integerglob,doubleglob,tieset,istartset,iendset,
			 ialset,ntie,&nasym,pslavsurf,pmastsurf,mortar,clearini,
			 ielprop,prop,&ne0,&kscale,iponoeln,inoeln,network,set,
			 nset);

	  /* UNSYM stage 2: the rank-1 damage correction

	         -sigma_eff (x) dD/d(eps)

	     is added on top of the symmetric operator that mafillsmasmain
	     has just mirrored into both halves of au.  Only integration
	     points that are actively softening carry a nonzero damjac, so
	     the extra work scales with the size of the process zone, not
	     with the model. */

	  if(damage_unsym_active==1){
	    damage_unsym_elems=0;
	    FORTRAN(mafilldamas,(co,kon,ipkon,lakon,&ne0,nactdof,jq,irow,
				 neq,nzs,au,ad,vold,mi,damage_damjac,nmpc,
				 &damage_unsym_elems,dam,damdamageini,
				 &damage_unsym_skip,&damage_unsym_adv,
				 &damage_unsym_hole,&damage_unsym_floor));
	    /* J-10: what the operator looks like on THIS iteration, not as a
	       running maximum.  A count that moves between iterations of a
	       SAME-LOAD solve means the operator is changing shape mid-Newton,
	       which no step size and no tolerance can repair. */
	    if((damage_unsym_census==1)&&(idamagereeq==1)){
	      printf("[DAMAGE TANGENT CENSUS] inc=%" ITGFORMAT " iit=%" ITGFORMAT
		     " assembled=%" ITGFORMAT " adv=%" ITGFORMAT
		     " gap=%" ITGFORMAT " hi=%" ITGFORMAT
		     " skip=%" ITGFORMAT "\n",
		     iinc,iit,damage_unsym_elems,damage_unsym_adv,
		     damage_unsym_hole,damage_unsym_floor,damage_unsym_skip);
	      fflush(stdout);
	    }
	    if(damage_unsym_adv>damage_unsym_advrep){
	      damage_unsym_advrep=damage_unsym_adv;
	      printf("[DAMAGE TANGENT HOLE] %" ITGFORMAT " element(s) with "
		     "ADVANCING damage carry no rank-1 term (new maximum); "
		     "%" ITGFORMAT " past initiation without one, "
		     "%" ITGFORMAT " assembled\n",
		     damage_unsym_adv,damage_unsym_skip,damage_unsym_elems);
	      fflush(stdout);
	    }

	    /* The number above conflates three populations and only one is a
	       defect.  dambase is the damage at the START OF THE INCREMENT, so
	       "advancing" there also covers every point that advanced earlier
	       in the increment and unloads elastically now - zero is correct
	       for those, and for points already on the residual-stiffness
	       floor, where dg/dD=0.  What is NOT correct is the band
	       D in [0.999, 1-gmin): resultsmech.f refuses the consistent
	       tangent at D>=0.999 while the floor only engages at 1-D<gmin,
	       so those elements carry g*C_ep with g in (1e-4,1e-3] - which
	       mpfd measures at rho=3443 with the WRONG SIGN on the dominant
	       entry already at D=0.994.  Report that band on its own so the
	       hypothesis can be falsified rather than argued. */
	    if(damage_unsym_hole>damage_unsym_holerep){
	      damage_unsym_holerep=damage_unsym_hole;
	      printf("[DAMAGE TANGENT GAP] inc=%" ITGFORMAT " %" ITGFORMAT
	    	 " element(s) in D=[0.999,1-gmin) carry g*C_ep with no rank-1 "
	    	 "term (new maximum); %" ITGFORMAT " on the residual-stiffness "
	    	 "floor where zero is correct; %" ITGFORMAT " advancing "
	    	 "without a term in total\n",
	    	 iinc,damage_unsym_hole,damage_unsym_floor,damage_unsym_adv);
	      fflush(stdout);
	    }
	    if(damage_unsym_elems>damage_unsym_report){
	      damage_unsym_report=damage_unsym_elems;
	      printf("[DAMAGE TANGENT UNSYM] rank-1 correction assembled "
		     "for %" ITGFORMAT " softening element(s) (new maximum)\n",
		     damage_unsym_elems);
	      fflush(stdout);
	    }
	  }
	}

	/* Per-node stiffness taken from the assembled diagonal.

	   This is the measure the topological rule in damdangle.f cannot
	   reach.  Node 501 keeps 3 of its original 8 elements and is only
	   the 63rd softest of 2387 nodes by mesh geometry - 62 softer nodes
	   cause no trouble at all - yet it is the one the solver throws
	   0.69 of a specimen length.  What separates it is the damage state
	   of the elements it still has, and ad[] carries g(D), so it reports
	   the stiffness the node actually has rather than the one its
	   connectivity suggests.

	   The reference is the node's OWN diagonal while the model was
	   intact, so the ratio is dimensionless and needs no median over the
	   mesh.  A constrained direction is excluded: a node held by a
	   boundary condition is never stranded. */

	/* The solver frees ad and au right after the factorisation, so a
	   probe that runs after results() would read an empty matrix -
	   which is exactly what the elastic control caught, reporting a
	   relative error of 1.0 because every coefficient read back was
	   zero.  Keep a copy while the operator still exists. */

	if((damage_fd_inc>0)&&(iinc>=damage_fd_inc)&&(iit>=damage_fd_it)){
	  if(damage_fd_ad==NULL){
	    NNEW(damage_fd_ad,double,neq[1]);
	    NNEW(damage_fd_au,double,(nasym+1)*nzs[1]);
	  }
	  memcpy(damage_fd_ad,ad,sizeof(double)*neq[1]);
	  memcpy(damage_fd_au,au,sizeof(double)*(nasym+1)*nzs[1]);
	}

	/* Addressed dump of one node.

	   Global counters have twice disagreed with each other here, so this
	   asks the question directly for a named node instead: every element
	   that touches it, its type, whether it is still assembled, and for a
	   cohesive element the damage at each integration point.  The
	   assembled diagonal is printed alongside, so the claim "this node has
	   lost its support" is backed by the operator rather than inferred
	   from connectivity. */

	if((damage_dump_node>0)&&(damage_dump_node<=*nk)&&
	   (iinc>=damage_dump_inc)){
	  damage_dump_n=damage_dump_node-1;
	  printf("[NODE DUMP] inc=%" ITGFORMAT " node %" ITGFORMAT
		 " at (%.4f, %.4f, %.4f)\n",
		 iinc,damage_dump_node,co[3*damage_dump_n],
		 co[3*damage_dump_n+1],co[3*damage_dump_n+2]);
	  for(idir=1;idir<=3;idir++){
	    k=nactdof[mt*damage_dump_n+idir];
	    if(k<=0){
	      printf("[NODE DUMP]   dof %" ITGFORMAT ": constrained\n",idir);
	    }else{
	      printf("[NODE DUMP]   dof %" ITGFORMAT ": ad=%.6e",idir,ad[k-1]);
	      if((damage_addiag0!=NULL)&&(damage_addiag0[damage_dump_n]>0.))
		printf("   intact reference %.6e   ratio %.4e",
		       damage_addiag0[damage_dump_n],
		       ad[k-1]/damage_addiag0[damage_dump_n]);
	      printf("\n");
	    }
	  }
	  damage_dump_nb=0;damage_dump_nu=0;
	  for(i=0;i<*ne;i++){
	    if(ipkon[i]==-1) continue;
	    damage_dump_alive=(ipkon[i]>=0)?1:0;
	    damage_dump_idx=(ipkon[i]>=0)?ipkon[i]:(-ipkon[i]-2);
	    if(damage_dump_idx<0) continue;
	    damage_dump_np=0;
	    if(lakon[8*i]=='C'){
	      damage_dump_np=(lakon[8*i+3]=='4')?4:0;
	    }else if(lakon[8*i]=='U'){
	      damage_dump_np=(ITG)((unsigned char)lakon[8*i+7]);
	      if((damage_dump_np<1)||(damage_dump_np>20)) damage_dump_np=0;
	    }
	    if(damage_dump_np<=0) continue;
	    damage_dump_hit=0;
	    for(j=0;j<damage_dump_np;j++){
	      if(kon[damage_dump_idx+j]-1==damage_dump_n) damage_dump_hit=1;
	    }
	    if(damage_dump_hit==0) continue;
	    if(lakon[8*i]=='C'){
	      if(damage_dump_alive) damage_dump_nb++;
	      printf("[NODE DUMP]   elem %-7" ITGFORMAT " C3D4  %-8s dam=%.6f\n",
		     i+1,damage_dump_alive?"alive":"DELETED",
		     (i<ne0)?dam[mi[0]*i]:-1.);
	    }else{
	      if(damage_dump_alive) damage_dump_nu++;
	      printf("[NODE DUMP]   elem %-7" ITGFORMAT " UC6   %-8s",
		     i+1,damage_dump_alive?"alive":"DELETED");
	      for(j=0;j<3;j++){
		damage_dump_dv=xstate[*nstate_*(mi[0]*i+j)+1];
		printf("  ip%" ITGFORMAT ": dvisc=%.6f g=%.6e",
		       j+1,damage_dump_dv,1.-damage_dump_dv);
	      }
	      printf("\n");
	    }
	  }
	  printf("[NODE DUMP]   live bulk=%" ITGFORMAT
		 "   live cohesive=%" ITGFORMAT "\n",
		 damage_dump_nb,damage_dump_nu);
	  /* The node's own 3x3 block of the assembled operator.
	
	     damage_addiag is the smallest AXIS-ALIGNED diagonal, which is only a
	     proxy for how compliant the node is: a node whose soft direction is
	     skew to x, y and z can carry three healthy diagonals while this block
	     has a small eigenvalue.  Three cohesive facets meeting at a node
	     define a plane whose normal is generally not an axis, so that is
	     exactly where the proxy should be weakest (E-79).  Printing the block
	     lets the eigenvalue be compared with the diagonal offline instead of
	     assumed.
	
	     Storage, from pardiso.c inputformat==3: the off-diagonal terms are a
	     full CSC of BOTH triangles, column by column, jq giving the 1-based
	     start of each column and irow the row.  A symmetric assembly stores
	     only the lower triangle, so it is mirrored. */
	  if((au!=NULL)&&(jq!=NULL)&&(irow!=NULL)){
	    ITG dmp_i,dmp_j,dmp_k,dmp_dof[3];
	    double dmp_blk[9];
	    for(dmp_i=0;dmp_i<9;dmp_i++) dmp_blk[dmp_i]=0.;
	    for(dmp_i=0;dmp_i<3;dmp_i++)
	      dmp_dof[dmp_i]=nactdof[mt*damage_dump_n+dmp_i+1];
	    for(dmp_j=0;dmp_j<3;dmp_j++){
	      if(dmp_dof[dmp_j]<=0) continue;
	      dmp_blk[dmp_j*3+dmp_j]=ad[dmp_dof[dmp_j]-1];
	      for(dmp_k=jq[dmp_dof[dmp_j]-1];dmp_k<=jq[dmp_dof[dmp_j]]-1;dmp_k++){
	        for(dmp_i=0;dmp_i<3;dmp_i++){
	          if(dmp_dof[dmp_i]<=0) continue;
	          if(irow[dmp_k-1]==dmp_dof[dmp_i]){
	            dmp_blk[dmp_i*3+dmp_j]=au[dmp_k-1];
	            /* au holds (nasym+1)*nzs[1] entries: the lower triangle in
	               [0,nzs) and, when nasym==1, its TRANSPOSE in [nzs,2nzs) -
	               the same layout the damping assembly uses a few hundred
	               lines below.  Reading only the first half produced a block
	               with a zero upper triangle, which is how this was caught. */
	            dmp_blk[dmp_j*3+dmp_i]=nasym?au[nzs[1]+dmp_k-1]:au[dmp_k-1];
	          }
	        }
	      }
	    }
	    printf("[NODE DUMP]   block3x3 nasym=%" ITGFORMAT,nasym);
	    for(dmp_i=0;dmp_i<9;dmp_i++) printf(" %.6e",dmp_blk[dmp_i]);
	    printf("\n");
	  }
	  fflush(stdout);
	}

	if(damage_stiff_probe>0){
	  /* [DAMSTATE] One owner for the judgement.  This block used to compute
	     the per-node diagonal, its intact reference and the AUTOSPC mask
	     inline; damstate.c now owns all three and the globals below are
	     views onto it, so every consumer - the displacement norm in
	     resultsini.c, the force norm here, the termination connectivity -
	     reads ONE decision instead of three separate ones.  Behaviour is
	     unchanged: damstate_update reproduces this loop exactly, including
	     that the minimum over the three dofs is a real minimum when a
	     negative diagonal is present, which its self test pins. */
	  if(damage_dstate.nk==0)
	    damstate_init(&damage_dstate,*nk,damage_spc_g,damage_spc_neg);
	  damstate_update(&damage_dstate,ad,nactdof,mt);
	  damage_addiag=damage_dstate.diag;
	  damage_addiag0=damage_dstate.diag0;
	  damage_addok=damage_dstate.ok;
	  if(damage_spc_g>0.){
	    damage_spc_mask=damage_dstate.dead;
	    damage_spc_nk=*nk;
	    damage_spc_count=damage_dstate.ndead;
	  }
	}

	iperturb[0]=iperturb_sav[0];
	iperturb[1]=iperturb_sav[1];

      }else{

	/* calculating the external loading 

	   This is only done once per increment. In reality, the
           external loading is a function of vold (specifically,
           the body forces and surface loading). This effect is
           neglected, since the increment size in dynamic explicit
           calculations is usually small */

	if((*mortar==-1)&&(masslesslinear==1)&&(iinc==2)){

	  /* check whether the distributed loading changes in this step 
	   (only for linear massless explicit dynamic calculations) */
	  
	  FORTRAN(checktempload,(iamload,nload,sideload,ibody,nbody,
				 &masslesslinear,&nloadrhs,&nbodyrhs,nam));

	  /* if no change: calculate the external force vector due to this
             loading only once at the start of the step */

	  if(masslesslinear==2){
	    NNEW(fextload,double,neq[1]);
	    nforcrhs=0;
	    rhsmain(co,nk,kon,ipkon,lakon,ne,
		    ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
		    &nforcrhs,nelemload,sideload,xloadact,nload,xbodyact,
		    ipobody,nbody,cgr,fextload,nactdof,&neq[1],
		    nmethod,ikmpc,ilmpc,
		    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
		    ielmat,ielorien,norien,orab,ntmat_,
		    t0,t1act,ithermal,iprestr,vold,iperturb,
		    iexpl,plicon,nplicon,plkcon,nplkcon,
		    npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
		    xbodyold,&reltime,veold,matname,mi,ikactmech,
		    &nactmech,ielprop,prop,sti,xstateini,xstate,nstate_,
		    ntrans,inotr,trab,fnext);
	  }
	}

	/* call to rhsmain in every increment: if the distributed
           loading does not change only point forces are taken into account */
	
	rhsmain(co,nk,kon,ipkon,lakon,ne,
	  	ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
	  	nforc,nelemload,sideload,xloadact,&nloadrhs,xbodyact,ipobody,
	  	&nbodyrhs,cgr,fext,nactdof,&neq[1],
	  	nmethod,ikmpc,ilmpc,
	  	elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	  	ielmat,ielorien,norien,orab,ntmat_,
	  	t0,t1act,ithermal,iprestr,vold,iperturb,
	  	iexpl,plicon,nplicon,plkcon,nplkcon,
	  	npmat_,ttime,&time,istep,&iinc,&dtime,physcon,ibody,
	  	xbodyold,&reltime,veold,matname,mi,ikactmech,
	  	&nactmech,ielprop,prop,sti,xstateini,xstate,nstate_,
	        ntrans,inotr,trab,fnext);
	//   for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}

	/* adding fextload due to distributed loading */

	if(masslesslinear==2){
	  for(i=0;i<neq[1];i++){fext[i]+=fextload[i];}
	}

      }
      

      /* calculating the damping matrix for implicit dynamic
         calculations */

      if((idamping==1)&&(*iexpl<=1)){

	/* Rayleigh damping */

	MNEW(adc,double,neq[1]);DOUMEMSET(adc,neq[0],neq[1],0.);
	for(k=0;k<neq[0];k++){
	  adc[k]=alpham*adb[k]+betam*ad[k];}
	if(nasym==0){
	  MNEW(auc,double,nzs[1]);DOUMEMSET(auc,nzs[0],nzs[1],0.);
	  for(k=0;k<nzs[0];k++){
	    auc[k]=alpham*aub[k]+betam*au[k];}
	}else{
	  NNEW(auc,double,2*nzs[1]);DOUMEMSET(auc,2*nzs[0],2*nzs[1],0.);
	  for(k=0;k<2*nzs[0];k++){
	    auc[k]=alpham*aub[k]+betam*au[k];}
	}
 
	/* dashpots and contact damping */

	FORTRAN(mafilldm,(co,nk,kon,ipkon,lakon,ne,nodeboun,
			  ndirboun,xbounact,nboun,
			  ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,
			  xforcact,
			  nforc,nelemload,sideload,xloadact,nload,xbodyact,
			  ipobody,nbody,cgr,
			  adc,auc,nactdof,icol,jq,irow,neq,nzl,nmethod,
			  ikmpc,ilmpc,ikboun,ilboun,
			  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
			  ielorien,norien,orab,ntmat_,
			  t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
			  nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
			  xstiff,npmat_,&dtime,matname,mi,ncmat_,
			  ttime,&time,istep,&iinc,ibody,clearini,mortar,
			  springarea,
			  pslavsurf,pmastsurf,&reltime,&nasym));
      }

      /* calculating the residual (RHS of equation system) */

      if(*mortar!=-1){
	calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux2,vold,
		     vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,
		     fini,islavnode,nslavnode,mortar,ntie,mi,
		     nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
		     &num_cpus);
      }else{
	NNEW(volddof,double,neq[0]);
	if(ncont!=0){NNEW(qb,double,neqtot);}
        massless(kslav,lslav,ktot,ltot,au,ad,auc,adc,jq,irow,neq,nzs,auw,jqw,
		 iroww,&nzsw,islavnode,nslavnode,nslavs,imastnode,nmastnode,
		 ntie,nactdof,mi,vold,volddof,veold,nk,fext,isolver,
		 &masslesslinear,co,springarea,&neqtot,qb,b,&dtime,aloc,fric,
		 iexpl,nener,ener,ne,&jqbi,&aubi,&irowbi,&jqib,&auib,&irowib,
		 &iclean,&iinc,fullgmatrix,fullr,alglob,&num_cpus,&ncont);
        if(masslesslinear==0){SFREE(ad);SFREE(au);} 
      }
      
      /*    for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}
	    for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}
	    for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	    for(k=0;k<neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k,ad[k]);}
	    for(k=0;k<nzs[1];++k){printf("au=%" ITGFORMAT ",%f\n",k,au[k]);}*/

      /* mortar contact */

      if(*mortar>1){
	
	/* trafo u -> util */
	
	premortar(nzs,&nzsc2,&auc2,&adc2,
		  &irowc2,&icolc2,&jqc2,&aubd,&irowbd,&jqbd,&aubdtil,
		  &irowbdtil,&jqbdtil,&aubdtil2,&irowbdtil2,&jqbdtil2,
		  &audd,&irowdd,&jqdd,&auddtil,&irowddtil,&jqddtil,
		  &auddtil2,&irowddtil2,&jqddtil2,&auddinv,&irowddinv,
		  &jqddinv,&jqtemp,&irowtemp,&icoltemp,nzstemp,&iit,
		  icol,irow,jq,ikboun,ilboun,ikmpc,ilmpc,
		  imastnode,nmastnode,co,nk,kon,ipkon,lakon,ne,stn,
		  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		  ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
		  iprestr,filab,eme,emn,een,iperturb,nactdof,&iout,qa,
		  vold,b,nodeboun,ndirboun,xbounact,xboun,nboun,ipompc,
		  nodempc,coefmpc,labmpc,nmpc,nmethod,neq,veold,accold,
		  &dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		  xstateini,xstiff,xstate,npmat_,matname,mi,&ielas,&icmd,
		  ncmat_,nstate_,stiini,vini,ener,enern,emeini,xstaten,
		  eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		  ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		  nelemload,nload,istep,&iinc,springarea,&reltime,&ne0,
		  xforc,nforc,thicke,shcon,nshcon,sideload,xload,xloadold,
		  &icfd,inomat,islavquadel,islavsurf,iponoels,inoels,
		  mortar,nslavnode,
		  islavnode,nslavs,ntie,aut,irowt,jqt,autinv,
		  irowtinv,jqtinv,tieset,
		  itiefac,&rhsi,au,ad,&f_cm,&f_cs,t1act,cam,&bet,&gam,epn,
		  xloadact,nodeforc,ndirforc,xforcact,xbodyact,ipobody,
		  nbody,cgr,nzl,sti,iexpl,mass,&buckling,&stiffness,
		  &intscheme,physcon,&coriolis,ibody,integerglob,
		  doubleglob,&nasym,&alpham,&betam,pslavsurf,
		  pmastsurf,clearini,ielprop,prop,islavact,cdn,&memmpc_,
		  &idamping,&iforbou,iperturb_sav,
		  itietri,cg,straight,koncont,energyini,energy,&kscale,
		  iponoeln,inoeln,nener,orname,network,typeboun,&num_cpus,
		  t0g,t1g,smscale,&mscalmethod,&nslavquadel,iponoel);
	
	/* calculating coupling matrices and embedding weak 
	   contact conditions */ 
      
	contactmortar(&ncont,ntie,tieset,nset,set,istartset,iendset,ialset,
		      itietri,lakon,ipkon,kon,koncont,ne,cg,straight,co,vold,
		      ielmat,elcon,istep,&iinc,&iit,ncmat_,ntmat_,&ne0,vini,
		      nmethod,neq,nzs,nactdof,itiefac,islavsurf,islavnode,
		      imastnode,nslavnode,nmastnode,ad,&au,b,&irow,icol,jq,
		      imastop,iponoels,inoels,&nzsc2,&auc2,adc2,&irowc2,jqc2,
		      islavact,gap,slavnor,slavtan,bhat,&irowbd,jqbd,&aubd,
		      &irowbdtil,jqbdtil,&aubdtil,&irowbdtil2,jqbdtil2,
		      &aubdtil2,&irowdd,jqdd,&audd,&irowddtil,jqddtil,&auddtil,
		      &irowddtil2,jqddtil2,&auddtil2,&irowddinv,jqddinv,
		      &auddinv,irowt,jqt,aut,irowtinv,jqtinv,
		      autinv,mi,ipe,ime,tietol,cstress,cstressini,
		      bp,nk,nboun,ndirboun,nodeboun,xbounact,nmpc,
		      ipompc,nodempc,coefmpc,ikboun,ilboun,ikmpc,ilmpc,
		      nslavspc,islavspc,nslavmpc,islavmpc,
		      nmastmpc,imastmpc,
		      pslavdual,islavactdof,islavtie,
		      plicon,nplicon,npmat_,nelcon,&dtime,islavnodeinv,&Bd,
		      &irowb,jqb,&Bdhelp,&irowbhelp,jqbhelp,&Dd,&irowd,jqd,
		      &Ddtil,&irowdtil,jqdtil,&Bdtil,&irowbtil,jqbtil,
		      &bet,cfsini,
		      &reltime,ithermal,plkcon,nplkcon);
	  
	nzs[0]=nzs[1];
	nzs[2]=nzs[1];
	symmetryflag=2;
	inputformat=3; 
      }

      /* storing the residuum in resold (for line search) */


      if((((*mortar==1)&&(iit!=1)&&(*ne-ne0>0)&&(*nmethod!=4))||
          ((damage_linesearch_mode==1)&&(damage_de12_enabled)&&
           (idamagereeq==0)&&(ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&
           (*ithermal<2)&&(*idrct==0)))){
	isiz=neq[1];cpypardou(resold,b,&isiz,&num_cpus);
      }
	  
      newstep=0;
      
      if(*nmethod==0){
	  
	/* error occurred in mafill: storing the geometry in frd format */
	  
	*nmethod=0;
	++*kode;
	NNEW(inum,ITG,*nk);ITGMEMSET(inum,0,*nk,1);
	if(strcmp1(&filab[1044],"ZZS")==0){
	  NNEW(neigh,ITG,40**ne);
	  MNEW(ipneigh,ITG,*nk);
	}
	  
	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,
	    ielprop,prop,sti,damn,&errn);

	if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);} 
#ifdef COMPANY
	FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
	SFREE(inum);
	if(nmethodold==0){FORTRAN(stopwithout201,());}else{FORTRAN(stop,());}
	  
      }
      
      /* implicit step (static or dynamic) */
      
      if(*iexpl<=1){
	if((*nmethod==4)&&(*mortar<2)){
	      
	  /* mechanical part */
	      
	  if(*ithermal!=2){
	    scal1=bet*dtime*dtime*(1.+alpha[0]);
	    for(k=0;k<neq[0];++k){
	      ad[k]=adb[k]+scal1*ad[k];
	    }
	    for(k=0;k<nzs[0];++k){
	      au[k]=aub[k]+scal1*au[k];
	    }
		  
	    /* upper triangle of asymmetric matrix */
		  
	    if(nasym>0){
	      for(k=nzs[2];k<nzs[2]+nzs[0];++k){
		au[k]=aub[k]+scal1*au[k];
	      }
	    }

	    /* damping */
		  
	    if(idamping==1){
	      scal1=gam*dtime*(1.+alpha[0]);
	      for(k=0;k<neq[0];++k){
		ad[k]+=scal1*adc[k];
	      }
	      for(k=0;k<nzs[0];++k){
		au[k]+=scal1*auc[k];
	      }
		      
	      /* upper triangle of asymmetric matrix */
		      
	      if(nasym>0){
		for(k=nzs[2];k<nzs[2]+nzs[0];++k){
		  au[k]+=scal1*auc[k];
		}
	      }
	    }

	  }
	      
	  /* thermal part */
	      
	  if(*ithermal>1){
	    for(k=neq[0];k<neq[1];++k){
	      ad[k]=adb[k]/dtime+ad[k];
	    }
	    for(k=nzs[0];k<nzs[1];++k){
	      au[k]=aub[k]/dtime+au[k];
	    }
		  
	    /* upper triangle of asymmetric matrix */
		  
	    if(nasym>0){
	      for(k=nzs[2]+nzs[0];k<nzs[2]+nzs[1];++k){
		au[k]=aub[k]/dtime+au[k];
	      }
	    }
	  }
	}
      
	/*	for(k=0;k<neq[1];++k){printf("fext=%" ITGFORMAT ",%f\n",k,fext[k]);}
	for(k=0;k<neq[1];++k){printf("f=%" ITGFORMAT ",%f\n",k,f[k]);}
	for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	for(k=0;k<neq[1];++k){printf("ad=%" ITGFORMAT ",%f\n",k,ad[k]);}
	for(k=0;k<nzs[1];++k){printf("au=%" ITGFORMAT ",%f\n",k,au[k]);}
	for(k=0;k<nzs[1];++k){printf("irow=%" ITGFORMAT ",%d\n",k,irow[k]);}
	for(k=0;k<neq[1]+1;++k){printf("jq=%" ITGFORMAT ",%d\n",k,jq[k]);}
	for(k=0;k<neq[1];++k){printf("icol=%" ITGFORMAT ",%d %d\n",k,icol[k],jq[k+1]-jq[k]);}*/
      

	/* Coupled dissipation control, part 1: capture the load direction.
	
	   f_hat = dR/dlambda needs no separate results() call.  At the first
	   iteration of an increment the prescribed values have just jumped
	   from lambda_n to lambda_n+dtheta and no correction has been applied
	   yet, so the residual standing in b is exactly the out-of-balance
	   that jump produced,
	
	       R_1 = -K_fp*u_hat*dtheta = f_hat*dtheta
	
	   and likewise P_1-P_n = k_pp*dtheta because du is still zero.  Both
	   fall out of the first iteration for free.
	
	   The staggered version failed here for a reason worth recording: it
	   could only move lambda once the residual was contracting, but at a
	   turning point the residual does not contract until lambda moves.
	   Solving both together removes that deadlock, because the constraint
	   becomes part of the linear system instead of a layer on top. */
	
	if((damage_diss_ctrl==2)&&(damage_diss_fhat!=NULL)){
	  if((iit==1)&&(dtheta>1.e-30)){
	    for(k=0;k<neq[1];k++) damage_diss_fhat[k]=b[k]/dtheta;
	    damage_diss_kpp=(damage_diss_p-damage_diss_pprev)/dtheta;
	    damage_diss_have=1;
	  }
	  if(damage_diss_have==1){
	    for(k=0;k<neq[1];k++) damage_diss_uf[k]=damage_diss_fhat[k];
	  }
	  if(damage_diss_probe==1){
	    printf("[DISS-GATE] it=%" ITGFORMAT " have=%" ITGFORMAT
	           " eng=%" ITGFORMAT " isolver=%" ITGFORMAT
	           " ithermal=%" ITGFORMAT " dtheta=%.4e\n",
	           iit,damage_diss_have,damage_diss_engaged,*isolver,
	           *ithermal,dtheta);
	    fflush(stdout);
	  }
	}



	/* ---- CCX_PATHFOLLOW_SOLVECHECK: keep the right-hand side -------
	   The method issues a SECOND solve with the same operator.  That is
	   only legitimate if the solver leaves ad/au intact and returns a
	   true solution both times.  Saving the rhs here lets both solves be
	   verified afterwards with one sparse mat-vec each. */

	if((pf_on==1)&&(getenv("CCX_PATHFOLLOW_SOLVECHECK")!=NULL)){
	  if(pf_rhs0==NULL){NNEW(pf_rhs0,double,neq[1]);NNEW(pf_y,double,neq[1]);}
	  isiz=neq[1];cpypardou(pf_rhs0,b,&isiz,&num_cpus);
	}

	/* ---- CCX_PATHFOLLOW: capture f_hat -----------------------------
	   b still holds fext-f = -R.  At the first iteration of an attempt
	   the ONLY thing that moved since the committed state is the
	   prescribed pattern, by pf_dlamjump, so -R = f_hat*pf_dlamjump and
	   the reference load vector conjugate to lambda comes out of the
	   iteration for free.  It is then frozen for the rest of the
	   increment, which is what makes the constraint exactly
	   differentiable. */

	/* ---- TOPOLOGY / RANK / RESIDUAL REPORT ------------------------
	   Runs at the first iteration of every attempt from CCX_TOPODIAG on.
	   b holds fext-f = -R and ad/au are assembled, so this is the real
	   operator and the real residual of the state the solver is about
	   to work on.  pardiso.c copies ad/au into its own array before
	   factorising, so the extra solve below cannot disturb the
	   factorisation the run then performs, and the scratch vector keeps
	   b untouched. */

	if((td_armed!=0)&&(td_from>0)&&(iinc>=td_from)&&(iit==1)&&
	   (*ithermal<2)){
	  double *td_x=NULL,*td_r0=NULL,*td_q=NULL,td_nx=0.,td_nb,td_pr;
	  ITG td_it,td_k,td_m,td_nm;
#define TD_NMODE 3

	  topodiag_run(&td_rep,td_comp,kon,ipkon,lakon,*ne,*nk,nactdof,mt,
	               nodeboun,ndirboun,*nboun,ipompc,nodempc,*nmpc,
	               ad,au,jq,irow,neq[1],nzs[0],b);
	  printf("[TOPODIAG] --- inc=%" ITGFORMAT " icutb=%" ITGFORMAT
	         " iter=%" ITGFORMAT " committed_state=%016llx ---\n",
	         iinc,icutb,iit,td_state);
	  topodiag_print(&td_rep,"assembled operator",iinc);

	  /* Softest mode by inverse iteration on the SAME operator, then
	     the residual's projection onto it.  If the residual lies along
	     a soft mode the obstruction is the topology; if it does not,
	     the topology is sound and the corrector is the problem.  That
	     is the whole point of the measurement. */

	  /* THREE soft directions, deflated against each other, not one.
	     A single inverse iteration converges to one vector, and a
	     residual can be orthogonal to that one while lying squarely in
	     a two- or three-dimensional soft SUBSPACE - which is exactly
	     what six rigid-body modes of a detached piece would look like.
	     The span projection below is the quantity that cannot be
	     fooled that way, and the random control says what "small"
	     means in 29501 dimensions. */

	  NNEW(td_x,double,neq[1]);
	  NNEW(td_r0,double,neq[1]);
	  NNEW(td_q,double,TD_NMODE*neq[1]);
	  isiz=neq[1];cpypardou(td_r0,b,&isiz,&num_cpus);
	  td_nm=0;
	  for(td_m=0;td_m<TD_NMODE;td_m++){
	    for(td_k=0;td_k<neq[1];td_k++)
	      td_x[td_k]=(double)(((td_k+7919*td_m)*2654435761U)%20011)
	                 /10005.-1.;
	    if(td_nm>0){
	      if(topodiag_deflate(td_x,td_q,td_nm,neq[1])==0) break;
	    }
	    for(td_it=0;td_it<3;td_it++){
	      td_nb=0.;
	      for(td_k=0;td_k<neq[1];td_k++) td_nb+=td_x[td_k]*td_x[td_k];
	      td_nb=sqrt(td_nb);
	      if(td_nb>0.) for(td_k=0;td_k<neq[1];td_k++) td_x[td_k]/=td_nb;
	      if(*isolver==0){
#ifdef SPOOLES
	        spooles(ad,au,adb,aub,&sigma,td_x,icol,irow,&neq[0],&nzs[0],
	                &symmetryflag,&inputformat,&nzs[2]);
#endif
	      }else if(*isolver==7){
#ifdef PARDISO
	        pardiso_main(ad,au,adb,aub,&sigma,td_x,icol,irow,&neq[0],
	                     &nzs[0],&symmetryflag,&inputformat,jq,&nzs[2],
	                     &nrhs);
#endif
	      }else{td_m=TD_NMODE;break;}
	      td_nx=0.;
	      for(td_k=0;td_k<neq[1];td_k++) td_nx+=td_x[td_k]*td_x[td_k];
	      td_nx=sqrt(td_nx);
	      if(td_nm>0) topodiag_deflate(td_x,td_q,td_nm,neq[1]);
	    }
	    if(td_m>=TD_NMODE) break;
	    if(topodiag_deflate(td_x,td_q,td_nm,neq[1])==0) break;
	    for(td_k=0;td_k<neq[1];td_k++) td_q[td_nm*neq[1]+td_k]=td_x[td_k];
	    td_nm++;
	    {
	      ITG td_bn=-1,td_bd=0,td_bc=-2,td_e,tb,tf;
	      double td_bv=0.;
	      for(td_k=0;td_k<*nk;td_k++){
	        for(td_it=1;td_it<mt;td_it++){
	          td_e=nactdof[mt*td_k+td_it];
	          if(td_e>0){
	            if(fabs(td_x[td_e-1])>td_bv){
	              td_bv=fabs(td_x[td_e-1]);td_bn=td_k+1;td_bd=td_it;
	              td_bc=td_comp[td_k];
	            }
	          }
	        }
	      }
	      topodiag_support(td_bn,kon,ipkon,lakon,*ne,&tb,&tf);
	      printf("[TOPODIAG]   soft mode %" ITGFORMAT
	             ": 1/sigma_min >= %.6e, peaks at node %" ITGFORMAT
	             " dir %" ITGFORMAT " (component %" ITGFORMAT
	             ", %" ITGFORMAT " live bulk element(s), %" ITGFORMAT
	             " live facet(s)), |cos(R,mode)|=%.3e\n",
	             td_nm,td_nx,td_bn,td_bd,td_bc,tb,tf,
	             topodiag_project(td_r0,td_x,neq[1]));
	    }
	  }
	  /* Keep the basis for the rest of the attempt: the decisive
	     question is not where the RESIDUAL points but where the
	     CORRECTION does.  With 1/sigma_min at 6e13 a residual component
	     along a soft mode of 1e-7 - which the cosine prints as zero -
	     becomes a correction of order one, and that is what a runaway
	     Newton step looks like. */

	  if((td_qkeep==NULL)||(td_qneq!=neq[1])){
	    if(td_qkeep!=NULL) SFREE(td_qkeep);
	    NNEW(td_qkeep,double,TD_NMODE*neq[1]);
	    td_qneq=neq[1];
	  }
	  isiz=TD_NMODE*neq[1];cpypardou(td_qkeep,td_q,&isiz,&num_cpus);
	  td_nmode=td_nm;

	  td_pr=topodiag_project_span(td_r0,td_q,td_nm,neq[1]);
	  for(td_k=0;td_k<neq[1];td_k++)
	    td_x[td_k]=(double)(((td_k+104729)*2246822519U)%20011)/10005.-1.;
	  printf("[TOPODIAG]   share of |R| in the span of the %" ITGFORMAT
	         " softest modes: %.6e   (random control %.6e)\n",
	         td_nm,td_pr,topodiag_project(td_r0,td_x,neq[1]));
	  {
	    ITG tb,tf;
	    topodiag_support(td_rep.resmaxnode,kon,ipkon,lakon,*ne,&tb,&tf);
	    printf("[TOPODIAG]   residual peak node %" ITGFORMAT
	           " has %" ITGFORMAT " live bulk element(s), %" ITGFORMAT
	           " live facet(s)\n",td_rep.resmaxnode,tb,tf);
	  }
	  SFREE(td_q);SFREE(td_r0);SFREE(td_x);
	  fflush(stdout);
	}

	if((pf_on==1)&&(iit==1)){

	  /* HOW FAR HAS THE FROZEN REFERENCE VECTOR DRIFTED?
	     f_hat is frozen at the first capture, and the corrector adds
	     dlambda*K^-1 f_hat to the displacement while adding dlambda to
	     lambda.  Those two are consistent only while f_hat is -dR/dlambda;
	     a scale error s makes the applied load-factor step 1/s of the one
	     the displacement correction assumed.  b/dlamjump IS the current
	     secant of -dR/dlambda, so the comparison costs one pass over the
	     vector and needs no extra residual evaluation.  Measured, not
	     assumed: on the mixed benchmark the ratio is still 0.997 at
	     increment 1500. */

	  if((pathfollow_have()==1)&&(fabs(pf_dlamjump)>1.e-8)){
	    const double *pffh=pathfollow_fhat();
	    double pfa=0.,pfb=0.,pfd=0.,pfs;
	    for(k=0;k<neq[1];k++){
	      pfs=b[k]/pf_dlamjump;
	      pfa+=pffh[k]*pffh[k];pfb+=pfs*pfs;pfd+=pffh[k]*pfs;
	    }
	    if((pfa>0.)&&(pfb>0.)){
	      pf_fhcos=pfd/sqrt(pfa*pfb);
	      pf_fhrat=sqrt(pfa/pfb);
	    }
	  }
	  pathfollow_capture(b,pf_dlamjump);
	  if(pathfollow_have()==1){
	    pathfollow_setPn(pf_project(pathfollow_fhat(),vini,pf_uref,
	                                nactdof,*nk,mt));
	  }
	}

	/* Stabilisation of detached pieces.

	   A piece that has come loose by FACE connectivity still touches the
	   structure at a node or an edge, so it keeps free rotational modes
	   and the operator is singular in those directions.  At the
	   m12_field wall three independent pointers - the softest mode
	   (1/sigma_min = 230 against 3.7 elastic), the dof whose Newton
	   correction equalled the whole increment, and the largest residual -
	   all sat on nodes carrying a one-element floating piece (E-56 and
	   its follow-up).

	   Deleting those pieces was implemented and measured: it turns R0,
	   R1, R2 and the UC6 PASS disk from PASS to FAIL (E-57).  So hold
	   them instead.  Only nodes whose every surviving element is
	   unreached are touched, which leaves every load-carrying dof alone,
	   and the added stiffness is a fraction of the mean diagonal so the
	   piece is held without being welded back on.

	   DEFAULT OFF.  CCX_DAMAGE_STABILISE=<alpha> switches it on. */

	if((damage_stab_alpha>0.)&&(damage_de12_enabled)&&(*ithermal<2)){
	  ITG nstabnode=0;
	  NNEW(damage_stab_node,ITG,*nk);
	  FORTRAN(damfloatstab,(ipkon,kon,lakon,&ne0,nk,nodeboun,nboun,
				ipompc,nodempc,nmpc,damage_stab_node,
				&nstabnode));
	  /* Second criterion, and the one that matters.

	     `damfloatstab` flags nodes unreachable by FACE connectivity.  It
	     does NOT flag a node whose single supporting element is still
	     face-attached to the body but is DEAD - and that is the measured
	     case.  On m14_fine, node 776 had exactly one live element, that
	     element stood at D = 1.0000, the node travelled 9.63 mm while the
	     other three nodes of the same element moved 0.44, and the
	     element's longest edge went from 0.0974 to 9.809 - a stretch of
	     101x on a 4 mm specimen.

	     A node whose entire live support sits at g = 1-D near gmin is very
	     nearly free, whatever the connectivity says.  Flag it on the
	     degradation, which is what E-24 recorded as missing from every
	     topology predicate in this branch.  The threshold is deliberately
	     severe: the whole support must be below 1% of its stiffness. */
	  {
	    double *gsup=NULL;
	    const double *dsrc=(damage_visc_eta>0.&&damage_damvisc!=NULL)?
	      damage_damvisc:dam;
	    ITG nn,jj,ndead=0;
	    NNEW(gsup,double,*nk);
	    for(nn=0;nn<*nk;nn++) gsup[nn]=-1.;
	    for(i=0;i<ne0;i++){
	      if(ipkon[i]<0) continue;
	      if(strcmp1(&lakon[8*i],"C3D4")!=0) continue;
	      {
		double dmx=0.,g;
		ITG usedam=((damage_visc_eta>0.&&damage_damvisc!=NULL)?0:1);
		for(jj=0;jj<mi[0];jj++){
		  /* dam holds 1+D, visc holds D - see damage_de13_mark_deadsole */
		  double dd=usedam?(dsrc[mi[0]*i+jj]-1.):dsrc[mi[0]*i+jj];
		  if(dd<0.) dd=0.;
		  if(dd>1.) dd=1.;
		  if(dd>dmx) dmx=dd;
		}
		g=1.-dmx; if(g<0.) g=0.;
		for(jj=0;jj<4;jj++){
		  nn=kon[ipkon[i]+jj]-1;
		  if((nn<0)||(nn>=*nk)) continue;
		  if(g>gsup[nn]) gsup[nn]=g;
		}
	      }
	    }
	    for(nn=0;nn<*nk;nn++){
	      if(gsup[nn]<0.) continue;              /* no live support at all */
	      if(gsup[nn]>=DAMAGE_STAB_GDEAD) continue;
	      if(damage_stab_node[nn]==0){
		damage_stab_node[nn]=1;
		nstabnode++;
		ndead++;
	      }
	    }
	    SFREE(gsup);
	    if(ndead>damage_stab_maxdead){
	      damage_stab_maxdead=ndead;
	      printf("[DAMAGE STABILISE] inc=%" ITGFORMAT " nodes held on DEAD "
		     "support only: %" ITGFORMAT " (g < %.1e) (new maximum)\n",
		     iinc,ndead,DAMAGE_STAB_GDEAD);
	      fflush(stdout);
	    }
	  }

	  if(nstabnode>0){
	    double dref=0.;
	    ITG ndref=0;
	    for(k=0;k<neq[1];k++){
	      if(ad[k]>0.){dref+=ad[k];ndref++;}
	    }
	    if(ndref>0){
	      dref/=(double)ndref;
	      ITG nstabdof=0;
	      for(i=0;i<*nk;i++){
		if(damage_stab_node[i]==0) continue;
		for(idir=1;idir<=3;idir++){
		  k=nactdof[mt*i+idir];
		  if(k>0){
		    ad[k-1]+=damage_stab_alpha*dref;
		    nstabdof++;
		  }
		}
	      }
	      if(nstabdof>damage_stab_maxdof){
		damage_stab_maxdof=nstabdof;
		printf("[DAMAGE STABILISE] inc=%" ITGFORMAT " nodes=%"
		       ITGFORMAT " dof=%" ITGFORMAT " alpha=%.3e "
		       "mean_diag=%.6e (new maximum)\n",
		       iinc,nstabnode,nstabdof,damage_stab_alpha,dref);
		fflush(stdout);
	      }
	    }
	  }
	  SFREE(damage_stab_node);
	}

        /* ---- [DAMAGE REG] positive diagonal regularization, level-3
           attempts only.  Same diagonal the stabiliser above edits, same
           moment: after the whole assembly, before any solver sees it,
           and the solver frees ad right after the factorisation so
           nothing persists.  Only the DIRECTION changes - b holds the
           residual, and checkconvergence decides on ram/cam/qa/uam from
           results()/calcresidual, none of which sees the shift. */

        if((damage_reg_on==1)&&(*ithermal<2)){
          ITG rneg=0,rzero=0,nsum=0,nfl=0,nneg2=0;
          double rmin=1.e300,rmax=-1.e300,rabs,dsum=0.,dmean,dfloor;
          double ssum=0.,smax=0.,sh;
          for(k=0;k<neq[1];k++){
            if(ad[k]<0.) rneg++;
            if(ad[k]==0.) rzero++;
            if(ad[k]<rmin) rmin=ad[k];
            if(ad[k]>rmax) rmax=ad[k];
            rabs=(ad[k]<0.)?-ad[k]:ad[k];
            if(rabs>0.){dsum+=rabs;nsum++;}
          }
          dmean=(nsum>0)?dsum/nsum:0.;
          /* Cheap guard.  NaN fails every comparison, so !(dmean>0.)
             catches NaN, zero and negative alike; the second test
             catches infinity.  On any of them the shift is not applied
             and the attempt proceeds exactly as the stock one. */
          if((!(dmean>0.))||(dmean>1.e300)){
            printf("[DAMAGE REG] mean|ad| is not a usable scale "
                   "(zero, NaN or infinite); regularization NOT applied, "
                   "the attempt falls back to stock%s","\n");
            fflush(stdout);
            damage_reg_on=0;
          }else{
            dfloor=1.e-6*dmean;
            /* D must be STRICTLY positive.  Plain |ad| is not: it
               vanishes where the diagonal is zero, and for ad<0 the
               shift gives |ad|*(lambda-1), so lambda=1 lands exactly ON
               zero.  Hence the floor and a ladder that runs past 1. */
            for(k=0;k<neq[1];k++){
              rabs=(ad[k]<0.)?-ad[k]:ad[k];
              if(rabs<dfloor){rabs=dfloor;nfl++;}
              sh=damage_reg_lambda*rabs;
              ssum+=sh; if(sh>smax) smax=sh;
              ad[k]+=sh;
            }
            for(k=0;k<neq[1];k++) if(ad[k]<0.) nneg2++;
            damage_reg_napply++;
            if(damage_reg_napply==1){
              printf("[DAMAGE REG] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                     " diagonal sign census BEFORE the shift - a "
                     "DIAGNOSTIC of the diagonal, NOT evidence about the "
                     "definiteness of K: n=%" ITGFORMAT " negative=%"
                     ITGFORMAT " zero=%" ITGFORMAT " min=%.6e max=%.6e "
                     "mean|ad|=%.6e%s",
                     iinc,iit,neq[1],rneg,rzero,rmin,rmax,dmean,"\n");
              printf("[DAMAGE REG] applied ad[k] += lambda*D[k] with "
                     "D[k]=max(|ad[k]|,%.6e), lambda=%.3e (ladder step %"
                     ITGFORMAT " of %" ITGFORMAT ").  ACTUAL SHIFT: "
                     "mean=%.6e max=%.6e floored_dof=%" ITGFORMAT
                     " ; diagonal negatives %" ITGFORMAT " -> %" ITGFORMAT
                     ".  Acceptance remains on the UNMODIFIED residual%s",
                     dfloor,damage_reg_lambda,damage_reg_level+1,
                     damage_reg_nlam,ssum/neq[1],smax,nfl,rneg,nneg2,"\n");
              fflush(stdout);
            }
          }
        }

        /* ---- [DAMAGE TR] the residual, before the solver overwrites b ----
           b holds fext-f here (calcresidual.c:48-52), i.e. exactly -R with
           R=f_int-f_ext.  phi=1/2|R|^2=1/2|b|^2 is sign-blind; the GRADIENT
           is not, so the sign is carried explicitly where d=J^T b is formed.
           Nothing between the calcresidual that filled b and this point
           writes b on the mortar<=1 path. */
        damage_dl_lc_due=0;
        if((damage_dl_lincheck>0)&&(iinc==damage_dl_lincheck)&&
           (iit>=damage_dl_lc_it)&&(iit<damage_dl_lc_it+damage_dl_lc_nit)&&
           (damage_dl_on==0)&&(damage_de12_enabled)&&(idamagereeq==0)&&
           (ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&
           (*idrct==0)&&(*mortar<=1)) damage_dl_lc_due=1;
        if((damage_dl_on==1)||(damage_dl_lc_due==1)){
          if(damage_dl_r0==NULL) NNEW(damage_dl_r0,double,neq[1]);
          isiz=neq[1];cpypardou(damage_dl_r0,b,&isiz,&num_cpus);
        }
        /* [DAMAGE CT] b holds fext-f = -R here; the solve overwrites it.
           The FULL read-before-write set is snapshotted HERE, at the current
           iterate, BEFORE the solve and before the stock results() applies
           z.  Taking it after that call mixed u_current with the history of
           u_current+z, so every finite difference was built on two
           different base states.  cam is RESTORED from this snapshot, never
           re-initialised: at this point it already holds the clean values
           the iteration top installed. */
        if(damage_ct_on==1){
          ITG cq;
          if(damage_ct_r0==NULL){
            NNEW(damage_ct_r0,double,neq[1]);
            NNEW(damage_ct_beps,double,neq[1]);
            NNEW(damage_ct_y,double,neq[1]);
            NNEW(damage_ct_z,double,neq[1]);
          }
          if(damage_ct_dam==NULL){
            NNEW(damage_ct_dam,double,mi[0]**ne);
            NNEW(damage_ct_visc,double,mi[0]**ne);
            if(*nstate_>0) NNEW(damage_ct_xs,double,*nstate_*mi[0]**ne);
            NNEW(damage_ct_jac,double,12*mi[0]**ne);
            NNEW(damage_ct_sgn,ITG,mi[0]*ne0);
          }
          isiz=neq[1];cpypardou(damage_ct_r0,b,&isiz,&num_cpus);
          isiz=mi[0]**ne;cpypardou(damage_ct_dam,dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_ct_visc,damage_damvisc,&isiz,&num_cpus);
          }
          if(*nstate_>0){
            isiz=*nstate_*mi[0]**ne;
            cpypardou(damage_ct_xs,xstate,&isiz,&num_cpus);
          }
          if(damage_damjac!=NULL){
            isiz=12*mi[0]**ne;
            cpypardou(damage_ct_jac,damage_damjac,&isiz,&num_cpus);
          }
          for(cq=0;cq<4;cq++) damage_ct_qa[cq]=qa[cq];
          for(cq=0;cq<5;cq++) damage_ct_cam[cq]=cam[cq];
          for(cq=0;cq<2;cq++) damage_ct_uam[cq]=uam[cq];
          damage_ct_lamsnap=damage_ct_lam;
          damage_ct_nfact++;
        }

	if(*isolver==0){
#ifdef SPOOLES
	  if(*ithermal<2){
	    spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
		    &symmetryflag,&inputformat,&nzs[2]);
	    
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    spooles(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
		    &sigma,&b[neq[0]],&icol[neq[0]],iruc,
		    &n1,&n2,&symmetryflag,&inputformat,&nzs[2]);
	  }else{
	    spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
		    &symmetryflag,&inputformat,&nzs[2]);
	  }
#else
	  printf(" *ERROR in nonlingeo: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if((*isolver==2)||(*isolver==3)){
	  if(symmetryflag==2){
	    if(*isolver==3){
	      printf(" *WARNING in nonlingeo: the iterative Cholesky solver");
	      printf(" cannot be used for asymmetric matrices.\nThe");
	      printf(" iterative scaling solver will be used instead\n\n");
	    }
	    NNEW(rwork,double,neq[1]);
	    NNEW(sol,double,neq[1]);
	    RENEW(au,double,2*nzs[1]+neq[1]);
	    isiz=neq[1];cpypardou(&au[2*nzs[1]],ad,&isiz,&num_cpus);
	    nelt=2*nzs[1]+neq[1];
	    lrgw=131+16*neq[1];
	    isym=0;
	    NNEW(rgwk,double,lrgw);
	    NNEW(igwk,ITG,20);
	    for(i=0;i<neq[1];i++){
	      rwork[i]=1./ad[i];}
	    FORTRAN(predgmres_struct,(&neq[1],b,sol,&nelt,irow,jq,au,
				      &isym,&itol,&tol,&itmax,&iter,
				      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
				      &ligw,rwork,iwork));
	    isiz=neq[1];cpypardou(b,sol,&isiz,&num_cpus);
	    SFREE(rgwk);SFREE(igwk);SFREE(rwork);SFREE(sol);
	  }else{
	    preiter(ad,&au,b,&icol,&irow,&neq[1],&nzs[1],isolver,iperturb);
	  }
	}
	else if(*isolver==4){
#ifdef SGI
	  if(symmetryflag==2){
	    printf(" *ERROR in nonlingeo: the SGI solver cannot be used for asymmetric matrices\n\n");
	    FORTRAN(stop,());
	  }
	  token=1;
	  if(*ithermal<2){
	    sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],token);
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    sgi_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
		     &sigma,&b[neq[0]],&icol[neq[0]],iruc,
		     &n1,&n2,token);
	  }else{
	    sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],token);
	  }
#else
	  printf(" *ERROR in nonlingeo: the SGI library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	  if(symmetryflag==2){
	    printf(" *ERROR in nonlingeo: the TAUCS solver cannot be used for asymmetric matrices\n\n");
	    FORTRAN(stop,());
	  }
	  tau(ad,&au,adb,aub,&sigma,b,icol,&irow,&neq[1],&nzs[1]);
#else
	  printf(" *ERROR in nonlingeo: the TAUCS library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==6){
#ifdef MATRIXSTORAGE
	  matrixstorage(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1],
			ntrans,inotr,trab,co,nk,nactdof,jobnamec,mi,ipkon,
			lakon,kon,ne,mei,nboun,nmpc,cs,mcs,ithermal,nmethod);
	  strcpy2(fneig,jobnamec,132);
	  strcat(fneig,".frd");
	  if((f1=fopen(fneig,"ab"))==NULL){
	    printf(" *ERROR in nonlingeo: cannot open frd file for writing...");
	    exit(0);
	  }
	  fprintf(f1," 9999\n");
	  fclose(f1);
	  FORTRAN(stopwithout201,());
#else
	  printf(" *ERROR in nonlingeo: the MATRIXSTORAGE library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	  if(*ithermal<2){
	    pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    NNEW(jqtherm,ITG,n1+1);
	    for(i=0;i<n1+1;i++){
	      jqtherm[i]=jq[neq[0]+i]-nzs[0];}
	    pardiso_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			 &sigma,&b[neq[0]],&icol[neq[0]],iruc,
			 &n1,&n2,&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	    SFREE(jqtherm);
	  }else{
	    pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }
#else
	  printf(" *ERROR in nonlingeo: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	else if(*isolver==8){
#ifdef PASTIX
	  if(*ithermal<2){
	    pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[0],&nzs[0],
			&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }else if((*ithermal==2)&&(uncoupled)){
	    n1=neq[1]-neq[0];
	    n2=nzs[1]-nzs[0];
	    NNEW(jqtherm,ITG,n1+1);
	    for(i=0;i<n1+1;i++){
	      jqtherm[i]=jq[neq[0]+i]-nzs[0];}
	    pastix_main(&ad[neq[0]],&au[nzs[0]],&adb[neq[0]],&aub[nzs[0]],
			&sigma,&b[neq[0]],&icol[neq[0]],iruc,
			&n1,&n2,&symmetryflag,&inputformat,jqtherm,&nzs[2],&nrhs);
	    SFREE(jqtherm);
	  }else{
	    pastix_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq[1],&nzs[1],
			&symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
	  }
#else
	  printf(" *ERROR in nonlingeo: the PASTIX library is not linked\n\n");
	  FORTRAN(stop,());
#endif
	}
	  

	/* Coupled dissipation control, part 2: second solve and the load
	   factor increment.
	
	       K du_R = R      (already done, du_R is in b)
	       K du_F = f_hat  (same matrix, second back substitution)
	       du     = du_R + dlambda*du_F
	
	   The constraint row for displacement control, with
	   dG = (P_n*lambda - P*lambda_n)/2 and g = dG - dtau, gives
	
	       dlambda = -[g + l_n/2*f_hat.du_R]
	                 /[l_n/2*f_hat.du_F + P_n/2 - l_n/2*k_pp]
	
	   PARDISO is called a second time on the same ad/au, so the
	   factorisation is repeated; correctness first, and the symbolic
	   analysis is already cached when CCX_PARDISO_REUSE_SYMBOLIC is set. */
	
	/* The constraint is unsatisfiable while the response is still
	   elastic: dG is identically zero there, so g = -dtau every
	   iteration and lambda is driven at the clamp forever.  That is
	   what the first coupled run did - residual 0.088 against an
	   average force of 0.023, diverging at increment 2.  Stay in
	   displacement control until the measured dissipation of an
	   accepted increment reaches a fraction of the target, then
	   engage and stay engaged. */

	if((damage_diss_ctrl==2)&&(damage_diss_have==1)&&
	   (damage_diss_engaged==1)&&
	   (*isolver==7)&&(*ithermal<2)){
#ifdef PARDISO
	  pardiso_main(ad,au,adb,aub,&sigma,damage_diss_uf,icol,irow,
	               &neq[0],&nzs[0],&symmetryflag,&inputformat,jq,
	               &nzs[2],&nrhs);
#endif
	  damage_diss_fr=0.;
	  damage_diss_ff=0.;
	  for(k=0;k<neq[1];k++){
	    damage_diss_fr+=damage_diss_fhat[k]*b[k];
	    damage_diss_ff+=damage_diss_fhat[k]*damage_diss_uf[k];
	  }
	  damage_diss_dgcur=0.5*(damage_diss_pprev*damage_diss_lamcur-
	                         damage_diss_p*damage_diss_lprev);
	  damage_diss_g=damage_diss_dgcur-damage_diss_target;
	  /* P is a function of u alone, so dg/dlambda is exactly 0.5*P_n and
	     every remaining lambda dependence already travels through du_F.
	     The earlier -0.5*lambda_n*k_pp counted dP/dlambda a second time,
	     which is what left the correction cancelling itself: the residual
	     sat at 0.043639 while the displacement correction was 5.7e-7. */
	  damage_diss_den=0.5*damage_diss_lprev*damage_diss_ff
	                 +0.5*damage_diss_pprev;
	  if(damage_diss_probe==1){
	    printf("[DISS-PROBE] it=%" ITGFORMAT " lam=%.6f dG=%.4e g=%.4e"
	           " fr=%.4e ff=%.4e den=%.4e\n",
	           iit,damage_diss_lamcur,damage_diss_dgcur,damage_diss_g,
	           damage_diss_fr,damage_diss_ff,damage_diss_den);
	    fflush(stdout);
	  }
	  if(fabs(damage_diss_den)>1.e-30){
	    damage_diss_dlam=-(damage_diss_g
	                       +0.5*damage_diss_lprev*damage_diss_fr)
	                     /damage_diss_den;
	    if(damage_diss_dlam>DAMAGE_DISS_DLAM*dthetaref)
	      damage_diss_dlam=DAMAGE_DISS_DLAM*dthetaref;
	    if(damage_diss_dlam<-DAMAGE_DISS_DLAM*dthetaref)
	      damage_diss_dlam=-DAMAGE_DISS_DLAM*dthetaref;
	    for(k=0;k<neq[1];k++)
	      b[k]+=damage_diss_dlam*damage_diss_uf[k];
	    damage_diss_lamcur+=damage_diss_dlam;
	    if(damage_diss_probe==1){
	      printf("[DISS-PROBE]      dlam=%.6e -> lam=%.6f\n",
	             damage_diss_dlam,damage_diss_lamcur);
	      fflush(stdout);
	    }
	    if(damage_arc==1){
	      if(damage_diss_lamcur<0.) damage_diss_lamcur=0.;
	    }else if(damage_diss_lamcur<=theta){
	      damage_diss_lamcur=theta+(*tmin);
	    }
	    if(damage_diss_lamcur>1.) damage_diss_lamcur=1.;
	    for(k=0;k<*nboun;k++){
	      xbounact[k]=xbounold[k]+
	        (xboun[k]-xbounold[k])*damage_diss_lamcur;
	    }
	  }
	}

	/* ---- CCX_PATHFOLLOW: close the bordered system -----------------
	   b now holds du_R.  Solve once more with the SAME operator for
	   du_F = K^-1 f_hat and close the second row.  The extra
	   factorisation is the honest price of staying solver agnostic: the
	   CalculiX solver wrappers factor and solve in one call, so a second
	   right-hand side cannot reuse the first factorisation without
	   changing them.  Correctness first; the cost is a constant factor
	   and is reported at the end of the step. */

	/* ---- where does the CORRECTION point? --------------------------
	   b now holds the ordinary Newton correction du_R = K^-1(-R), before
	   any constraint, line search or trust region touches it.  If its
	   norm is dominated by the soft subspace measured at iteration 1,
	   the runaway steps are the ill-conditioning being amplified, and no
	   amount of step-size control fixes that.  If it is not, the steps
	   are large for a constitutive reason. */

	if((td_armed!=0)&&(td_from>0)&&(iinc>=td_from)&&(td_qkeep!=NULL)&&
	   (td_qneq==neq[1])&&(td_nmode>0)&&(*ithermal<2)){
	  double tdn=0.,tds;
	  ITG tdk;
	  for(tdk=0;tdk<neq[1];tdk++) tdn+=b[tdk]*b[tdk];
	  tdn=sqrt(tdn);
	  tds=topodiag_project_span(b,td_qkeep,td_nmode,neq[1]);
	  printf("[TOPODIAG]   correction inc=%" ITGFORMAT " iter=%" ITGFORMAT
	         " |du|2=%.6e  share in the %" ITGFORMAT
	         " softest modes: %.6e\n",iinc,iit,tdn,td_nmode,tds);
	  fflush(stdout);
	}

	/* Reset here, not inside the block below: the line search has to be
	   able to tell "this iteration applied a constraint step" from "the
	   value left over by an earlier one". */

	pf_applied=0;pf_dlamit=0.;pf_lam0it=pf_lam;

	if((pf_on==1)&&(pathfollow_have()==1)&&(*ithermal<2)){
	  const double *pf_fh=pathfollow_fhat();
	  for(k=0;k<neq[1];k++) pf_uf[k]=pf_fh[k];
	  if(*isolver==0){
#ifdef SPOOLES
	    spooles(ad,au,adb,aub,&sigma,pf_uf,icol,irow,&neq[0],&nzs[0],
		    &symmetryflag,&inputformat,&nzs[2]);
#endif
	  }else if(*isolver==7){
#ifdef PARDISO
	    pardiso_main(ad,au,adb,aub,&sigma,pf_uf,icol,irow,&neq[0],&nzs[0],
			 &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#endif
	  }
	  pf_nstep++;

	  /* ff is recorded on every iterate, engaged or not, so that the
	     tangent predictor has a stiffness the moment the constraint
	     takes over.  Only the APPLICATION is gated on engagement. */

	  if((pf_rhs0!=NULL)&&(getenv("CCX_PATHFOLLOW_SOLVECHECK")!=NULL)){
	    double pfr1=0.,pfr2=0.,pfn1=0.,pfn2=0.,pft;
	    ITG pfi,pfone=1;
	    const double *pffh2=pathfollow_fhat();
	    for(pfi=0;pfi<neq[1];pfi++) pf_y[pfi]=0.;
	    FORTRAN(op,(b,pf_y,ad,au,jq,irow,&pfone,&neq[0]));
	    for(pfi=0;pfi<neq[0];pfi++){
	      pft=pf_y[pfi]-pf_rhs0[pfi];
	      pfr1+=pft*pft;pfn1+=pf_rhs0[pfi]*pf_rhs0[pfi];
	    }
	    for(pfi=0;pfi<neq[1];pfi++) pf_y[pfi]=0.;
	    FORTRAN(op,(pf_uf,pf_y,ad,au,jq,irow,&pfone,&neq[0]));
	    for(pfi=0;pfi<neq[0];pfi++){
	      pft=pf_y[pfi]-pffh2[pfi];
	      pfr2+=pft*pft;pfn2+=pffh2[pfi]*pffh2[pfi];
	    }
	    printf("[PF-SOLVE] it=%" ITGFORMAT " |K*duR-rhs|/|rhs|=%.3e "
	           "|K*duF-fhat|/|fhat|=%.3e\n",iit,
	           sqrt(pfr1)/(sqrt(pfn1)+1.e-300),
	           sqrt(pfr2)/(sqrt(pfn2)+1.e-300));
	    fflush(stdout);
	  }

/* Evaluate CalculiX's own residual at (vold, LAMV) and copy fext-f into
   DST.  b is zeroed first so results() leaves the displacement where it
   is: the probe must read the residual AT a state, never move it.

   iout=-1, not 0.  resultsini gates the correction on iout>-1 and the
   prescribed-dof update on abs(iout)<2, so -1 sets the prescribed dofs
   from xbounact (which IS "evaluating R at this lambda") while leaving the
   free dofs alone - and, crucially, it does not request the output fields.
   iout=0 does request them, and results() then writes inum/een/emn/epn,
   which are not allocated on this path: measured as a SIGSEGV inside
   memset called from nonlingeo. */

#define PF_LIN_RESID(LAMV,DST) do{                                        \
    ITG _k,_io=iout;                                                      \
    /* Every evaluation starts from the identical base state.  Without    \
       this the probe measures its own leakage: idempotency (evaluating   \
       R at lambda0 twice) came out at 1.7e+05 relative near the limit    \
       point, which is exactly the spurious "discontinuity" the sweep     \
       then reported.  vold is NOT restored here - the caller perturbs it \
       on purpose. */                                                     \
    isiz=*nstate_*mi[0]**ne;cpypardou(xstate,pf_sxs,&isiz,&num_cpus);     \
    isiz=27*mi[0]**ne;cpypardou(xstiff,pf_sxst,&isiz,&num_cpus);          \
    isiz=neq[1];cpypardou(f,pf_sf,&isiz,&num_cpus);                       \
    isiz=mt**nk;cpypardou(fn,pf_sfn,&isiz,&num_cpus);                     \
    isiz=6*mi[0]**ne;cpypardou(stx,pf_sstx,&isiz,&num_cpus);              \
    if((dam!=NULL)&&(pf_sdam!=NULL)){                                     \
      isiz=mi[0]**ne;cpypardou(dam,pf_sdam,&isiz,&num_cpus);}             \
    for(_k=0;_k<4;_k++) qa[_k]=pf_sqa[_k];                                \
    for(_k=0;_k<5;_k++) cam[_k]=pf_scam[_k];                              \
    for(_k=0;_k<*nboun;_k++)                                              \
      xbounact[_k]=xbounold[_k]+(xboun[_k]-xbounold[_k])*(LAMV);          \
    for(_k=0;_k<neq[1];_k++) b[_k]=0.;                                    \
    iout=-1;                                                              \
    isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);                        \
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,                      \
            elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,         \
            ielorien,norien,orab,ntmat_,t0,t1act,ithermal,                \
            prestr,iprestr,filab,eme,emn,een,iperturb,                    \
            f,fn,nactdof,&iout,qa,vold,b,nodeboun,                        \
            ndirboun,xbounact,nboun,ipompc,                               \
            nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,        \
            accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,           \
            plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,matname,    \
            mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,     \
            ener,enern,emeini,xstaten,eei,enerini,cocon,ncocon,set,       \
            nset,istartset,iendset,ialset,nprint,prlab,prset,qfx,qfn,     \
            trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,ilmpc,istep,     \
            &iinc,springarea,&reltime,&ne0,thicke,shcon,nshcon,           \
            sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,            \
            pmastsurf,mortar,islavact,cdn,islavnode,nslavnode,ntie,       \
            clearini,islavsurf,ielprop,prop,energyini,energy,&kscale,     \
            iponoeln,inoeln,nener,orname,network,ipobody,xbodyact,        \
            ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,     \
            t0g,t1g,islavquadel,aut,irowt,jqt,&mortartrafoflag,           \
            &intscheme,physcon,dam,damn,iponoel);                         \
    calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux2,vold,vini,       \
                 &dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,      \
                 fini,islavnode,nslavnode,mortar,ntie,mi,nzs,&nasym,      \
                 &idamping,veold,adc,auc,cvini,cv,&alpham,&num_cpus);     \
    isiz=neq[1];cpypardou((DST),b,&isiz,&num_cpus);                       \
    iout=_io;                                                             \
  }while(0)

	  /* ================= CCX_PATHFOLLOW_LINCHECK =====================
	     The decisive gate.  A synthetic self test proves the formulas of
	     pathfollow.c; it proves nothing about how u, lambda, xbounact, b
	     and K are actually wired together inside CalculiX.  This measures
	     that wiring, at the REAL base point - after prediction() has
	     extrapolated and after the operator has been assembled - by
	     finite differences of the residual CalculiX itself computes.

	     Checked, from one base point, restoring the full trial state
	     between evaluations, for a sweep of decreasing eps:

	       q_FD  = [R(u,lambda+eps) - R(u,lambda)]/eps
	       row 1 : [R(u+eps*p, lambda+eps*dl) - R(u,lambda)]/eps
	               against  K*p + q_FD*dl

	     and the frozen reference vector f_hat against -q_FD, which is
	     the claim "a constant reference load reproduces the actual
	     dR/dlambda" that must not be assumed. */

	  if((pf_lincheck>0)&&(iinc==pf_lincheck)&&(iit==1)){
	    ITG pq,pfone=1,pj,pnode,pi;
	    double peps,plam0,pdl,pn1,pn2,pnq,pt,pfhq,pnf,pkp,pql,pqc;
	    double pepsv[5]={1.e-3,1.e-4,1.e-5,1.e-6,1.e-7};

	    if(pf_p==NULL){
	      NNEW(pf_p,double,neq[1]);NNEW(pf_r0,double,neq[1]);
	      NNEW(pf_r1,double,neq[1]);NNEW(pf_q,double,neq[1]);
	      NNEW(pf_sv,double,mt**nk);
	      NNEW(pf_sxs,double,*nstate_*mi[0]**ne);
	      if(dam!=NULL) NNEW(pf_sdam,double,mi[0]**ne);
	      NNEW(pf_sf,double,neq[1]);NNEW(pf_sfn,double,mt**nk);
	      NNEW(pf_sstx,double,6*mi[0]**ne);
	      NNEW(pf_sxb,double,*nboun);
	      /* pf_y is the mat-vec scratch; it is otherwise allocated only
	         by SOLVECHECK, and LINCHECK dereferenced it as NULL. */
	      if(pf_y==NULL) NNEW(pf_y,double,neq[1]);
	      NNEW(pf_sxst,double,27*mi[0]**ne);
	      NNEW(pf_lhs,double,neq[1]);NNEW(pf_rhs,double,neq[1]);
	    }

	    /* --- snapshot everything the probe is about to disturb --- */
	    isiz=mt**nk;cpypardou(pf_sv,vold,&isiz,&num_cpus);
	    isiz=*nstate_*mi[0]**ne;cpypardou(pf_sxs,xstate,&isiz,&num_cpus);
	    if((dam!=NULL)&&(pf_sdam!=NULL)){
	      isiz=mi[0]**ne;cpypardou(pf_sdam,dam,&isiz,&num_cpus);}
	    isiz=neq[1];cpypardou(pf_sf,f,&isiz,&num_cpus);
	    isiz=mt**nk;cpypardou(pf_sfn,fn,&isiz,&num_cpus);
	    isiz=6*mi[0]**ne;cpypardou(pf_sstx,stx,&isiz,&num_cpus);
	    isiz=*nboun;cpypardou(pf_sxb,xbounact,&isiz,&num_cpus);
	    isiz=27*mi[0]**ne;cpypardou(pf_sxst,xstiff,&isiz,&num_cpus);
	    for(k=0;k<4;k++) pf_sqa[k]=qa[k];
	    for(k=0;k<5;k++) pf_scam[k]=cam[k];
	    isiz=neq[1];cpypardou(pf_p,b,&isiz,&num_cpus);   /* direction p */
	    plam0=pf_lam;
	    pdl=(fabs(pf_dlampred)>0.)?pf_dlampred:1.e-3;

	    /* p is left as du_R, the Newton correction itself, and is NOT
	       normalised.  Rescaling it to unit norm makes eps*p a
	       perturbation of order 1e-3 in displacement, which is larger
	       than the cohesive d0 = 4e-4 and therefore no longer probes the
	       branch the solver is on; the sweep then stopped converging for
	       a reason that has nothing to do with the tangent.  du_R is the
	       direction the method actually takes, which is the one worth
	       verifying. */
	    pt=sqrt(damage_fd_dot(pf_p,pf_p,neq[0]));

	    printf("\n[PF-LIN] increment %" ITGFORMAT ", base point AFTER "
	           "prediction(); lambda=%.10f  |p|=%.4e  dlambda_dir=%.4e\n",
	           iinc,plam0,sqrt(damage_fd_dot(pf_p,pf_p,neq[1])),pdl);

	    /* R at the base point */
	    PF_LIN_RESID(plam0,pf_r0);
	    pn1=sqrt(damage_fd_dot(pf_r0,pf_r0,neq[0]));
	    printf("[PF-LIN]   |R0| = %.6e   neq0=%" ITGFORMAT " neq1=%"
	           ITGFORMAT "\n",pn1,neq[0],neq[1]);

	    /* bisection: is ONE evaluation self-contaminating, or is it the
	       perturbed ones that leave something behind? */
	    PF_LIN_RESID(plam0,pf_r1);
	    pn2=0.;
	    for(k=0;k<neq[0];k++){pt=pf_r1[k]-pf_r0[k];pn2+=pt*pt;}
	    printf("[PF-LIN]   immediate re-evaluation at lambda0: "
	           "|dR|/|R0| = %.6e\n",(pn1>0.)?sqrt(pn2)/pn1:sqrt(pn2));
	    fflush(stdout);

	    printf("[PF-LIN]   eps        |q_FD|      "
	           "cos(f_hat,-q_FD)  |f_hat|/|q_FD|   row1 rel.err\n");
	    for(pq=0;pq<10;pq++){
	      if(pq==5){
	        /* Second sweep with dlambda = 0.  This removes the load term
	           entirely and tests ONLY whether the assembled operator is
	           the derivative of the residual CalculiX computes, which is
	           a statement about the element tangent and has nothing to do
	           with path following. */
	        pdl=0.;
	        printf("[PF-LIN]   ---- dlambda = 0: pure tangent test, "
	               "[R(u+eps*p)-R(u)]/eps  against  K*p ----\n");
	      }
	      peps=pepsv[pq%5];

	      /* Re-anchor the base residual for every eps.  A single
	         evaluation is idempotent (measured: |dR|/|R0| = 0 exactly),
	         but the perturbed ones must not be allowed to bias the next
	         difference. */
	      PF_LIN_RESID(plam0,pf_r0);

	      /* ORDER MATTERS.  pf_q was originally computed first and read
	         last, and it did not survive the calls in between: |q| taken
	         at the norm loop disagreed with the |q_FD| printed from the
	         same buffer moments earlier (5.3e+03 against 8.1e-11).  The
	         perturbed evaluation and the mat-vec are therefore done
	         FIRST, and q is rebuilt immediately before it is used, so
	         nothing can run between its definition and its use. */

	      /* (a) full perturbation, in u and lambda together */
	      for(pi=0;pi<*nk;pi++){
	        for(pj=1;pj<mt;pj++){
	          pnode=nactdof[mt*pi+pj];
	          if(pnode>0) vold[mt*pi+pj]+=peps*pf_p[pnode-1];
	        }
	      }
	      PF_LIN_RESID(plam0+peps*pdl,pf_r1);
	      isiz=mt**nk;cpypardou(vold,pf_sv,&isiz,&num_cpus);
	      for(k=0;k<neq[0];k++) pf_lhs[k]=(pf_r1[k]-pf_r0[k])/peps;

	      /* (b) K*p */
	      for(k=0;k<neq[1];k++) pf_y[k]=0.;
	      FORTRAN(op,(pf_p,pf_y,ad,au,jq,irow,&pfone,&neq[0]));

	      /* (c) q_FD, built last */
	      PF_LIN_RESID(plam0+peps,pf_r1);
	      for(k=0;k<neq[0];k++) pf_q[k]=(pf_r1[k]-pf_r0[k])/peps;
	      pnq=sqrt(damage_fd_dot(pf_q,pf_q,neq[0]));

	      pfhq=0.;pnf=0.;
	      if(pathfollow_have()==1){
	        const double *pfh=pathfollow_fhat();
	        pfhq=damage_fd_dot(pfh,pf_q,neq[0]);
	        pnf=sqrt(damage_fd_dot(pfh,pfh,neq[0]));
	      }

	      pn2=0.;pn1=0.;pkp=0.;pql=0.;pqc=0.;
	      for(k=0;k<neq[0];k++) pf_rhs[k]=-pf_y[k]+pdl*pf_q[k];
	      for(k=0;k<neq[0];k++){
	        pt=pf_lhs[k]-pf_rhs[k];  pn2+=pt*pt;
	        pn1+=pf_lhs[k]*pf_lhs[k];
	        pkp+=pf_y[k]*pf_y[k];
	        pql+=pf_rhs[k]*pf_rhs[k];
	        pqc+=pf_q[k]*pf_q[k];
	      }
	      printf("[PF-LIN]   %.1e  %.6e  %+.10f   %12.6e   %.6e"
	             "   |lhs|=%.4e |Kp|=%.4e |rhs|=%.4e  |q|inloop=%.6e\n",
	             peps,pnq,(pnf*pnq>0.)?pfhq/(pnf*pnq):0.,
	             (pnq>0.)?pnf/pnq:0.,
	             (pn1>0.)?sqrt(pn2)/sqrt(pn1):sqrt(pn2),
	             sqrt(pn1),sqrt(pkp),sqrt(pql),sqrt(pqc));
	      fflush(stdout);
	    }

	    /* Idempotency of the probe itself.  Everything above is only
	       meaningful if evaluating R at the SAME point twice returns the
	       same vector; if it does not, the sweep is measuring the probe's
	       own leakage rather than dR/dlambda. */

	    PF_LIN_RESID(plam0,pf_r1);
	    pn1=0.;pn2=0.;
	    for(k=0;k<neq[0];k++){
	      pt=pf_r1[k]-pf_r0[k];pn2+=pt*pt;pn1+=pf_r0[k]*pf_r0[k];
	    }
	    printf("[PF-LIN]   idempotency: |R(lam0) again - R0|/|R0| = %.6e"
	           "   %s\n",(pn1>0.)?sqrt(pn2)/sqrt(pn1):sqrt(pn2),
	           (sqrt(pn2)<=1.e-10*sqrt(pn1))?"clean":"PROBE LEAKS STATE");

	    /* --- restore --- */
	    isiz=mt**nk;cpypardou(vold,pf_sv,&isiz,&num_cpus);
	    isiz=*nstate_*mi[0]**ne;cpypardou(xstate,pf_sxs,&isiz,&num_cpus);
	    if((dam!=NULL)&&(pf_sdam!=NULL)){
	      isiz=mi[0]**ne;cpypardou(dam,pf_sdam,&isiz,&num_cpus);}
	    isiz=neq[1];cpypardou(f,pf_sf,&isiz,&num_cpus);
	    isiz=mt**nk;cpypardou(fn,pf_sfn,&isiz,&num_cpus);
	    isiz=6*mi[0]**ne;cpypardou(stx,pf_sstx,&isiz,&num_cpus);
	    isiz=*nboun;cpypardou(xbounact,pf_sxb,&isiz,&num_cpus);
	    isiz=27*mi[0]**ne;cpypardou(xstiff,pf_sxst,&isiz,&num_cpus);
	    for(k=0;k<4;k++) qa[k]=pf_sqa[k];
	    for(k=0;k<5;k++) cam[k]=pf_scam[k];
	    isiz=neq[1];cpypardou(b,pf_p,&isiz,&num_cpus);
	    pf_lam=plam0;
	    printf("[PF-LIN] state restored\n\n");
	    fflush(stdout);
	  }
#undef PF_LIN_RESID

	  pathfollow_measure(pf_uf);
	  pf_applied=0;
	  if(pf_engaged==1){

	    /* the trial displacement is READ FROM THE MODEL, not
	       accumulated; see the design note in pathfollow.c */

	    if(pf_codmode>=1){

	      /* Mode 1 (CCX_PATHFOLLOW_COD): phi is ABSOLUTE, the mean normal
	         separation measured from the reference configuration.
	         Mode 2 (CCX_CRACK_CONTROL): phi is INCREMENTAL, c^T(u-u_n),
	         because the functional itself is refrozen every attempt and
	         only its increment is meaningful across a refreeze.

	         c only couples plus/minus node pairs, so in both cases the
	         projection is exactly the controlled opening. */

	      pf_lam0it=pf_lam;
	      pf_dlamit=0.;
	      pf_cu=pf_project(pathfollow_cod_c(),vold,
	                       (pf_codmode==2)?vini:pf_uref,nactdof,*nk,mt);
	      pf_applied=pathfollow_cod_step(b,pf_uf,pf_cu,&pf_lam,pf_clip,
	                                     &pf_g,&pf_dlam,&pf_reason);
	      pf_dg=pf_cu;
	      if(pf_applied==1){
	        pf_dlamit=pf_lam-pf_lam0it;
	        for(k=0;k<*nboun;k++){
	          xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*pf_lam;
	        }
	      }
	      goto pf_after_step;
	    }
	    pf_pdu=pf_project(pathfollow_fhat(),vold,vini,nactdof,*nk,mt);

	    /* At the first iteration put lambda exactly on the constraint,
	       in closed form, so that the coupled iteration does not open
	       with a residual the stock extrapolation put there. */

	    /* Measured: this does NOT rescue the descending branch (1 accepted
	       increment past engagement against 2 without it), so it is
	       opt-in and off by default.  It is kept because it isolates one
	       hypothesis cleanly - the opening constraint residual is not
	       what stops the run. */

	    if((iit==1)&&(getenv("CCX_PATHFOLLOW_PROJECT")!=NULL)){
	      if(pathfollow_project_lambda(pf_pdu,&pf_lam)==1){
	        for(k=0;k<*nboun;k++){
	          xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*pf_lam;
	        }
	      }
	    }
	    pf_applied=pathfollow_step(b,pf_uf,pf_pdu,&pf_lam,pf_clip,
	                               &pf_dg,&pf_g,&pf_dlam,&pf_reason);
	  }
	  if(pf_applied==1){
	    for(k=0;k<*nboun;k++){
	      xbounact[k]=xbounold[k]+(xboun[k]-xbounold[k])*pf_lam;
	    }
	  }
	 pf_after_step:
	  if(getenv("CCX_PATHFOLLOW_PROBE")!=NULL){
	    printf("[PATHFOLLOW] it=%" ITGFORMAT " lambda=%.8f dlambda=%+.4e "
	           "dG=%.6e g=%.4e applied=%" ITGFORMAT " reason=%" ITGFORMAT
	           "\n",iit,pf_lam,pf_applied?pf_dlam:0.,pf_dg,pf_g,
	           pf_applied,pf_reason);
	    fflush(stdout);
	  }
	}



	/* Locate the softest mode of the assembled operator.

	   At the stall an unloading of 0.5% moves the structure by 2e-2 in a
	   corner far from the crack, and |du|/|r| has grown three orders from
	   its mid-run value.  That is a soft mode, and because it survives a
	   change of load it belongs to the topology, not to lambda - which is
	   why walking lambda downwards does not rescue the run.

	   Inverse iteration finds it: x <- K^-1 x, renormalised, converges to
	   the smallest-singular-value direction, and ||x||/||b|| after the
	   first solve estimates 1/sigma_min.  Reporting where that vector
	   lives says which part of the model is the mechanism, instead of
	   guessing at candidates one at a time.

	   The factorisation is already in PARDISO's hands, so the extra
	   solves are back-substitutions. */

	/* [WALLDIAG] Each deflated inverse iteration is a full factorisation
	   and solve.  Armed by hand on ONE named increment that is affordable;
	   armed by the load-factor gate it would run on every iteration of
	   every increment from the wall on, so under that gate it is taken
	   once per attempt.  CCX_DAMAGE_NULLVEC on its own is unchanged. */
	if((damage_null_inc>0)&&(iinc>=damage_null_inc)&&
	   ((damage_wall_armed==0)||(iit==1))&&
	   (*isolver==7)&&(*ithermal<2)&&(*mortar<=1)){
	  NNEW(damage_null_x,double,neq[1]);
	  for(k=0;k<neq[1];k++){
	    damage_null_seed=(1103515245*damage_null_seed+12345)&0x7fffffff;
	    damage_null_x[k]=(double)(damage_null_seed%20001-10000)/10000.;
	  }
	  for(damage_null_it=0;damage_null_it<damage_null_nit;damage_null_it++){
	    damage_null_nb=0.;
	    for(k=0;k<neq[1];k++)
	      damage_null_nb+=damage_null_x[k]*damage_null_x[k];
	    damage_null_nb=sqrt(damage_null_nb);
	    if(damage_null_nb>0.)
	      for(k=0;k<neq[1];k++) damage_null_x[k]/=damage_null_nb;
#ifdef PARDISO
	    pardiso_main(ad,au,adb,aub,&sigma,damage_null_x,icol,irow,
			 &neq[0],&nzs[0],&symmetryflag,&inputformat,jq,
			 &nzs[2],&nrhs);
#endif
	    damage_null_nx=0.;
	    for(k=0;k<neq[1];k++)
	      damage_null_nx+=damage_null_x[k]*damage_null_x[k];
	    damage_null_nx=sqrt(damage_null_nx);
	    if((damage_null_it<3)||(damage_null_it==damage_null_nit-1))
	    printf("[DAMAGE NULLVEC] inc=%" ITGFORMAT " iteration %" ITGFORMAT
		   " 1/sigma_min >= %.6e\n",iinc,damage_null_it+1,
		   damage_null_nx);
	  }
	  /* where does the mode live */
	  damage_null_amax=0.;
	  for(i=0;i<*nk;i++){
	    for(idir=1;idir<=3;idir++){
	      k=nactdof[mt*i+idir];
	      if(k<=0) continue;
	      if(fabs(damage_null_x[k-1])>damage_null_amax)
		damage_null_amax=fabs(damage_null_x[k-1]);
	    }
	  }
	  printf("[DAMAGE NULLVEC] nodes carrying the mode "
		 "(|component| > 0.2 of the maximum):\n");
	  damage_null_cnt=0;
	  for(i=0;i<*nk;i++){
	    damage_null_nn=0.;
	    for(idir=1;idir<=3;idir++){
	      k=nactdof[mt*i+idir];
	      if(k<=0) continue;
	      if(fabs(damage_null_x[k-1])>damage_null_nn)
		damage_null_nn=fabs(damage_null_x[k-1]);
	    }
	    if((damage_null_amax>0.)&&(damage_null_nn>0.2*damage_null_amax)){
	      damage_null_cnt++;
	      if(damage_null_cnt<=20){
		printf("[DAMAGE NULLVEC]   node %" ITGFORMAT
		       " at (%.3f, %.3f, %.3f) amplitude %.4f\n",
		       i+1,co[3*i],co[3*i+1],co[3*i+2],
		       damage_null_nn/damage_null_amax);
	      }
	    }
	  }
	  printf("[DAMAGE NULLVEC] %" ITGFORMAT " node(s) above the "
		 "threshold out of %" ITGFORMAT "\n",damage_null_cnt,*nk);
	  fflush(stdout);
	  SFREE(damage_null_x);
	  FORTRAN(stopwithout201,());
	}

	if(*mortar<=1){
	  if(isensitivity){
	    SFREE(adcpy);MNEW(adcpy,double,neq[1]);
	    SFREE(aucpy);MNEW(aucpy,double,(nasym+1)*nzs[1]);
	    isiz=neq[1];cpypardou(adcpy,ad,&isiz,&num_cpus);
	    isiz=(nasym+1)*nzs[1];cpypardou(aucpy,au,&isiz,&num_cpus);
	  }
	  /* ---- [DAMAGE TR] everything the dogleg needs from J, taken while
	     J is still allocated.  Two lines below, ad and au are freed, and
	     every trial residual evaluation happens after that - so a dogleg
	     that wanted a matrix-vector product per trial could not have one.
	     It does not need one:

	         J p_N = r0        (exactly: p_N is what the solver returned)
	         J d   = w         (computed here, once)

	     and every step this trust region can propose is p = pa*d + pb*p_N,
	     so J p = pa*w + pb*r0 is a two-scalar combination of vectors that
	     already exist.  No matrix, no second factorisation and no second
	     solve for the whole trial loop.

	     Storage (add_sm_st_as.f:29-56, mastruct.c:795-820, cross-checked
	     against opas.f:34-46): jq and irow are 1-based; column c (0-based)
	     owns au[jq[c]-1 .. jq[c+1]-2]; for slot k in that column with
	     r=irow[k]-1 > c,  au[k] = J(r,c)  and  au[nzs[2]+k] = J(c,r);
	     ad[i] = J(i,i).  sigma is 0 everywhere in nonlingeo, so the
	     factorised operator is (ad,au) with no shift.

	     THE TRANSPOSE IS PROVED, NOT ASSERTED.  d = J^T r0 and J p_N = r0
	     give dot(d,p_N) = r0^T J J^-1 r0 = |r0|^2 identically.  Swap the
	     two halves of au, or use J where J^T is meant, and the identity
	     fails at once.  It is checked on every armed iteration and the
	     mechanism REFUSES TO ARM when it is off by more than 1e-3. */
	  if(((damage_dl_on==1)||(damage_dl_lc_due==1))&&
	     (damage_dl_r0!=NULL)&&(*mortar<=1)){
	    if((nasym!=1)||(symmetryflag!=2)||(*ithermal>=2)||
	       (nzs[2]!=nzs[1])||(neq[0]!=neq[1])){
	      if(damage_dl_have>=0){
	        printf("[DAMAGE TR] REFUSING TO ARM: the assembled operator is "
	               "not the one this construction was proved on (nasym=%"
	               ITGFORMAT " symmetryflag=%" ITGFORMAT " ithermal=%"
	               ITGFORMAT " neq0=%" ITGFORMAT " neq1=%" ITGFORMAT
	               " nzs1=%" ITGFORMAT " nzs2=%" ITGFORMAT ").  The wall "
	               "goes to the original stock stop.%s",
	               nasym,symmetryflag,*ithermal,neq[0],neq[1],
	               nzs[1],nzs[2],"\n");
	        fflush(stdout);
	      }
	      damage_dl_have=(damage_dl_on==1)?-1:damage_dl_have;
	      damage_dl_on=0;damage_dl_lc_due=0;
	    }else{
	      ITG dk,dc,dr;
	      double daden,dad;
	      if(damage_dl_d==NULL){
	        NNEW(damage_dl_d,double,neq[1]);
	        NNEW(damage_dl_w,double,neq[1]);
	        NNEW(damage_dl_pn,double,neq[1]);
	      }
	      /* d = J^T r0 */
	      for(dk=0;dk<neq[1];dk++)
	        damage_dl_d[dk]=ad[dk]*damage_dl_r0[dk];
	      for(dc=0;dc<neq[1];dc++){
	        for(dk=jq[dc]-1;dk<jq[dc+1]-1;dk++){
	          dr=irow[dk]-1;
	          damage_dl_d[dc]+=au[dk]*damage_dl_r0[dr];
	          damage_dl_d[dr]+=au[nzs[2]+dk]*damage_dl_r0[dc];
	        }
	      }
	      /* w = J d   (the same loop with the two halves exchanged) */
	      for(dk=0;dk<neq[1];dk++)
	        damage_dl_w[dk]=ad[dk]*damage_dl_d[dk];
	      for(dc=0;dc<neq[1];dc++){
	        for(dk=jq[dc]-1;dk<jq[dc+1]-1;dk++){
	          dr=irow[dk]-1;
	          damage_dl_w[dr]+=au[dk]*damage_dl_d[dc];
	          damage_dl_w[dc]+=au[nzs[2]+dk]*damage_dl_d[dr];
	        }
	      }
	      isiz=neq[1];cpypardou(damage_dl_pn,b,&isiz,&num_cpus);
	      /* [WALLDIAG] MASKSTEP.  pm = p_N with every component that
	         belongs to an AUTOSPC-masked node zeroed, and wm = J pm by the
	         same proved loop.  Built HERE because ad and au are freed
	         before the linearisation check runs; pm and wm are the only
	         things that survive to it.  Nothing on a solution path reads
	         either. */
	      if((damage_wall_maskstep!=0)&&(damage_spc_mask!=NULL)&&
	         (damage_spc_nk>=*nk)){
	        ITG mi_,mj_,mk_;
	        if(damage_dl_pm==NULL){
	          NNEW(damage_dl_pm,double,neq[1]);
	          NNEW(damage_dl_wm,double,neq[1]);
	        }
	        for(dk=0;dk<neq[1];dk++) damage_dl_pm[dk]=damage_dl_pn[dk];
	        for(mi_=0;mi_<*nk;mi_++){
	          if(damage_spc_mask[mi_]==0) continue;
	          for(mj_=1;mj_<mt;mj_++){
	            mk_=nactdof[mt*mi_+mj_];
	            if(mk_>0) damage_dl_pm[mk_-1]=0.;
	          }
	        }
	        for(dk=0;dk<neq[1];dk++)
	          damage_dl_wm[dk]=ad[dk]*damage_dl_pm[dk];
	        for(dc=0;dc<neq[1];dc++){
	          for(dk=jq[dc]-1;dk<jq[dc+1]-1;dk++){
	            dr=irow[dk]-1;
	            damage_dl_wm[dr]+=au[dk]*damage_dl_pm[dc];
	            damage_dl_wm[dc]+=au[nzs[2]+dk]*damage_dl_pm[dr];
	          }
	        }
	        damage_dl_npm2=0.;
	        for(dk=0;dk<neq[1];dk++)
	          damage_dl_npm2+=damage_dl_pm[dk]*damage_dl_pm[dk];
	      }
	      damage_dl_nb2=0.;damage_dl_nd2=0.;damage_dl_nw2=0.;
	      damage_dl_npn2=0.;damage_dl_dtpn=0.;
	      for(dk=0;dk<neq[1];dk++){
	        damage_dl_nb2+=damage_dl_r0[dk]*damage_dl_r0[dk];
	        damage_dl_nd2+=damage_dl_d[dk]*damage_dl_d[dk];
	        damage_dl_nw2+=damage_dl_w[dk]*damage_dl_w[dk];
	        damage_dl_npn2+=damage_dl_pn[dk]*damage_dl_pn[dk];
	        damage_dl_dtpn+=damage_dl_d[dk]*damage_dl_pn[dk];
	      }
	      /* how far from symmetric the operator actually IS - measured,
	         because the dogleg is only worth building if J^T != J */
	      damage_dl_asym=0.;daden=0.;
	      for(dk=0;dk<nzs[1];dk++){
	        dad=fabs(au[dk]-au[nzs[2]+dk]);
	        if(dad>damage_dl_asym) damage_dl_asym=dad;
	        if(fabs(au[dk])>daden) daden=fabs(au[dk]);
	      }
	      if(daden>0.) damage_dl_asym/=daden;
	      damage_dl_ident=(damage_dl_nb2>0.)?
	                       damage_dl_dtpn/damage_dl_nb2:0.;
	      if((damage_dl_nb2<=0.)||(damage_dl_nw2<=0.)||
	         (damage_dl_nd2<=0.)||(fabs(damage_dl_ident-1.)>1.e-3)){
	        printf("[DAMAGE TR] TRANSPOSE-CHECK FAILED at inc=%" ITGFORMAT
	               " iter=%" ITGFORMAT ": dot(J^T R,p_N)/|R|^2 = %.12e, "
	               "and it must be 1.  Either the transpose is not a "
	               "transpose or the solve is not a solve.  The dogleg "
	               "DISARMS and the wall goes to the original stock stop.%s",
	               iinc,iit,damage_dl_ident,"\n");
	        fflush(stdout);
	        damage_dl_have=-1;damage_dl_on=0;damage_dl_lc_due=0;
	      }else{
	        damage_dl_tc=damage_dl_nd2/damage_dl_nw2;
	        damage_dl_have=1;
	        if(damage_dl_on==1) damage_dl_nfact++;
	        if(damage_dl_banner==0){
	          damage_dl_banner=1;
	          printf("[DAMAGE TR] TRANSPOSE-CHECK inc=%" ITGFORMAT " iter=%"
	                 ITGFORMAT " dot(J^T R,p_N)/|R|^2 = %.12e (exact value "
	                 "1; a J used in place of J^T does not satisfy it).  "
	                 "OPERATOR ASYMMETRY max|au_L-au_U|/max|au_L| = %.6e, "
	                 "so J^T R is genuinely not J R here.%s",
	                 iinc,iit,damage_dl_ident,damage_dl_asym,"\n");
	          fflush(stdout);
	        }
	      }
	    }
	  }
	  SFREE(ad);SFREE(au);
	} 
      }
      
      /* explicit dynamic step */
      
      else{
	if(((mscalmethod==0)||(mscalmethod==2))&&(*mortar!=-1)){
	  //    for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	  if(*ithermal!=2){
	    isiz=neq[0];divparll(b,adb,&isiz,&num_cpus);
	  }
	  //      for(k=0;k<neq[1];++k){printf("b=%" ITGFORMAT ",%f\n",k,b[k]);}
	  if(*ithermal>1){
	    for(k=neq[0];k<neq[1];++k){
	      b[k]=b[k]*dtime/adb[k];
	    }
	  }
	}
	else{
	  if(*ithermal!=2){
	    if(*isolver==0){
#ifdef SPOOLES
	      spooles_solve(b,&neq[0]);
#endif
	    }
	    else if(*isolver==4){
#ifdef SGI
	      sgi_solve(b,token);
#endif
	    }
	    else if(*isolver==5){
#ifdef TAUCS
	      tau_solve(b,&neq[0]);
#endif
	    }
	    else if(*isolver==7){
#ifdef PARDISO
	      pardiso_solve(b,&neq[0],&symmetryflag,&inputformat,&nrhs);
#endif
	    }
	    else if(*isolver==8){
#ifdef PASTIX
	      pastix_solve(b,&neq[0],&symmetryflag,&nrhs);
#endif
	    }
	    if(*mortar==-1){
	      if(ncont!=0){
		if(iinc==1){
		  for(i=0;i<neqtot;i++){
		    k=floor(ltot[i]/10);
		    l=ltot[i]-10*k;
		    b[ktot[i]-1]=veold[mt*(k-1)+l];
		  }
		} else{

		  /* determine the velocity in the contact nodes */
	      
		  for(i=0;i<neqtot;++i){
		    b[ktot[i]-1]=(qb[i]-volddof[ktot[i]-1])/(dtime);
		  }
		}
		SFREE(qb);
	      }
	      SFREE(volddof);
	    }
	  }
	  if(*ithermal>1){
	    for(k=neq[0];k<neq[1];++k){
	      b[k]=b[k]*dtime/adb[k];
	    }
	  }
	}
      }
      
      /* mortar */

      if(*mortar>1){	    
  
	/* restoring the structure of the original stiffness
	   matrix */

	for(i=0;i<3;i++){
	  nzs[i]=nzstemp[i];}
	for (i=0;i<neq[1];i++){jq[i]=jqtemp[i];icol[i]=icoltemp[i];}
	jq[neq[1]]=jqtemp[neq[1]];
	for (i=0;i<nzs[1];i++){irow[i]=irowtemp[i];}
	SFREE(jqtemp);SFREE(irowtemp);SFREE(icoltemp);

	/* trafo util->u , calculate cstress and update active set  */

	stressmortar(bhat,adc2,auc2,jqc2,irowc2,neq,gap,b,islavact,irowddinv,
		     jqddinv,auddinv,irowt,jqt,aut,irowtinv,
		     jqtinv,autinv,ntie,nslavnode,islavnode,nmastnode,
		     imastnode,slavnor,slavtan,nactdof,&iflagact,cstress,
		     cstressini,mi,cdisp,f_cs,f_cm,&iit,&iinc,vold,vini,bp,
		     nk,nboun,ndirboun,nodeboun,xboun,
		     nmpc,
		     ipompc,nodempc,coefmpc,
		     nslavmpc,islavmpc,
		     tieset,elcon,tietol,ncmat_,ntmat_,plicon,nplicon,npmat_,
		     nelcon,&dtime,cfs,cfm,islavnodeinv,Bd,irowb,jqb,Dd,
		     irowd,jqd,Ddtil,irowdtil,jqdtil,Bdtil,irowbtil,jqbtil,
		     nmethod,&bet,ithermal,
		     iperturb,labmpc,cam,veold,accold,&gam,
		     cfsini,cfstil,plkcon,nplkcon,filab,f,fn,qa,nprint,prlab,
		     xforc,nforc,iponoel);
	  
	SFREE(auc2);SFREE(adc2);SFREE(irowc2);SFREE(icolc2);SFREE(jqc2);
	SFREE(au);SFREE(ad);	  
      }

      /* calculating the displacements, stresses and forces */
      
      MNEW(v,double,mt**nk);
      isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
      
      NNEW(stx,double,6*mi[0]**ne);
      MNEW(fn,double,mt**nk);

      /* for massless explicit dynamics without energy
         calculation only the displacements have to be calculated

         this does not work for explicit dynamics without massless
         contact since in that case f and fini have to be calculated  */
      
      if((*mortar==-1)&&(masslesslinear>0)&&(*nener==0)){
	resultsini(nk,v,ithermal,filab,iperturb,f,fn,
		   nactdof,&iout,qa,vold,b,nodeboun,ndirboun,
		   xboun,nboun,ipompc,nodempc,coefmpc,labmpc,nmpc,nmethod,cam,
		   neq,veold,accold,&bet,&gam,&dtime,mi,vini,nprint,prlab,
		   &intpointvarm,&calcul_fn,&calcul_f,&calcul_qa,&calcul_cauchy,
		   &ikin,&intpointvart,typeboun,&num_cpus,mortar,nener,iponoeln,
		   network);
      }else{
	if(ne1d2d==1)NNEW(inum,ITG,*nk);
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		ndirboun,xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		&reltime,&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
		inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavquadel,aut,irowt,jqt,&mortartrafoflag,
		&intscheme,physcon,dam,damn,iponoel);
	if(ne1d2d==1)SFREE(inum);
      }

      /* implicit dynamics (Matteo Pacher) */

      if((*ne!=ne0)&&(*nmethod==4)&&(*ithermal<2)&&(*iexpl<=1)){
	FORTRAN(storecontactprop,(ne,&ne0,lakon,kon,ipkon,mi,ielmat,elcon,
				  mortar,adblump,nactdof,springarea,ncmat_,
				  ntmat_,stx,&temax));
      }

      /* Reaction conjugate to the prescribed pattern u_p = lambda*xboun.
         Taken on the main Newton path, where fn is the array just filled
         by results() and resultsini.c has set calcul_fn=1 for NLGEOM.
         Two earlier placements were wrong and both showed it plainly:
         after SFREE(fn) the sum was exactly zero, and at the output calls
         it was stale enough to give a negative dissipation increment. */

      /* [LOADPATH] capture the grip reaction HERE, where fn is the array
         results() has just filled - the same reason the dissipation sum
         below is taken here.  The census prints it; it must not read fn. */
      if(damage_lp_ready) loadpath_setreaction(&damage_lp,fn,mt);

      if(damage_diss_report==1){
        damage_diss_p=0.;
        for(i=0;i<*nboun;i++){
          if((ndirboun[i]<1)||(ndirboun[i]>mi[1])) continue;
          damage_diss_p+=fn[mt*(nodeboun[i]-1)+ndirboun[i]]*xboun[i];
        }
      }

      if((damage_fd_inc>0)&&(iinc>=damage_fd_inc)&&(iit>=damage_fd_it)){

        NNEW(damage_fd_vsav,double,mt**nk);
        NNEW(damage_fd_fp,double,mt**nk);
        NNEW(damage_fd_fm,double,mt**nk);
        memcpy(damage_fd_vsav,v,sizeof(double)*mt**nk);

        /* columns are taken at the nodes of the most damaged element: a
           tangent that is right in the elastic bulk and wrong in the
           process zone is precisely the error every earlier test would
           have missed */

        damage_fd_el=-1;damage_fd_dmax=-1.;
        for(i=0;i<ne0;i++){
          if(ipkon[i]<0) continue;
          if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
          if(dam[mi[0]*i]>damage_fd_dmax){
            damage_fd_dmax=dam[mi[0]*i];damage_fd_el=i;
          }
        }
        if(damage_fd_el<0){
          printf("[STRUCT-FD] no active C3D4 element to probe\n");
        }else{
          printf("[STRUCT-FD] inc=%" ITGFORMAT " iter=%" ITGFORMAT
      	   " element=%" ITGFORMAT " dam=%.6f h=%.3e\n",
      	   iinc,iit,damage_fd_el+1,damage_fd_dmax,damage_fd_h);
          printf("[STRUCT-FD] %-6s %-5s %-14s %-14s %-12s %s\n",
      	   "node","dir","max|K_fd|","max|K_asm-K_fd|","rel.err",
      	   "worst row");
          fflush(stdout);
        }

        damage_fd_ncol=0;
        for(damage_fd_j=0;(damage_fd_j<4)&&(damage_fd_el>=0);damage_fd_j++){
          damage_fd_node=kon[ipkon[damage_fd_el]+damage_fd_j]-1;
          if((damage_fd_node<0)||(damage_fd_node>=*nk)) continue;

          for(idir=1;idir<=3;idir++){
            damage_fd_col=nactdof[mt*damage_fd_node+idir];
            if(damage_fd_col<=0) continue;

            /* central difference on the internal force */

            for(damage_fd_s=0;damage_fd_s<2;damage_fd_s++){
      	memcpy(v,damage_fd_vsav,sizeof(double)*mt**nk);
      	v[mt*damage_fd_node+idir]+=
      	  (damage_fd_s==0)?damage_fd_h:-damage_fd_h;
      	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
      		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
      		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
      		prestr,iprestr,filab,eme,emn,een,iperturb,
      		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
      		ndirboun,xbounact,nboun,ipompc,
      		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
      		accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
      		plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,
      		matname,mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,
      		ikboun,ilboun,ener,enern,emeini,xstaten,eei,enerini,
      		cocon,ncocon,set,nset,istartset,iendset,ialset,nprint,
      		prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,
      		nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,&ne0,
      		thicke,shcon,nshcon,sideload,xloadact,xloadold,&icfd,
      		inomat,pslavsurf,pmastsurf,mortar,islavact,cdn,
      		islavnode,nslavnode,ntie,clearini,islavsurf,ielprop,
      		prop,energyini,energy,&kscale,iponoeln,inoeln,nener,
      		orname,network,ipobody,xbodyact,ibody,typeboun,itiefac,
      		tieset,smscale,&mscalmethod,nbody,t0g,t1g,islavquadel,
      		aut,irowt,jqt,&mortartrafoflag,&intscheme,physcon,dam,
      		damn,iponoel);
      	if(damage_fd_s==0)
      	  memcpy(damage_fd_fp,fn,sizeof(double)*mt**nk);
      	else
      	  memcpy(damage_fd_fm,fn,sizeof(double)*mt**nk);
            }

            /* compare the measured column against the stored one */

            damage_fd_amax=0.;damage_fd_emax=0.;damage_fd_worst=-1;
            damage_fd_nbad=0;damage_fd_wa=0.;damage_fd_wf=0.;
            for(k=0;k<*nk;k++){
      	for(damage_fd_d=1;damage_fd_d<=3;damage_fd_d++){
      	  damage_fd_row=nactdof[mt*k+damage_fd_d];
      	  if(damage_fd_row<=0) continue;
      	  damage_fd_num=(damage_fd_fp[mt*k+damage_fd_d]
      			 -damage_fd_fm[mt*k+damage_fd_d])
      	                /(2.*damage_fd_h);
      	  damage_fd_asm=damage_fd_coeff(damage_fd_ad,damage_fd_au,
                                                jq,irow,nzs,
      					damage_fd_row,damage_fd_col);
      	  if(fabs(damage_fd_num)>damage_fd_amax)
      	    damage_fd_amax=fabs(damage_fd_num);
      	  if(fabs(damage_fd_asm-damage_fd_num)>damage_fd_emax){
      	    damage_fd_emax=fabs(damage_fd_asm-damage_fd_num);
      	    damage_fd_worst=k+1;
      	    damage_fd_wa=damage_fd_asm;
      	    damage_fd_wf=damage_fd_num;
      	    damage_fd_wd=damage_fd_d;
      	  }
      	  if(fabs(damage_fd_asm-damage_fd_num)>1.e-4*damage_fd_amax)
      	    damage_fd_nbad++;
      	}
            }
            printf("[STRUCT-FD] %-6" ITGFORMAT " %-5" ITGFORMAT
      	     " %-13.5e %-12.4e worst row %" ITGFORMAT "/%" ITGFORMAT
      	     " K_fd=%-13.5e K_asm=%-13.5e ratio=%-9.4f bad=%"
      	     ITGFORMAT "\n",
      	     damage_fd_node+1,idir,damage_fd_amax,
      	     (damage_fd_amax>0.)?damage_fd_emax/damage_fd_amax:-1.,
      	     damage_fd_worst,damage_fd_wd,damage_fd_wf,damage_fd_wa,
      	     (fabs(damage_fd_wf)>0.)?damage_fd_wa/damage_fd_wf:0.,
      	     damage_fd_nbad);
            fflush(stdout);
            damage_fd_ncol++;
          }
        }

        memcpy(v,damage_fd_vsav,sizeof(double)*mt**nk);
        SFREE(damage_fd_vsav);SFREE(damage_fd_fp);SFREE(damage_fd_fm);
        printf("[STRUCT-FD] %" ITGFORMAT " columns compared; stopping, the "
      	 "probe perturbed v repeatedly and this run is diagnostic "
      	 "only\n",damage_fd_ncol);
        fflush(stdout);
        FORTRAN(stopwithout201,());
      }


      /* Dissipation-controlled load factor.

         Shrinking the step cannot pass a limit point: with the target
         tightened threefold the stall moved by 0.0006 in step time, so at
         that lambda no equilibrium exists on the branch being tracked and
         the load factor itself has to become an unknown.

         This is the staggered form of the constraint rather than a fully
         coupled bordered solve.  Each iteration still solves K du = R at
         the current lambda, then moves lambda to drive

             g = dG - dtau = 0,    dG = (P_n*lambda - P*lambda_n)/2

         with the slope taken by secant from the previous iteration.  It
         converges more slowly than the coupled form, but it needs no
         second right-hand side and no f_hat vector, and it has the one
         property that matters here: lambda is free to decrease.

         xbounact is rebuilt directly, because tempload composes it as
         xbounold + (xboun-xbounold)*reltime and resultsini.c writes that
         value straight into v.  theta stays untouched until the increment
         is accepted; only then is dtheta set so checkconvergence lands on
         the lambda that was actually reached. */

      /* The constraint may only be driven by an equilibrium reaction.  In
         the coupled form that is automatic, because lambda and u are
         solved together and both hold at convergence.  Staggered, the
         early iterations carry a reaction that is not in equilibrium yet,
         and driving lambda with it produced a negative dG - impossible for
         a dissipation - and threw the load factor about near the turning
         point.  So lambda is held until the residual has contracted for
         two iterations, and a non-positive dG is treated as "not yet
         meaningful" rather than as a constraint violation. */

      damage_diss_ok=0;
      if((damage_diss_ctrl==1)&&(iit>=3)&&(damage_diss_init==1)&&
         (ram[0]<ram1[0])&&(ram1[0]<ram2[0])){
        damage_diss_ok=1;
      }

      if(damage_diss_ok==1){
        damage_diss_dgcur=0.5*(damage_diss_pprev*damage_diss_lamcur-
                               damage_diss_p*damage_diss_lprev);
        if(damage_diss_dgcur<=0.) damage_diss_ok=0;
      }

      if(damage_diss_ok==1){
        damage_diss_g=damage_diss_dgcur-damage_diss_target;

        if(damage_diss_have==1){
          damage_diss_slope=(damage_diss_dgcur-damage_diss_dgold)/
            (damage_diss_lamcur-damage_diss_lamold);
        }else{
          /* first pass: the only slope estimate available is the secant
             through the converged state */
          damage_diss_slope=0.5*damage_diss_pprev;
        }
        if(fabs(damage_diss_slope)<1.e-30) damage_diss_slope=
          (damage_diss_slope<0.)?-1.e-30:1.e-30;

        damage_diss_dgold=damage_diss_dgcur;
        damage_diss_lamold=damage_diss_lamcur;
        damage_diss_have=1;

        damage_diss_dlam=-damage_diss_g/damage_diss_slope;

        /* bound the move so a bad secant cannot throw lambda across the
           whole step */
        if(damage_diss_dlam>DAMAGE_DISS_DLAM*dthetaref)
          damage_diss_dlam=DAMAGE_DISS_DLAM*dthetaref;
        if(damage_diss_dlam<-DAMAGE_DISS_DLAM*dthetaref)
          damage_diss_dlam=-DAMAGE_DISS_DLAM*dthetaref;

        damage_diss_lamcur+=damage_diss_dlam;
        if(damage_arc==1){
          if(damage_diss_lamcur<0.) damage_diss_lamcur=0.;
        }else if(damage_diss_lamcur<theta){
          damage_diss_lamcur=theta;
        }
        if(damage_diss_lamcur>1.) damage_diss_lamcur=1.;

        for(i=0;i<*nboun;i++){
          xbounact[i]=xbounold[i]+
            (xboun[i]-xbounold[i])*damage_diss_lamcur;
        }
      }

      /* updating the external work (only for dynamic calculations) */

      if((*nmethod==4)&&(*ithermal<2)&&(*nener==1)){
	allwk=allwkini;
	worparll(&allwk,fnext,&mt,fnextini,v,vini,nk,&num_cpus);

        /* Work due to damping forces (cv and cvini) --> MPADD */

	if(idamping==1){
	  dampwk=dampwkini;
	  dam1parll(&mt,nactdof,aux2,v,vini,nk,&num_cpus);
	  dam2parll(&dampwk,cv,cvini,aux2,&neq[0],&num_cpus);
	}
        /* Damping forces --> MPADD */
      }

      /* line search (only for static surface-to-surface penalty contact)
         and not in the first iteration */

      if((*mortar==1)&&(iit!=1)&&(*ne-ne0>0)&&(*nmethod!=4)){

	SFREE(v);SFREE(stx);SFREE(fn);
      
	/* calculating the residual */
      
	NNEW(res,double,neq[1]);
	calcresidual(nmethod,neq,res,fext,f,iexpl,nactdof,aux2,vold,vini,
		     &dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
		     islavnode,nslavnode,mortar,ntie,mi,nzs,&nasym,
		     &idamping,veold,adc,auc,cvini,cv,&alpham,&num_cpus);

	/* calculating the line search factor */

	sum1=0.;sum2=0.;
	for(i=0;i<neq[1];i++){
	  sum1+=b[i]*resold[i];
	  sum2+=b[i]*res[i];
	}
	SFREE(res);

	if(fabs(sum1-sum2)<1.e-30){
	  flinesearch=1.;
	}else{
	  flinesearch=sum1/(sum1-sum2);
	  if(flinesearch>smaxls){
	    flinesearch=smaxls;
	  }else if(flinesearch<sminls){
	    flinesearch=sminls;
	  }
	}
	printf("line search factor=%f\n\n",flinesearch);

	/* update the solution */

	for(i=0;i<neq[1];i++){
	  b[i]*=flinesearch;}
      
	MNEW(v,double,mt**nk);
	isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
	  
	NNEW(stx,double,6*mi[0]**ne);
	MNEW(fn,double,mt**nk);
	  
	if(ne1d2d==1)NNEW(inum,ITG,*nk);
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		ndirboun,xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
		&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
		emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
		iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
		fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		&reltime,&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
		inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavquadel,aut,irowt,jqt,&mortartrafoflag,
		&intscheme,physcon,dam,damn,iponoel);
	if(ne1d2d==1)SFREE(inum);

	/* Structural finite-difference check of the ASSEMBLED Jacobian.

	   mpfd validates the 6x6 material tangent and says it is right to
	   1e-3 even on a shear path.  Nothing validates what the element
	   loop and add_sm_st_as then make of it, and that is the only step
	   left between a correct constitutive law and the operator PARDISO
	   actually factorises.

	   K(i,j) = d fn_i / d u_j is measured by central differences on the
	   internal force, one column per perturbed degree of freedom, and
	   compared against the coefficient read back out of the sparse
	   structure.  Columns are taken at the nodes of the most damaged
	   element, because a tangent that is right in the elastic bulk and
	   wrong in the process zone is exactly the failure that would
	   survive every test written so far.

	   Storage convention, from add_sm_st_as.f: ad(i) is the diagonal;
	   K(i,j) with i>j sits in column j of au; with i<j it sits in
	   column i offset by nzs(3).

	   The probe perturbs v and calls results repeatedly, so the run is
	   stopped as soon as it reports.  This is a diagnostic build, not a
	   production one. */

      }

      /* BK3 adaptive damage line search.

         Evaluate the full Newton trial first.  A line search is admissible
         only in a contact-free, static physical increment, only while
         progressive damage is actually increasing from the committed
         baseline, and only when the mechanical infinity norm of the residual
         grows by more than DAMAGE_LINESEARCH_GROWTH.  Contracting Newton
         iterations therefore pay no extra constitutive evaluation and are
         never damped.  Terminal same-load topology solves stay on NC2. */



      /* ---- [DAMAGE TR LINCHECK] the sign convention, measured -----------

         Everything the trust region does rests on three conventions that are
         easy to get backwards and impossible to see in the output:

             b (before the solve) = fext - f = -R,   R = f_int - f_ext
             p_N = +b (after the solve),  so  J p_N = r0 = -R
             d   = J^T r0 = -g,  the steepest-DESCENT direction for phi

         This block measures them instead of asserting them.  On a healthy
         increment - the state must be smooth for the test to mean anything -
         it walks eps down a ladder and compares the residual actually
         returned by results()/calcresidual() against the linear model the
         dogleg uses:

             R(u + eps*p)  ->  R(u) + eps*J p        i.e.
             res(u+eps*p)  ->  r0 - eps*(J p)        (res = -R throughout)

         DEFECT = |res(u+eps p) - (r0 - eps Jp)| / (eps |Jp|)  must fall like
         O(eps).  It is run twice: once with p = p_N, where J p = r0 exactly
         and the ratio |res|/((1-eps)|r0|) must go to 1 - this pins the sign
         of b, of R and of the Newton correction together - and once with
         p = d scaled to |p_N|, where J p = w.  The second pass is the one
         that pins J^T: d came out of the transpose loop and w out of the
         forward loop, so a swapped pair cannot pass both.

         It ends by evaluating p_N, which is the state the unprobed code
         would have had, and puts cam/qa/uam back first. */

      if((damage_dl_lc_due==1)&&(damage_dl_have==1)){
        /* [WALLDIAG] The ladder has to reach BELOW the step the search
           actually takes, or it cannot tell a consistent tangent with a
           small radius of validity from an inconsistent one: the second
           wall's best rung is alpha~0.004, and the old ladder stopped at
           0.0156.  1/2^k down to k=14 puts six rungs below 0.004 while the
           defect is still far above the cancellation floor, which for
           |J p|=|r0| sits at about 1e-16/eps. */
        static const double lcE[15]={1.,0.5,0.25,0.125,0.0625,0.03125,
                                     0.015625,0.0078125,0.00390625,
                                     0.001953125,0.0009765625,
                                     0.00048828125,0.000244140625,
                                     0.0001220703125,0.00006103515625};
        ITG lnst,lii,ljj,lpass,lact,lspc,lmpcd;
        double le,lsc,lnum,lden,lnjp,lnres,lmv,ldv,linf,linf0;
        double lqas[4],luams[2];
        double *lp=NULL,*ljp=NULL;

        lnst=*nstate_;
        if(damage_dl_res==NULL){
          NNEW(damage_dl_res,double,neq[1]);
          NNEW(damage_dl_dam,double,mi[0]**ne);
          NNEW(damage_dl_visc,double,mi[0]**ne);
          if(lnst>0) NNEW(damage_dl_xs,double,lnst*mi[0]**ne);
        }
        NNEW(lp,double,neq[1]);
        NNEW(ljp,double,neq[1]);

        lact=0;lspc=0;lmpcd=0;
        for(ljj=0;ljj<*nk;ljj++){
          for(lii=1;lii<mt;lii++){
            if(nactdof[mt*ljj+lii]>0) lact++;
            else if(nactdof[mt*ljj+lii]==0) lspc++;
            else lmpcd++;   /* negative: MPC-dependent, and in this tree
                               SPC-eliminated slots land here as well */
          }
        }
        printf("[DAMAGE TR LINCHECK] inc=%" ITGFORMAT " iter=%" ITGFORMAT
               " REDUCED SPACE: neq[0]=%" ITGFORMAT " neq[1]=%" ITGFORMAT
               "; mechanical slots %" ITGFORMAT " = active (nactdof>0) %"
               ITGFORMAT " + excluded, nactdof==0 %" ITGFORMAT
               " + excluded, nactdof<0 %" ITGFORMAT "; nboun=%" ITGFORMAT
               " nmpc=%" ITGFORMAT " nzs[1]=%" ITGFORMAT " nzs[2]=%" ITGFORMAT
               ".  active == neq[1] and excluded == nboun+MPC-dependent is "
               "the proof that b, f, fext, ad and au share ONE post-SPC/MPC "
               "numbering (nactdof>0, 1-based), so every vector the trust "
               "region forms already lives in the reduced space and no "
               "constrained degree of freedom is ever stepped.%s",
               iinc,iit,neq[0],neq[1],3*(*nk),lact,lspc,lmpcd,*nboun,*nmpc,
               nzs[1],nzs[2],"\n");
        printf("[DAMAGE TR LINCHECK] transpose identity dot(J^T R,p_N)/|R|^2 "
               "= %.12e (exact 1); operator asymmetry "
               "max|au_L-au_U|/max|au_L| = %.6e; |p_N|=%.6e |p_C|=%.6e "
               "|R|2=%.6e%s",
               damage_dl_ident,damage_dl_asym,sqrt(damage_dl_npn2),
               damage_dl_tc*sqrt(damage_dl_nd2),sqrt(damage_dl_nb2),"\n");
        fflush(stdout);

        isiz=mi[0]**ne;cpypardou(damage_dl_dam,dam,&isiz,&num_cpus);
        if(damage_damvisc!=NULL){
          isiz=mi[0]**ne;
          cpypardou(damage_dl_visc,damage_damvisc,&isiz,&num_cpus);
        }
        if((lnst>0)&&(damage_dl_xs!=NULL)){
          isiz=lnst*mi[0]**ne;
          cpypardou(damage_dl_xs,xstate,&isiz,&num_cpus);
        }
        for(ljj=0;ljj<4;ljj++) lqas[ljj]=qa[ljj];
        for(ljj=0;ljj<2;ljj++) luams[ljj]=uam[ljj];

        /* [WALLDIAG] The active set AT THE CURRENT ITERATE.  xstate, dam and
           stx here belong to u, because the iteration's own results() call
           built them and nothing has stepped yet.  Every rung below is
           compared against this one census, so a transition count is the
           number of integration points the step moved across a branch. */

        if(damage_wall_cat==NULL) NNEW(damage_wall_cat,ITG,mi[0]**ne);
        damage_ray_census(damage_wall_cat,xstate,xstateini,dam,damdamageini,
                          damage_damvisc,stx,ipkon,lakon,ne0,mi[0],*nstate_);
        {
          ITG wj,wnadv=0,wnlive=0;
          double wpinf=0.;
          for(wj=0;wj<mi[0]*ne0;wj++){
            if(damage_wall_cat[wj]&DAMCAT_USOFT) wnlive++;
            if(damage_wall_cat[wj]&DAMCAT_UADV) wnadv++;
          }
          for(wj=0;wj<neq[1];wj++)
            if(fabs(damage_dl_pn[wj])>wpinf) wpinf=fabs(damage_dl_pn[wj]);
          printf("[WALLDIAG] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                 " base state: UC6 points past initiation %" ITGFORMAT
                 ", of which ADVANCING (deff>dmax0) %" ITGFORMAT
                 "; |p_N|inf=%.6e |p_N|2=%.6e |R|2=%.6e%s",
                 iinc,iit,wnlive,wnadv,wpinf,sqrt(damage_dl_npn2),
                 sqrt(damage_dl_nb2),"\n");
          damage_wall_where("residual",damage_dl_r0,neq[1],nactdof,mt,*nk,5,
                            ipkon,kon,lakon,xstate,stx,dam,*ne,ne0,mi[0],
                            *nstate_);
          damage_wall_where("correction",damage_dl_pn,neq[1],nactdof,mt,*nk,
                            5,ipkon,kon,lakon,xstate,stx,dam,*ne,ne0,mi[0],
                            *nstate_);
          /* [WALLDIAG] STIFFNESS AT THE RESIDUAL PEAK.  The convergence test
             checkconvergence() applies is on max|R| over the mechanical
             block, not on |R|2, so the dof that decides the run is the peak
             one.  This reports, for the five largest residual dofs, the
             assembled diagonal AGAINST that node's own intact value and the
             displacement |R|/k that would be needed to null the residual
             locally.  A peak sitting on a node whose diagonal has collapsed,
             needing a displacement far larger than anything physical, is a
             different object from a peak on a healthy node, and only the
             second is a convergence problem in the ordinary sense. */
          {
            ITG *wsn=NULL,*wsd=NULL,wi,wj,wk,wt,wbest;
            double wa;
            NNEW(wsn,ITG,neq[1]);NNEW(wsd,ITG,neq[1]);
            for(wi=0;wi<neq[1];wi++){wsn[wi]=-1;wsd[wi]=0;}
            for(wi=0;wi<*nk;wi++)
              for(wj=1;wj<mt;wj++){
                wk=nactdof[mt*wi+wj];
                if((wk>0)&&(wk<=neq[1])){wsn[wk-1]=wi;wsd[wk-1]=wj;}
              }
            for(wt=0;wt<5;wt++){
              wbest=-1;wa=-1.;
              for(wi=0;wi<neq[1];wi++){
                if(wsn[wi]<0) continue;
                if(fabs(damage_dl_r0[wi])>wa){wa=fabs(damage_dl_r0[wi]);wbest=wi;}
              }
              if(wbest<0) break;
              wi=wsn[wbest];
              printf("[WALLDIAG]   Rpeak #%" ITGFORMAT ": node %" ITGFORMAT
                     " dir %" ITGFORMAT " R=%.6e  addiag=%.6e addiag0=%.6e "
                     "ratio=%.6e spc_masked=%" ITGFORMAT " need_du=|R|/k=%.6e"
                     "%s",wt+1,wi+1,wsd[wbest],damage_dl_r0[wbest],
                     (damage_addiag!=NULL)?damage_addiag[wi]:0.,
                     (damage_addiag0!=NULL)?damage_addiag0[wi]:0.,
                     ((damage_addiag!=NULL)&&(damage_addiag0!=NULL)&&
                      (damage_addiag0[wi]>0.))?
                       damage_addiag[wi]/damage_addiag0[wi]:-1.,
                     ((damage_spc_mask!=NULL)&&(damage_spc_nk>wi))?
                       damage_spc_mask[wi]:-1,
                     ((damage_addiag!=NULL)&&(damage_addiag[wi]>0.))?
                       fabs(damage_dl_r0[wbest])/damage_addiag[wi]:-1.,"\n");
              wsn[wbest]=-1;
            }
            SFREE(wsn);SFREE(wsd);
          }
          fflush(stdout);
        }

        {
        ITG lpassn=2;
        if((damage_wall_maskstep!=0)&&(damage_dl_pm!=NULL)) lpassn=3;
        for(lpass=0;lpass<lpassn;lpass++){
          if(lpass==0){
            for(ljj=0;ljj<neq[1];ljj++){
              lp[ljj]=damage_dl_pn[ljj];
              ljp[ljj]=damage_dl_r0[ljj];
            }
            printf("[DAMAGE TR LINCHECK] pass 1: p = p_N, so J p = r0 "
                   "EXACTLY.  res(u+eps p)/((1-eps)|r0|) must go to 1 and "
                   "the defect to zero like O(eps).%s","\n");
          }else if(lpass==1){
            lsc=sqrt(damage_dl_npn2/damage_dl_nd2);
            for(ljj=0;ljj<neq[1];ljj++){
              lp[ljj]=lsc*damage_dl_d[ljj];
              ljp[ljj]=lsc*damage_dl_w[ljj];
            }
            printf("[DAMAGE TR LINCHECK] pass 2: p = %.6e * d (scaled to "
                   "|p_N|), so J p = %.6e * w.  d comes from the TRANSPOSE "
                   "loop and w from the forward loop, so a swapped pair "
                   "cannot pass this.%s",lsc,lsc,"\n");
          }else{
            /* [WALLDIAG] MASKSTEP.  p = p_N with the AUTOSPC-masked nodes'
               components zeroed.  UNSCALED on purpose: the question is not
               how this direction behaves at |p_N|, it is what the step
               actually does once the collapsed-diagonal dofs are taken out
               of it, so eps=1 here means "the rest of the Newton step". */
            for(ljj=0;ljj<neq[1];ljj++){
              lp[ljj]=damage_dl_pm[ljj];
              ljp[ljj]=damage_dl_wm[ljj];
            }
            printf("[DAMAGE TR LINCHECK] pass 3 (MASKSTEP): p = p_N with the "
                   "%" ITGFORMAT " AUTOSPC-masked nodes zeroed.  |p|=%.6e "
                   "against |p_N|=%.6e, so the mask removes %.4f%% of "
                   "|p_N|^2.  If the best rung here beats pass 1's best rung "
                   "(and the ladder's accepted alpha), the step is unusable "
                   "BECAUSE it is spent on collapsed-diagonal nodes; if it "
                   "does not, that hypothesis is dead.%s",
                   damage_spc_count,sqrt(damage_dl_npm2),
                   sqrt(damage_dl_npn2),
                   (damage_dl_npn2>0.)?
                     100.*(1.-damage_dl_npm2/damage_dl_npn2):0.,"\n");
          }
          fflush(stdout);
          lnjp=0.;
          for(ljj=0;ljj<neq[1];ljj++) lnjp+=ljp[ljj]*ljp[ljj];
          lnjp=sqrt(lnjp);
          for(lii=0;lii<15;lii++){
            le=lcE[lii];
          /* [DAMAGE TR LINCHECK] cam starts clean for every probe, for the
             same reason as in the trust-region block. */
          for(ljj=0;ljj<3;ljj++) cam[ljj]=0.;
          for(ljj=3;ljj<5;ljj++) cam[ljj]=0.5;
          isiz=mi[0]**ne;cpypardou(dam,damage_dl_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_dl_visc,&isiz,&num_cpus);
          }
          if((lnst>0)&&(damage_dl_xs!=NULL)){
            isiz=lnst*mi[0]**ne;
            cpypardou(xstate,damage_dl_xs,&isiz,&num_cpus);
          }
          SFREE(v);SFREE(stx);SFREE(fn);
          for(ljj=0;ljj<neq[1];ljj++) b[ljj]=le*lp[ljj];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_dl_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);

            damage_dl_neval++;
            lnum=0.;lnres=0.;linf=0.;linf0=0.;
            if(damage_wall_def==NULL) NNEW(damage_wall_def,double,neq[1]);
            for(ljj=0;ljj<neq[1];ljj++){
              lmv=damage_dl_r0[ljj]-le*ljp[ljj];
              ldv=damage_dl_res[ljj]-lmv;
              damage_wall_def[ljj]=ldv;
              lnum+=ldv*ldv;
              lnres+=damage_dl_res[ljj]*damage_dl_res[ljj];
            }
            /* [WALLDIAG] The line search judges a trial by max|res| over the
               MECHANICAL block neq[0] (nonlingeo.c, damage_linesearch_*norm),
               while the direction it damps is only guaranteed to descend the
               2-norm merit - the transpose identity dot(J^T R,p_N)=|R|^2 is
               what makes that guarantee exact.  The two norms are therefore
               measured on the SAME rung here, or the disagreement between
               them stays an inference. */
            for(ljj=0;ljj<neq[0];ljj++){
              if(fabs(damage_dl_res[ljj])>linf) linf=fabs(damage_dl_res[ljj]);
              if(fabs(damage_dl_r0[ljj])>linf0) linf0=fabs(damage_dl_r0[ljj]);
            }
            lnum=sqrt(lnum);lnres=sqrt(lnres);
            lden=le*lnjp;
            printf("[DAMAGE TR LINCHECK] pass %" ITGFORMAT " eps=%.9f "
                   "|res(u+eps p)|=%.12e  linear-model defect=%.6e  "
                   "defect/eps=%.6e%s",lpass+1,le,lnres,
                   (lden>0.)?lnum/lden:0.,
                   (lden>0.)?lnum/(lden*le):0.,"\n");
            printf("[WALLDIAG] pass %" ITGFORMAT " eps=%.9f  TWO NORMS: "
                   "|res|2/|r0|2=%.12f  |res|inf/|r0|inf=%.12f  "
                   "(|r0|2=%.6e |r0|inf=%.6e)%s",lpass+1,le,
                   (damage_dl_nb2>0.)?lnres/sqrt(damage_dl_nb2):0.,
                   (linf0>0.)?linf/linf0:0.,sqrt(damage_dl_nb2),linf0,"\n");
            if(lpass==0)
              printf("[DAMAGE TR LINCHECK]          ratio "
                     "|res|/((1-eps)|r0|) = %.12f  (must go to 1)%s",
                     ((1.-le)>0.)?lnres/((1.-le)*sqrt(damage_dl_nb2)):0.,
                     "\n");
            {
              ITG wtot;
              wtot=damage_wall_setdiff(damage_wall_cat,xstate,xstateini,dam,
                                       damdamageini,damage_damvisc,stx,
                                       ipkon,lakon,ne0,mi[0],*nstate_,
                                       damage_wall_nb);
              printf("[WALLDIAG] pass %" ITGFORMAT " eps=%.9f  active-set "
                     "transitions %" ITGFORMAT " = UC6 loading/unloading %"
                     ITGFORMAT " + UC6 initiation %" ITGFORMAT
                     " + UC6 viscous %" ITGFORMAT " + UC6 failure %"
                     ITGFORMAT " + UC6 tension/compression %" ITGFORMAT
                     " + bulk plastic %" ITGFORMAT " + bulk initiation %"
                     ITGFORMAT " + bulk damage growth %" ITGFORMAT "%s",
                     lpass+1,le,wtot,damage_wall_nb[6],damage_wall_nb[3],
                     damage_wall_nb[4],damage_wall_nb[5],damage_wall_nb[7],
                     damage_wall_nb[0],damage_wall_nb[1],damage_wall_nb[2],
                     "\n");
            }
            /* [WALLDIAG] IS THE MISSING TERM THE DAMAGE RANK-1 TERM?

               J = K0 + E, where E is exactly what mafilldamas.f adds from
               damjac.  Assemble E ALONE into a zeroed pair and apply it to
               p_N.  If the true operator is A = K0 + (1+c)E - that is, if
               the rank-1 damage term is the right shape and the wrong size -
               then (J-A)p = -c*E*p, so the measured defect must be
               ANTI-PARALLEL to E*p and |defect|/(eps*|E p|) must equal |c|.
               A cosine of zero says the missing term is a different term,
               and no scaling of this one can supply it. */
            if(lii==14){
              ITG e1i,e1c,e1k,e1r;
              double *e1ad=NULL,*e1au=NULL,*e1y=NULL,e1n=0.,e1d=0.,e1w=0.;
              ITG e1nd,e1sk,e1ad2,e1ho,e1fl;
              NNEW(e1ad,double,neq[1]);
              NNEW(e1au,double,(nasym+1)*nzs[1]);
              NNEW(e1y,double,neq[1]);
              e1nd=0;e1sk=0;e1ad2=0;e1ho=0;e1fl=0;
              FORTRAN(mafilldamas,(co,kon,ipkon,lakon,&ne0,nactdof,jq,irow,
                                   neq,nzs,e1au,e1ad,vold,mi,damage_damjac,
                                   nmpc,&e1nd,dam,damdamageini,
                                   &e1sk,&e1ad2,&e1ho,&e1fl));
              for(e1k=0;e1k<neq[1];e1k++)
                e1y[e1k]=e1ad[e1k]*damage_dl_pn[e1k];
              for(e1c=0;e1c<neq[1];e1c++){
                for(e1k=jq[e1c]-1;e1k<jq[e1c+1]-1;e1k++){
                  e1r=irow[e1k]-1;
                  e1y[e1r]+=e1au[e1k]*damage_dl_pn[e1c];
                  e1y[e1c]+=e1au[nzs[2]+e1k]*damage_dl_pn[e1r];
                }
              }
              for(e1i=0;e1i<neq[1];e1i++){
                e1n+=e1y[e1i]*e1y[e1i];
                e1d+=damage_wall_def[e1i]*damage_wall_def[e1i];
                e1w+=e1y[e1i]*damage_wall_def[e1i];
              }
              e1n=sqrt(e1n);e1d=sqrt(e1d);
              /* Is the rank-1 term small because dD/d(eps) is small, or
                 because the assembly loses it?  damjac slots 1..6 are the
                 effective stress and 7..12 are dD/d(eps); censusing both
                 separates "the derivative is tiny" from "the derivative is
                 there and the operator does not carry it". */
              {
                ITG qi,qn=0,qz=0;
                double qmax=0.,qsum=0.,tmax=0.;
                for(qi=0;qi<ne0;qi++){
                  ITG qk;double qa2=0.,qt2=0.;
                  if(ipkon[qi]<0) continue;
                  if(lakon[8*qi]!='C') continue;
                  for(qk=6;qk<12;qk++)
                    qa2+=damage_damjac[12*mi[0]*qi+qk]
                        *damage_damjac[12*mi[0]*qi+qk];
                  for(qk=0;qk<6;qk++)
                    qt2+=damage_damjac[12*mi[0]*qi+qk]
                        *damage_damjac[12*mi[0]*qi+qk];
                  if((qa2<=0.)&&(qt2<=0.)) continue;
                  qn++;
                  qa2=sqrt(qa2);qt2=sqrt(qt2);
                  if(qa2<=1.e-12) qz++;
                  if(qa2>qmax) qmax=qa2;
                  if(qt2>tmax) tmax=qt2;
                  qsum+=qa2;
                }
                printf("[WALLDIAG]   damjac census over %" ITGFORMAT
                       " elements with any entry: |dD/d(eps)| max=%.6e "
                       "mean=%.6e, of which %" ITGFORMAT
                       " have it BELOW 1e-12 (assembled anyway, because the "
                       "skip test sums the stress slots too); |sigma_eff| "
                       "max=%.6e%s",qn,qmax,(qn>0)?qsum/qn:0.,qz,tmax,"\n");
              }
              printf("[WALLDIAG]   damage rank-1 term applied to p_N: "
                     "elements %" ITGFORMAT ", |E p|=%.6e, |defect|/eps=%.6e,"
                     "  cos(defect,E p)=%+.6f  |defect|/(eps|E p|)=%.6f%s",
                     e1nd,e1n,e1d/le,
                     ((e1n>0.)&&(e1d>0.))?e1w/(e1n*e1d):0.,
                     (e1n>0.)?e1d/(le*e1n):0.,"\n");
              SFREE(e1ad);SFREE(e1au);SFREE(e1y);
            }
            if(lii==14){
              damage_wall_split("linear-model defect",damage_wall_def,neq[1],
                                nactdof,mt,*nk,ipkon,kon,lakon,dam,xstate,
                                *ne,ne0,mi[0],*nstate_);
              damage_wall_split("residual r0",damage_dl_r0,neq[1],nactdof,mt,
                                *nk,ipkon,kon,lakon,dam,xstate,*ne,ne0,mi[0],
                                *nstate_);
              damage_wall_split("Newton step p_N",damage_dl_pn,neq[1],nactdof,
                                mt,*nk,ipkon,kon,lakon,dam,xstate,*ne,ne0,
                                mi[0],*nstate_);
              damage_wall_where("defect",damage_wall_def,neq[1],nactdof,mt,
                                *nk,5,ipkon,kon,lakon,xstate,stx,dam,*ne,ne0,
                                mi[0],*nstate_);
            }
            fflush(stdout);
          }
        }
        }

        /* leave the full Newton step, exactly as the unprobed code would */
        for(ljj=0;ljj<4;ljj++) qa[ljj]=lqas[ljj];
        for(ljj=0;ljj<2;ljj++) uam[ljj]=luams[ljj];
        le=1.;
        for(ljj=0;ljj<neq[1];ljj++) lp[ljj]=damage_dl_pn[ljj];
          /* [DAMAGE TR LINCHECK] cam starts clean for every probe, for the
             same reason as in the trust-region block. */
          for(ljj=0;ljj<3;ljj++) cam[ljj]=0.;
          for(ljj=3;ljj<5;ljj++) cam[ljj]=0.5;
          isiz=mi[0]**ne;cpypardou(dam,damage_dl_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_dl_visc,&isiz,&num_cpus);
          }
          if((lnst>0)&&(damage_dl_xs!=NULL)){
            isiz=lnst*mi[0]**ne;
            cpypardou(xstate,damage_dl_xs,&isiz,&num_cpus);
          }
          SFREE(v);SFREE(stx);SFREE(fn);
          for(ljj=0;ljj<neq[1];ljj++) b[ljj]=le*lp[ljj];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_dl_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);

        lnres=0.;
        for(ljj=0;ljj<neq[1];ljj++)
          lnres+=damage_dl_res[ljj]*damage_dl_res[ljj];
        printf("[DAMAGE TR LINCHECK] restored to the full Newton step; "
               "|res| there = %.12e.  Compare with the eps=1 line of pass 1: "
               "equal means every probe rolled the mutable state back "
               "completely.%s",sqrt(lnres),"\n");
        fflush(stdout);
        SFREE(lp);SFREE(ljp);
        damage_dl_lc_due=0;
      }


      /* ============ [DAMAGE CT] bordered corrector ==================
         One factorisation per iteration, two solves against it:
             K z = -R = b        (b as the solver left it)
             K y = -q            (-q straight from two b vectors)
             den = c_lambda + c_u^T y
             dlambda = (-c - c_u^T z)/den ,   du = z + y*dlambda
         The factorisation has ALREADY been taken from the base state above,
         so the perturbation evaluations below cannot contaminate it: the
         cached LU is what pardiso_solve uses for y.
         c(u,lambda) = m.(delta - delta_c) - ds is exactly linear in u with
         a frozen m, so c is driven to zero every iteration by construction
         and its value is reported, never assumed. */

      damage_ct_used=0;
      if((damage_ct_on==1)&&(idamagereeq==0)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)&&
         (*mortar<=1)){
        ITG ctj,ctk,ctnst,ctbad=0,ctnsw=0;
        damage_ct_used=1;
        double ctlam,ctc,ctcuz,ctcuy,ctden,ctdlam,ctrho,ctyinf;
        double ctdl[3],ctrm[9],ctsh[3],ctnrm,ctb0;

        ctnst=*nstate_;
        damage_ct_it++;damage_ct_ncorr++;
        /* z is what the solver returned in b.  The base snapshot was taken
           BEFORE the solve, at this iterate; it is only restored here. */
        isiz=neq[1];cpypardou(damage_ct_z,b,&isiz,&num_cpus);

        /* ---- q by transactional finite difference at the CURRENT u ----
           arc_r0 and the whole read-before-write set were captured BEFORE
           the solve, at this iterate, so every probe starts from ONE base
           state.  The ladder runs on a step's first corrector iteration and
           whenever a probe reports a category switch; the ACCEPTED q is the
           refined q(eps_final). */
        {
          ITG cthalve,ctacc=0,cte,ctp;
          double cteps,ctqn,ctdn,ctbase;
          cteps=damage_ct_eps;
          if(damage_ct_qh==NULL) NNEW(damage_ct_qh,double,neq[1]);

          /* unperturbed baseline through the probe path: it fixes the
             category map AND proves the base state is reproduced */
          ctlam=damage_ct_lam;
          for(ctj=0;ctj<neq[1];ctj++) b[ctj]=0.;
          for(ctj=0;ctj<4;ctj++) qa[ctj]=damage_ct_qa[ctj];
          for(ctj=0;ctj<5;ctj++) cam[ctj]=damage_ct_cam[ctj];
          for(ctj=0;ctj<2;ctj++) uam[ctj]=damage_ct_uam[ctj];
          isiz=mi[0]**ne;cpypardou(dam,damage_ct_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_ct_visc,&isiz,&num_cpus);
          }
          if((ctnst>0)&&(damage_ct_xs!=NULL)){
            isiz=ctnst*mi[0]**ne;
            cpypardou(xstate,damage_ct_xs,&isiz,&num_cpus);
          }
          if((damage_damjac!=NULL)&&(damage_ct_jac!=NULL)){
            isiz=12*mi[0]**ne;
            cpypardou(damage_damjac,damage_ct_jac,&isiz,&num_cpus);
          }
          for(ctj=0;ctj<*nboun;ctj++)
            xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
          SFREE(v);SFREE(stx);SFREE(fn);
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ct_beps,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          damage_ct_neval++;
          damage_evt_sign(stx,ipkon,lakon,ne0,mi[0],damage_ct_sgn);
          ctbase=0.;
          for(ctj=0;ctj<neq[1];ctj++)
            if(fabs(damage_ct_beps[ctj]-damage_ct_r0[ctj])>ctbase)
              ctbase=fabs(damage_ct_beps[ctj]-damage_ct_r0[ctj]);
          if(damage_ct_it<=1){
            printf("[DAMAGE CT FD] inc=%" ITGFORMAT " baseline reproduced: "
                   "max|R(base)-arc_r0|=%.6e (must be 0; a non-zero value "
                   "means the probe path does not start from the iterate the "
                   "residual was taken at)%s",iinc,ctbase,"\n");
            fflush(stdout);
          }

          for(cthalve=0;cthalve<5;cthalve++){
            ctlam=damage_ct_lam+cteps;
            for(ctj=0;ctj<neq[1];ctj++) b[ctj]=0.;
          for(ctj=0;ctj<4;ctj++) qa[ctj]=damage_ct_qa[ctj];
          for(ctj=0;ctj<5;ctj++) cam[ctj]=damage_ct_cam[ctj];
          for(ctj=0;ctj<2;ctj++) uam[ctj]=damage_ct_uam[ctj];
          isiz=mi[0]**ne;cpypardou(dam,damage_ct_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_ct_visc,&isiz,&num_cpus);
          }
          if((ctnst>0)&&(damage_ct_xs!=NULL)){
            isiz=ctnst*mi[0]**ne;
            cpypardou(xstate,damage_ct_xs,&isiz,&num_cpus);
          }
          if((damage_damjac!=NULL)&&(damage_ct_jac!=NULL)){
            isiz=12*mi[0]**ne;
            cpypardou(damage_damjac,damage_ct_jac,&isiz,&num_cpus);
          }
          for(ctj=0;ctj<*nboun;ctj++)
            xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
          SFREE(v);SFREE(stx);SFREE(fn);
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ct_beps,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          damage_ct_neval++;

            ctnsw=damage_evt_flips(stx,ipkon,lakon,ne0,mi[0],damage_ct_sgn,
                                   &cte,&ctp);
            for(ctj=0;ctj<neq[1];ctj++)
              damage_ct_y[ctj]=(damage_ct_beps[ctj]-damage_ct_r0[ctj])/cteps;
            if(ctnsw!=0){
              printf("[DAMAGE CT FD] inc=%" ITGFORMAT " it=%" ITGFORMAT
                     " eps=%.3e REJECTED: %" ITGFORMAT " UC6 category "
                     "switch(es) between R(lambda) and R(lambda+eps)%s",
                     iinc,damage_ct_it,cteps,ctnsw,"\n");
              cteps*=0.5;continue;
            }
            if((damage_ct_epsok==1)&&(cthalve==0)){ctacc=1;break;}
            if(cthalve>0){
              ctqn=0.;ctdn=0.;
              for(ctj=0;ctj<neq[1];ctj++){
                ctqn+=damage_ct_y[ctj]*damage_ct_y[ctj];
                ctdn+=(damage_ct_y[ctj]-damage_ct_qh[ctj])
                     *(damage_ct_y[ctj]-damage_ct_qh[ctj]);
              }
              ctqn=sqrt(ctqn);ctdn=sqrt(ctdn);
              printf("[DAMAGE CT FD] inc=%" ITGFORMAT " it=%" ITGFORMAT
                     " eps=%.3e |q|=%.6e rel.change vs 2eps=%.6e %s%s",
                     iinc,damage_ct_it,cteps,ctqn,(ctqn>0.)?ctdn/ctqn:0.,
                     ((ctqn>0.)&&(ctdn/ctqn<0.1))?
                     "ACCEPT (refined q used)":"halve again","\n");
              if((ctqn>0.)&&(ctdn/ctqn<0.1)){
                ctacc=1;damage_ct_eps=cteps;damage_ct_epsok=1;break;
              }
            }
            isiz=neq[1];cpypardou(damage_ct_qh,damage_ct_y,&isiz,&num_cpus);
            cteps*=0.5;
          }
          fflush(stdout);
          if(ctacc==0) ctbad=9;
        }
#ifdef PARDISO
        pardiso_solve(damage_ct_y,&neq[0],&symmetryflag,&inputformat,&nrhs);
#endif

        /* ---- the constraint at the current iterate ----
           delta_c is anchored HERE, at the first corrector iteration of each
           continuation step, from the same vold the corrector itself reads.
           Anchoring it at arming instead left c off by 8e-10 against
           ds=9.5e-12 on the short deck even though the ring endpoint and
           vold agreed exactly at arming: the state moves between the load
           build and the first corrector call.  With this anchor c = -ds
           holds to machine precision at the start of every step, by
           construction rather than by argument. */
        damage_ct_kin(co,kon,ipkon[damage_ct_elem],vold,mt,damage_ct_ip,
                      ctdl,ctrm,ctsh);
        if(damage_ct_newstep==1){
          damage_ct_newstep=0;
          for(ctj=0;ctj<3;ctj++) damage_ct_dc[ctj]=ctdl[ctj];
          damage_ct_lamc=damage_ct_lam;
          damage_ct_epsok=0;
        }
        ctc=damage_ct_m[0]*(ctdl[0]-damage_ct_dc[0])
           +damage_ct_m[1]*(ctdl[1]-damage_ct_dc[1])
           +damage_ct_m[2]*(ctdl[2]-damage_ct_dc[2])-damage_ct_ds;
        ctcuz=0.;ctcuy=0.;ctnrm=0.;ctyinf=0.;
        for(ctj=0;ctj<damage_ct_nsupp;ctj++){
          ctk=damage_ct_supp[ctj];
          if(ctk<0) continue;
          ctcuz+=damage_ct_w[ctj]*damage_ct_z[ctk];
          ctcuy+=damage_ct_w[ctj]*damage_ct_y[ctk];
          ctnrm+=damage_ct_w[ctj]*damage_ct_w[ctj];
        }
        ctnrm=sqrt(ctnrm);
        ctb0=0.;
        for(ctj=0;ctj<damage_ct_nsupp;ctj++){
          ctk=damage_ct_supp[ctj];
          if(ctk<0) continue;
          ctb0+=damage_ct_y[ctk]*damage_ct_y[ctk];
        }
        ctb0=sqrt(ctb0);
        for(ctj=0;ctj<neq[1];ctj++)
          if(fabs(damage_ct_y[ctj])>ctyinf) ctyinf=fabs(damage_ct_y[ctj]);

        if(damage_ct_bordered(damage_ct_clam,ctcuz,ctcuy,ctc,
                              &ctden,&ctdlam)==0) ctbad=1;
        ctrho=damage_ct_rhoden(damage_ct_clam,ctnrm,ctb0,ctden);
        damage_ct_den=ctden;damage_ct_rho=ctrho;damage_ct_cprev=ctc;

        if((ctbad==0)&&(ctrho<damage_ct_rhomin)) ctbad=2;
        if((ctbad==0)&&
           (fabs(ctdlam)>damage_ct_clim*damage_ct_lamref)) ctbad=3;
        if((ctbad==0)&&
           (fabs(ctdlam)*ctyinf>damage_ct_ulim*damage_ct_duref)) ctbad=4;
        if((ctbad==0)&&(ctnsw!=0)) ctbad=5;
        if((ctbad==0)&&(damage_ct_neval>=damage_ct_maxeval)) ctbad=6;
        if((ctbad==0)&&(damage_ct_nfact>=damage_ct_maxfact)) ctbad=7;
        if((ctbad==0)&&(damage_ct_it>damage_ct_maxcorr)) ctbad=8;

        printf("[DAMAGE CT] inc=%" ITGFORMAT " it=%" ITGFORMAT " step=%"
               ITGFORMAT " lambda=%.12e ds=%.6e c=%.6e den=%.6e rho_den=%.4e "
               "dlambda=%.6e |y|inf=%.4e cuz=%.4e cuy=%.4e switches=%"
               ITGFORMAT " evals=%" ITGFORMAT " fact=%" ITGFORMAT " %s%s",
               iinc,damage_ct_it,damage_ct_step,damage_ct_lam,damage_ct_ds,
               ctc,ctden,ctrho,ctdlam,ctyinf,ctcuz,ctcuy,ctnsw,
               damage_ct_neval,damage_ct_nfact,
               (ctbad==0)?"OK":
               ((ctbad==1)?"REFUSE:degenerate-den":
               ((ctbad==2)?"REFUSE:rho_den":
               ((ctbad==3)?"REFUSE:dlambda-bound":
               ((ctbad==4)?"REFUSE:du-bound":
               ((ctbad==5)?"REFUSE:category-switch-in-FD":
               ((ctbad==6)?"REFUSE:eval-budget":
               ((ctbad==7)?"REFUSE:fact-budget":
               ((ctbad==9)?"REFUSE:FD-eps-ladder":
                           "REFUSE:corrector-limit")))))))),
               "\n");
        fflush(stdout);

        if(ctbad!=0){
          /* A refusal attributable to THIS control point blacklists it for
             the epoch and allows the next wall to try the next candidate;
             a budget or FD refusal is not the point's fault and disarms the
             mechanism outright.  The epoch clears on proven progress. */
          if((ctbad==2)||(ctbad==3)||(ctbad==4)){
            if(damage_ct_bl==NULL) NNEW(damage_ct_bl,ITG,64);
            if(damage_ct_nbl<64){
              damage_ct_bl[damage_ct_nbl++]=4*damage_ct_elem+damage_ct_ip;
              printf("[DAMAGE CT] control point element %" ITGFORMAT " ip %"
                     ITGFORMAT " BLACKLISTED for this epoch (%" ITGFORMAT
                     " blacklisted); the next wall may try another "
                     "candidate.%s",damage_ct_elem+1,damage_ct_ip+1,
                     damage_ct_nbl,"\n");
              damage_ct_refused=0;    /* a re-arm is allowed */
            }else{
              damage_ct_refused=1;
            }
          }else{
            damage_ct_refused=1;
          }
          damage_ct_on=0;
          if(damage_ct_refused==1) damage_ct_partial=1;
          damage_ct_lam=damage_ct_lamc;
          ctlam=damage_ct_lamc;
          for(ctj=0;ctj<neq[1];ctj++) b[ctj]=damage_ct_z[ctj];
          for(ctj=0;ctj<4;ctj++) qa[ctj]=damage_ct_qa[ctj];
          for(ctj=0;ctj<5;ctj++) cam[ctj]=damage_ct_cam[ctj];
          for(ctj=0;ctj<2;ctj++) uam[ctj]=damage_ct_uam[ctj];
          isiz=mi[0]**ne;cpypardou(dam,damage_ct_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_ct_visc,&isiz,&num_cpus);
          }
          if((ctnst>0)&&(damage_ct_xs!=NULL)){
            isiz=ctnst*mi[0]**ne;
            cpypardou(xstate,damage_ct_xs,&isiz,&num_cpus);
          }
          if((damage_damjac!=NULL)&&(damage_ct_jac!=NULL)){
            isiz=12*mi[0]**ne;
            cpypardou(damage_damjac,damage_ct_jac,&isiz,&num_cpus);
          }
          for(ctj=0;ctj<*nboun;ctj++)
            xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
          SFREE(v);SFREE(stx);SFREE(fn);
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ct_beps,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          damage_ct_neval++;

        }else{
          damage_ct_lam+=ctdlam;
          ctlam=damage_ct_lam;
          for(ctj=0;ctj<neq[1];ctj++)
            b[ctj]=damage_ct_z[ctj]+damage_ct_y[ctj]*ctdlam;
          for(ctj=0;ctj<4;ctj++) qa[ctj]=damage_ct_qa[ctj];
          for(ctj=0;ctj<5;ctj++) cam[ctj]=damage_ct_cam[ctj];
          for(ctj=0;ctj<2;ctj++) uam[ctj]=damage_ct_uam[ctj];
          isiz=mi[0]**ne;cpypardou(dam,damage_ct_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_ct_visc,&isiz,&num_cpus);
          }
          if((ctnst>0)&&(damage_ct_xs!=NULL)){
            isiz=ctnst*mi[0]**ne;
            cpypardou(xstate,damage_ct_xs,&isiz,&num_cpus);
          }
          if((damage_damjac!=NULL)&&(damage_ct_jac!=NULL)){
            isiz=12*mi[0]**ne;
            cpypardou(damage_damjac,damage_ct_jac,&isiz,&num_cpus);
          }
          for(ctj=0;ctj<*nboun;ctj++)
            xbounact[ctj]=xbounold[ctj]+(xboun[ctj]-xbounold[ctj])*ctlam;
          SFREE(v);SFREE(stx);SFREE(fn);
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ct_beps,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          damage_ct_neval++;

        }
      }

      /* ---- [DAMAGE TR] TRUST-REGION DOGLEG on the original residual -----

         Armed only as rescue level 3, i.e. only after the Rescue2 levels 1
         and 2 have both failed on the same wall, and only where the measured
         failure is a step-LENGTH failure.  At s3rad inc=569 BK3 reports, on
         all four identical attempts,

             res_old=1.942365e-03  res_full=2.520626e-02  lambda=0.100000
             trials=3  contracted=0  res_damped=2.001658e-03

         so the full Newton step overshoots 13x and the most damped step BK3
         is permitted to take still leaves the residual ABOVE where it
         started, because DAMAGE_LINESEARCH_MIN floors lambda at 0.1.  A
         trust region has no such floor, and it is not confined to the Newton
         direction.

             phi(u) = 1/2|R|^2,   R = f_int - f_ext = -r0
             p_N : J p_N = r0            (what PARDISO returned)
             d   = J^T r0 = -g           (steepest descent, formed as J^T)
             p_C = (|d|^2/|J d|^2) d
             p   = pa*d + pb*p_N         (Newton, Cauchy, or the blend)
             J p = pa*w + pb*r0          (so the model needs no matrix)
             rho = [phi(u)-phi(u+p)] / [phi(u) - 1/2|R+Jp|^2]

         Transactionality follows the accepted BT block: dam, damvisc and
         xstate are snapshotted once and restored before EVERY trial, and the
         loop always ends by re-evaluating the point it is going to keep, so
         what survives is that point own state and not a residue of the last
         trial.  cam, qa and uam are restored as well - cam[0] is a running
         MAXIMUM (resultsini.c:77-79) and uam[0]=max(uam[0],cam[0]) is never
         reset inside a step (nonlingeo.c) - so without this the rejected
         trials would permanently inflate the displacement criterion.  The
         kept point re-evaluation is compared against the value the trial
         pass measured for it: they must agree, and a mismatch is reported
         as an IMPURITY instead of being absorbed.

         Convergence of the INCREMENT is untouched: checkconvergence still
         decides on ram/cam/qa/qam from the unmodified residual. */

      damage_dl_used=0;
      if((damage_dl_on==1)&&(damage_dl_have==1)&&(damage_de12_enabled)&&
         (damage_ct_on==0)&&(damage_ct_used==0)&&
         (idamagereeq==0)&&(ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&
         (*ithermal<2)&&(*idrct==0)&&(*mortar<=1)){
        ITG tnst,tii,tjj,tacc,tkind,tbnd,tkkind;
        double tpa,tpb,tnrm,tphi,tl2,tinf,tpred,tared,trho;
        double tpcn,tpnn,tdold;
        double tkpa,tkpb,tkphi,tknrm,tkrho;
        double tqas[4],tuams[2];
        const char *tname;

        tnst=*nstate_;
        if(damage_dl_res==NULL){
          NNEW(damage_dl_res,double,neq[1]);
          NNEW(damage_dl_dam,double,mi[0]**ne);
          NNEW(damage_dl_visc,double,mi[0]**ne);
          if(tnst>0) NNEW(damage_dl_xs,double,tnst*mi[0]**ne);
        }

        tpnn=sqrt(damage_dl_npn2);
        tpcn=damage_dl_tc*sqrt(damage_dl_nd2);
        damage_dl_phi0=0.5*damage_dl_nb2;
        tname="NEWTON";tkind=1;tii=0;
        if(damage_dl_delta<=0.){
          damage_dl_delta=damage_dl_d0fac*tpnn;
          damage_dl_dmax=1.e3*tpnn;
          if(damage_dl_delta<=0.){damage_dl_delta=1.;damage_dl_dmax=1.e3;}
          printf("[DAMAGE TR] ARMED inc=%" ITGFORMAT " iter=%" ITGFORMAT
                 " attempt=%" ITGFORMAT ": |p_N|=%.6e |p_C|=%.6e |R|2=%.6e "
                 "phi=%.6e Delta0=%.6e.  BK3 is bypassed for this attempt; "
                 "the correction is chosen by the trust region.%s",
                 iinc,iit,icutb+1,tpnn,tpcn,sqrt(damage_dl_nb2),
                 damage_dl_phi0,damage_dl_delta,"\n");
          fflush(stdout);
        }

        /* snapshot the trial-derived state and the convergence bookkeeping */
        isiz=mi[0]**ne;cpypardou(damage_dl_dam,dam,&isiz,&num_cpus);
        if(damage_damvisc!=NULL){
          isiz=mi[0]**ne;
          cpypardou(damage_dl_visc,damage_damvisc,&isiz,&num_cpus);
        }
        if((tnst>0)&&(damage_dl_xs!=NULL)){
          isiz=tnst*mi[0]**ne;
          cpypardou(damage_dl_xs,xstate,&isiz,&num_cpus);
        }
        /* uam is the running maximum correction over the WHOLE step and is
           updated only after this block, so the value captured here is the
           value BEFORE the rejected full Newton step.  Restoring it before
           the final evaluation is what keeps a rejected step out of the
           displacement criterion for the rest of the step. */
        for(tjj=0;tjj<4;tjj++) tqas[tjj]=qa[tjj];
        for(tjj=0;tjj<2;tjj++) tuams[tjj]=uam[tjj];

        tacc=-1;tkpa=0.;tkpb=1.;tkphi=0.;tknrm=tpnn;tkrho=0.;tkkind=1;
        for(tii=0;tii<damage_dl_maxtrial;tii++){
          if(damage_dl_neval>=damage_dl_maxeval) break;
          tdold=damage_dl_delta;

          /* ---- the dogleg step for the current radius.  Same function
             the geometry self-test exercised before the run started; a
             degenerate or non-finite input returns 0 and is a refusal, never
             a step of zero length. ---- */
          tkind=damage_dl_pick(damage_dl_delta,damage_dl_nd2,damage_dl_nw2,
                               damage_dl_npn2,damage_dl_dtpn,&tpa,&tpb,&tnrm);
          if(tkind==0){
            printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   " trial=%" ITGFORMAT ": the step construction REFUSED a "
                   "degenerate or non-finite input (Delta=%.6e |d|^2=%.6e "
                   "|Jd|^2=%.6e |p_N|^2=%.6e dot=%.6e).  Nothing is accepted; "
                   "the full Newton step is restored and the stock divergence "
                   "and cutback machinery takes over unchanged.%s",
                   iinc,iit,tii+1,damage_dl_delta,damage_dl_nd2,
                   damage_dl_nw2,damage_dl_npn2,damage_dl_dtpn,"\n");
            fflush(stdout);
            break;
          }
          tname=(tkind==1)?"NEWTON":((tkind==2)?"CAUCHY":"DOGLEG");
          tbnd=(tnrm>=0.99*damage_dl_delta)?1:0;
          tpred=damage_dl_pred(tpa,tpb,damage_dl_nb2,damage_dl_nd2,
                               damage_dl_nw2);
          if(tpred<=0.){
            printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   " trial=%" ITGFORMAT " %s Delta=%.6e |p|=%.6e pred=%.6e "
                   "<= 0 - the model itself promises nothing, shrinking "
                   "without evaluating%s",
                   iinc,iit,tii+1,tname,damage_dl_delta,tnrm,tpred,"\n");
            fflush(stdout);
            damage_dl_delta=0.25*((tnrm>0.)?tnrm:damage_dl_delta);
            if(damage_dl_delta<1.e-12*damage_dl_dmax) break;
            continue;
          }

          /* ---- transactional trial evaluation of the ORIGINAL residual --- */
          /* [DAMAGE TR] cam starts from the STOCK CLEAN state, never from a
             snapshot: cam[0] is a running maximum, and any snapshot taken
             after the full-Newton results() already carries that step's
             correction.  Same initialisation the iteration top uses. */
          for(tjj=0;tjj<3;tjj++) cam[tjj]=0.;
          for(tjj=3;tjj<5;tjj++) cam[tjj]=0.5;
          isiz=mi[0]**ne;cpypardou(dam,damage_dl_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_dl_visc,&isiz,&num_cpus);
          }
          if((tnst>0)&&(damage_dl_xs!=NULL)){
            isiz=tnst*mi[0]**ne;
            cpypardou(xstate,damage_dl_xs,&isiz,&num_cpus);
          }
          SFREE(v);SFREE(stx);SFREE(fn);
          for(tjj=0;tjj<neq[1];tjj++)
            b[tjj]=tpa*damage_dl_d[tjj]+tpb*damage_dl_pn[tjj];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_dl_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          damage_dl_neval++;
          tphi=0.;tinf=0.;
          for(tjj=0;tjj<neq[1];tjj++){
            tphi+=damage_dl_res[tjj]*damage_dl_res[tjj];
            if(fabs(damage_dl_res[tjj])>tinf) tinf=fabs(damage_dl_res[tjj]);
          }
          tl2=sqrt(tphi);tphi*=0.5;

          tared=damage_dl_phi0-tphi;
          trho=tared/tpred;
          if(trho>1.e-4){
            tacc=tii;tkpa=tpa;tkpb=tpb;tkphi=tphi;tknrm=tnrm;
            tkrho=trho;tkkind=tkind;
          }
          if(trho<0.25){
            damage_dl_delta=0.25*tnrm;
          }else if((trho>0.75)&&(tbnd==1)){
            damage_dl_delta=2.*tnrm;
            if(damage_dl_delta>damage_dl_dmax)
              damage_dl_delta=damage_dl_dmax;
          }
          printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT " trial=%"
                 ITGFORMAT " %s Delta=%.6e |p|=%.6e boundary=%" ITGFORMAT
                 " |R|2=%.6e |R|inf=%.6e phi=%.6e pred=%.6e ared=%.6e "
                 "rho=%.6e %s Delta %.6e -> %.6e%s",
                 iinc,iit,tii+1,tname,tdold,tnrm,tbnd,tl2,tinf,tphi,tpred,
                 tared,trho,(trho>1.e-4)?"ACCEPT":"REJECT",tdold,
                 damage_dl_delta,"\n");
          fflush(stdout);
          if(tacc>=0) break;
          if(damage_dl_delta<1.e-12*damage_dl_dmax){
            printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   ": the radius COLLAPSED to %.6e - no step of any length "
                   "along either leg reduces phi%s",
                   iinc,iit,damage_dl_delta,"\n");
            fflush(stdout);
            break;
          }
        }

        /* ---- leave the point that is kept, and prove the transaction ----
           cam/qa/uam go back first, so that only the kept step own
           evaluation contributes to the convergence bookkeeping. */
        for(tjj=0;tjj<4;tjj++) qa[tjj]=tqas[tjj];
        for(tjj=0;tjj<2;tjj++) uam[tjj]=tuams[tjj];
        if(tacc<0){
          tkpa=0.;tkpb=1.;tkkind=0;
          damage_dl_nfail++;
          damage_dl_nrej+=tii;
        }else{
          damage_dl_nacc++;
          damage_dl_nrej+=tacc;
          if(tkkind==1) damage_dl_nnewt++;
          else if(tkkind==2) damage_dl_ncau++;
          else damage_dl_ndog++;
        }
        tpa=tkpa;tpb=tkpb;
          /* [DAMAGE TR] cam starts from the STOCK CLEAN state, never from a
             snapshot: cam[0] is a running maximum, and any snapshot taken
             after the full-Newton results() already carries that step's
             correction.  Same initialisation the iteration top uses. */
          for(tjj=0;tjj<3;tjj++) cam[tjj]=0.;
          for(tjj=3;tjj<5;tjj++) cam[tjj]=0.5;
          isiz=mi[0]**ne;cpypardou(dam,damage_dl_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_dl_visc,&isiz,&num_cpus);
          }
          if((tnst>0)&&(damage_dl_xs!=NULL)){
            isiz=tnst*mi[0]**ne;
            cpypardou(xstate,damage_dl_xs,&isiz,&num_cpus);
          }
          SFREE(v);SFREE(stx);SFREE(fn);
          for(tjj=0;tjj<neq[1];tjj++)
            b[tjj]=tpa*damage_dl_d[tjj]+tpb*damage_dl_pn[tjj];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_dl_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          damage_dl_neval++;
          tphi=0.;tinf=0.;
          for(tjj=0;tjj<neq[1];tjj++){
            tphi+=damage_dl_res[tjj]*damage_dl_res[tjj];
            if(fabs(damage_dl_res[tjj])>tinf) tinf=fabs(damage_dl_res[tjj]);
          }
          tl2=sqrt(tphi);tphi*=0.5;

        printf("[DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT " KEEP %s "
               "|p|=%.6e phi_kept=%.6e phi_start=%.6e rho=%.6e "
               "Delta_next=%.6e evals=%" ITGFORMAT " armed_factorisations=%"
               ITGFORMAT " rollback_purity=%.3e%s",
               iinc,iit,
               (tkkind==0)?"FULL-NEWTON (nothing accepted; the stock "
               "divergence and cutback machinery takes over unchanged)":
               ((tkkind==1)?"NEWTON":((tkkind==2)?"CAUCHY":"DOGLEG")),
               tknrm,tphi,damage_dl_phi0,tkrho,damage_dl_delta,
               damage_dl_neval,damage_dl_nfact,
               ((tacc>=0)&&(tkphi>0.))?fabs(tphi-tkphi)/tkphi:0.,"\n");
        fflush(stdout);
        if((tacc>=0)&&(tkphi>0.)&&(fabs(tphi-tkphi)>1.e-12*tkphi)){
          printf("[DAMAGE TR] *WARNING IMPURITY: re-evaluating the kept step "
                 "gave phi=%.17e against %.17e in the trial pass.  The trial "
                 "is then NOT a pure function of the step, and every rho "
                 "above compared different physical states.%s",
                 tphi,tkphi,"\n");
          fflush(stdout);
        }
        if(cam[0]>tknrm*(1.+1.e-9)+1.e-30){
          printf("*ERROR [DAMAGE TR] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                 ": cam[0]=%.12e exceeds |p_kept|_2=%.12e.  cam[0] is the "
                 "largest component of the correction actually applied and "
                 "cannot exceed its 2-norm, so the state left behind does "
                 "NOT belong to the step that was kept and every "
                 "displacement criterion from here on is meaningless.  "
                 "Stopping the experiment.%s",iinc,iit,cam[0],tknrm,"\n");
          fflush(stdout);FORTRAN(stop,());
        }
        damage_dl_used=1;

        if((damage_dl_neval>=damage_dl_maxeval)||
           (damage_dl_nfact>=damage_dl_maxfact)){
          printf("[DAMAGE TR] BUDGET EXHAUSTED (evaluations %" ITGFORMAT
                 "/%" ITGFORMAT ", armed factorisations %" ITGFORMAT "/%"
                 ITGFORMAT ").  The trust region switches OFF; the rest of "
                 "this attempt runs the stock path and the wall, if it "
                 "returns, goes to the original stock stop.%s",
                 damage_dl_neval,damage_dl_maxeval,damage_dl_nfact,
                 damage_dl_maxfact,"\n");
          fflush(stdout);
          damage_dl_on=0;
        }
      }

      damage_linesearch_applied=0;
      if((damage_linesearch_mode==1)&&(damage_de12_enabled)&&
         (damage_dl_used==0)&&(damage_dl_on==0)&&(damage_ct_on==0)&&
         (damage_ct_used==0)&&
         (idamagereeq==0)&&(ncont==0)&&(*nmethod!=4)&&(*nmethod!=5)&&
         (*ithermal<2)&&(*idrct==0)){
        NNEW(res,double,neq[1]);
        calcresidual(nmethod,neq,res,fext,f,iexpl,nactdof,aux2,vold,vini,
                     &dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,
                     fini,islavnode,nslavnode,mortar,ntie,mi,nzs,&nasym,
                     &idamping,veold,adc,auc,cvini,cv,&alpham,&num_cpus);

        damage_linesearch_oldnorm=0.;
        damage_linesearch_fullnorm=0.;
        for(i=0;i<neq[0];i++){
          if(fabs(resold[i])>damage_linesearch_oldnorm)
            damage_linesearch_oldnorm=fabs(resold[i]);
          if(fabs(res[i])>damage_linesearch_fullnorm)
            damage_linesearch_fullnorm=fabs(res[i]);
        }

        damage_linesearch_active=0;
        damage_linesearch_nsoft=0;
        damage_linesearch_maxdd=0.;
        if((isfinite(damage_linesearch_oldnorm))&&
           (isfinite(damage_linesearch_fullnorm))&&
           (damage_linesearch_fullnorm>
              DAMAGE_LINESEARCH_GROWTH*damage_linesearch_oldnorm)){
          damage_linesearch_active=damage_de12_trial_softening(
              dam,damdamageini,ipkon,lakon,ielmat,mi[2],ndmcon,dmcon,
              *ndmat_,*ntmat_,ne0,mi[0],&damage_linesearch_nsoft,
              &damage_linesearch_maxdd);
        }

        if(damage_linesearch_active){
          sum1=0.;sum2=0.;
          for(i=0;i<neq[0];i++){
            sum1+=b[i]*resold[i];
            sum2+=b[i]*res[i];
          }

          if((!isfinite(sum1))||(!isfinite(sum2))||
             (fabs(sum1-sum2)<1.e-30)){
            flinesearch=DAMAGE_LINESEARCH_FALLBACK;
          }else{
            flinesearch=sum1/(sum1-sum2);
            if((!isfinite(flinesearch))||(flinesearch<=0.))
              flinesearch=DAMAGE_LINESEARCH_FALLBACK;
          }
          if(flinesearch>DAMAGE_LINESEARCH_MAX)
            flinesearch=DAMAGE_LINESEARCH_MAX;
          if(flinesearch<damage_ls_min)
            flinesearch=damage_ls_min;

          NNEW(damage_linesearch_step,double,neq[1]);
          isiz=neq[1];cpypardou(damage_linesearch_step,b,&isiz,&num_cpus);
          damage_linesearch_contracted=0;
          lsladder_start(&damage_lsl,flinesearch,
                         damage_linesearch_oldnorm,damage_ls_min,0.5,
                         damage_ls_trials,damage_ls_legacy);

          /* ---- BACKTRACKING LADDER PROBE (CCX_DAMAGE_LS_PROBE) --------
             The line search backtracks 1.0 -> 0.5 -> 0.1 and stops: three
             trials with a floor of 0.1.  At the increments that fail it
             reports contracted=0, i.e. even a tenth of the step makes the
             residual worse - by 14000x at increment 346.  Two things can
             cause that and they need different answers:

               the DIRECTION is usable and the floor is simply too coarse
               for a damage transition, in which case some smaller alpha
               contracts;

               the direction is wrong, in which case no alpha does.

             This walks a finer ladder and PRINTS the residual at each
             alpha.  It decides nothing: the normal trial loop below runs
             afterwards, unchanged, and its own final evaluation is what
             the solver keeps.  alpha=1 is evaluated TWICE so the log
             carries its own proof that a trial evaluation is
             reproducible - dam and damvisc are rebuilt from the committed
             baseline on every results() call, and this is what shows it
             rather than assuming it. */

          if(damage_ls_probe!=0){
            static const double lsp[9]={1.,1.,0.5,0.2,0.1,0.03,0.01,
                                        0.003,0.001};
            ITG lpi;
            double lpr,lp0=-1.;
            for(lpi=0;lpi<9;lpi++){
              SFREE(v);SFREE(stx);SFREE(fn);
              for(i=0;i<neq[1];i++)
                b[i]=lsp[lpi]*damage_linesearch_step[i];
              MNEW(v,double,mt**nk);
              isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
              NNEW(stx,double,6*mi[0]**ne);
              MNEW(fn,double,mt**nk);
              if(ne1d2d==1)NNEW(inum,ITG,*nk);
              results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
                  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
                  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
                  prestr,iprestr,filab,eme,emn,een,iperturb,
                  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
                  ndirboun,xbounact,nboun,ipompc,
                  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
                  accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
                  plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,
                  matname,mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,
                  ikboun,ilboun,ener,enern,emeini,xstaten,eei,enerini,
                  cocon,ncocon,set,nset,istartset,iendset,ialset,nprint,
                  prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,nelemload,
                  nload,ikmpc,ilmpc,istep,&iinc,springarea,&reltime,&ne0,
                  thicke,shcon,nshcon,sideload,xloadact,xloadold,&icfd,
                  inomat,pslavsurf,pmastsurf,mortar,islavact,cdn,
                  islavnode,nslavnode,ntie,clearini,islavsurf,ielprop,
                  prop,energyini,energy,&kscale,iponoeln,inoeln,nener,
                  orname,network,ipobody,xbodyact,ibody,typeboun,itiefac,
                  tieset,smscale,&mscalmethod,nbody,t0g,t1g,islavquadel,
                  aut,irowt,jqt,&mortartrafoflag,&intscheme,physcon,dam,
                  damn,iponoel);
              if(ne1d2d==1)SFREE(inum);
              calcresidual(nmethod,neq,res,fext,f,iexpl,nactdof,aux2,vold,
                           vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                           fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                           nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,
                           &alpham,&num_cpus);
              lpr=0.;
              for(i=0;i<neq[0];i++) if(fabs(res[i])>lpr) lpr=fabs(res[i]);
              if(lpi==0) lp0=lpr;
              printf("[LS-PROBE] inc=%" ITGFORMAT " attempt=%" ITGFORMAT
                     " iter=%" ITGFORMAT " alpha=%.4f |R|inf=%.6e"
                     " (res_old=%.6e, ratio=%.3e)%s\n",
                     iinc,icutb+1,iit,lsp[lpi],lpr,
                     damage_linesearch_oldnorm,
                     lpr/(damage_linesearch_oldnorm+1.e-300),
                     ((lpi==1)&&(lp0>=0.))?
                     ((fabs(lpr-lp0)<=1.e-12*(lp0+1.e-300))?
                      "  [repeat of alpha=1: REPRODUCIBLE]":
                      "  [repeat of alpha=1: NOT REPRODUCIBLE]"):"");
            }
            fflush(stdout);
          }

          for(damage_linesearch_trial=1;
              damage_linesearch_trial<=damage_ls_trials;
              damage_linesearch_trial++){
            SFREE(v);SFREE(stx);SFREE(fn);
            for(i=0;i<neq[1];i++)
              b[i]=flinesearch*damage_linesearch_step[i];

            /* The bordered step is ONE Newton step in (u,lambda).
               Damping only its displacement half leaves the prescribed
               dofs at the undamped load factor, so the trial state is
               not on the ray the line search thinks it is contracting
               along - and the constraint it was solved to satisfy is not
               satisfied at any point of that ray.  Contract both. */

            if((pf_on==1)&&(pf_applied==1)){
              pf_lam=pf_lam0it+flinesearch*pf_dlamit;
              for(i=0;i<*nboun;i++)
                xbounact[i]=xbounold[i]+(xboun[i]-xbounold[i])*pf_lam;
            }

            MNEW(v,double,mt**nk);
            isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
            NNEW(stx,double,6*mi[0]**ne);
            MNEW(fn,double,mt**nk);

            if(ne1d2d==1)NNEW(inum,ITG,*nk);
            results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
                elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
                ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
                prestr,iprestr,filab,eme,emn,een,iperturb,
                f,fn,nactdof,&iout,qa,vold,b,nodeboun,
                ndirboun,xbounact,nboun,ipompc,
                nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
                &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
                xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
                &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
                emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
                iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
                fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
                &reltime,&ne0,thicke,shcon,nshcon,
                sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
                mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
                islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
                inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
                itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
                islavquadel,aut,irowt,jqt,&mortartrafoflag,
                &intscheme,physcon,dam,damn,iponoel);
            if(ne1d2d==1)SFREE(inum);

            calcresidual(nmethod,neq,res,fext,f,iexpl,nactdof,aux2,vold,
                         vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                         fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                         nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,
                         &alpham,&num_cpus);
            damage_linesearch_dampednorm=0.;
            for(i=0;i<neq[0];i++){
              if(fabs(res[i])>damage_linesearch_dampednorm)
                damage_linesearch_dampednorm=fabs(res[i]);
            }

            /* The ladder decides: accept, descend, or give up.  It is a
               separate unit (lsladder.c) with its own regression test,
               because both of the things it used to get wrong - the rungs
               it could reach and which rung it handed back - were the
               cause of the recorded wall. */

            damage_lsr=lsladder_step(&damage_lsl,
                                     damage_linesearch_dampednorm);
            if(damage_lsr==1){
              damage_linesearch_contracted=1;
              break;
            }
            if(damage_lsr<0) break;
            flinesearch=damage_lsl.alpha;
          }
          if(damage_linesearch_trial>damage_ls_trials)
            damage_linesearch_trial=damage_ls_trials;

          /* Nothing contracted.  Falling out with the LAST trial in b is
             the defect: the last trial is the floor, and the floor is not
             the best of what was measured.  Re-evaluate at the best alpha
             instead, so the step the solver takes is never worse than the
             best step this search actually saw. */

          if((damage_linesearch_contracted==0)&&
             (fabs(lsladder_final(&damage_lsl)-flinesearch)>1.e-12)){
            flinesearch=lsladder_final(&damage_lsl);
            SFREE(v);SFREE(stx);SFREE(fn);
            for(i=0;i<neq[1];i++)
              b[i]=flinesearch*damage_linesearch_step[i];
            MNEW(v,double,mt**nk);
            isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
            NNEW(stx,double,6*mi[0]**ne);
            MNEW(fn,double,mt**nk);
            if(ne1d2d==1)NNEW(inum,ITG,*nk);
            results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
                elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
                ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
                prestr,iprestr,filab,eme,emn,een,iperturb,
                f,fn,nactdof,&iout,qa,vold,b,nodeboun,
                ndirboun,xbounact,nboun,ipompc,
                nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
                &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
                xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
                &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
                emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
                iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
                fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
                &reltime,&ne0,thicke,shcon,nshcon,
                sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
                mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
                islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
                inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
                itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
                islavquadel,aut,irowt,jqt,&mortartrafoflag,
                &intscheme,physcon,dam,damn,iponoel);
            if(ne1d2d==1)SFREE(inum);
            calcresidual(nmethod,neq,res,fext,f,iexpl,nactdof,aux2,vold,
                         vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                         fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                         nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,
                         &alpham,&num_cpus);
            damage_linesearch_dampednorm=damage_lsl.bestres;
            damage_ls_bestused++;
          }
          SFREE(damage_linesearch_step);
          damage_linesearch_applied=1;
        }
        SFREE(res);
      }

      /* ---- [DAMAGE RAY] residual scan along the Newton direction --------
         Runs only where BK3 cannot: idamagereeq==1, the same-load solve after
         a topology change.  b still holds the Newton direction here; below,
         calcresidual overwrites b with the residual, so the direction is
         saved and restored explicitly.

         Two-point test first (alpha=0 and alpha=1) on every re-equilibration
         iteration - cheap, and alpha=0 is the reference this pass has no
         other way to obtain, resold being stale under this gate.  The full
         ray is walked only when the full step GREW the norm, which is the
         fatal signature, and at most damage_ray_max times.

         The walk ends on alpha=1, so v/stx/fn/dam/xstate are left in exactly
         the state the unprobed code would have produced. */
      damage_ray_incok=1;
      if(damage_ray_ninc>0){
        damage_ray_incok=0;
        for(i=0;i<damage_ray_ninc;i++)
          if(iinc==damage_ray_inc[i]) damage_ray_incok=1;
      }
      if((damage_ray_probe)&&(damage_ray_incok==1)&&
         (idamagereeq==1)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
        /* 0 and 1 first; then the ladder from 1/1024 up, ending on the
           repeated 0.25 (purity pair, indices 10 and 12) and on alpha=1 so
           the state left behind is the unprobed one. */
        static const double rayA[]={0.,1.,
                                    0.0009765625,0.001953125,0.00390625,
                                    0.0078125,0.015625,0.03125,0.0625,0.125,
                                    0.25,0.5,0.25,1.};
        ITG rnpt,rii,rjj,rfull,rflip,rnpl[14],rnuc[14];
        double ra,rinf[14],rl2[14],rptr[14],rn2[14],rs;

        if(damage_ray_p==NULL){
          NNEW(damage_ray_p,double,neq[1]);
          NNEW(damage_ray_res,double,neq[1]);
        }
        isiz=neq[1];cpypardou(damage_ray_p,b,&isiz,&num_cpus);

        rs=(qam[0]>0.)?qam[0]:1.;
        rfull=0; rflip=0;
        rnpt=2;                       /* alpha=0 and alpha=1 only, at first */
        for(rii=0;rii<rnpt;rii++){
          ra=rayA[rii];
          SFREE(v);SFREE(stx);SFREE(fn);
          for(rjj=0;rjj<neq[1];rjj++) b[rjj]=ra*damage_ray_p[rjj];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
              &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
              xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
              &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
              emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
              iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
              fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
              &reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
              mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
              islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
              inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
              itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
              islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ray_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,
                       nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,
                       &alpham,&num_cpus);
          rinf[rii]=0.; rl2[rii]=0.; rptr[rii]=0.;
          for(rjj=0;rjj<neq[0];rjj++){
            if(fabs(damage_ray_res[rjj])>rinf[rii])
              rinf[rii]=fabs(damage_ray_res[rjj]);
            rl2[rii]+=damage_ray_res[rjj]*damage_ray_res[rjj];
            rptr[rii]+=damage_ray_p[rjj]*damage_ray_res[rjj];
          }
          rn2[rii]=sqrt(rl2[rii]);
          rl2[rii]=rn2[rii]/rs;
          damage_ray_tally(xstate,xstateini,dam,damdamageini,damage_damvisc,
                           stx,ipkon,lakon,ne0,mi[0],*nstate_,
                           &rnpl[rii],&rnuc[rii]);

          /* alpha=0 is the reference: keep the WHOLE residual vector and the
             per-integration-point branch map.  Three scalars agreeing with the
             linear model prove nothing about the tangent VECTOR - that was an
             overreach in J-14 and it is corrected here by measuring the
             componentwise defect E(alpha)=R(alpha)-(1-alpha)R(0). */
          if(rii==0){
            if(damage_ray_r0==NULL) NNEW(damage_ray_r0,double,neq[1]);
            isiz=neq[1];cpypardou(damage_ray_r0,damage_ray_res,&isiz,&num_cpus);
            if(damage_ray_cat==NULL) NNEW(damage_ray_cat,ITG,mi[0]*ne0);
            damage_ray_census(damage_ray_cat,xstate,xstateini,dam,
                              damdamageini,damage_damvisc,stx,ipkon,lakon,
                              ne0,mi[0],*nstate_);
          }

          /* The fatal signature is decided on the RESIDUAL ALONE.  The
             budget must not enter this test: conflating "no overshoot" with
             "overshoot, but the walk budget is spent" makes the label lie,
             and it lied - the first build capped at 8 walks, all 8 were
             consumed by inc<=152, and inc=154 then printed flip=no for a
             step that had not been tested at all. */
          if((rii==1)&&(rinf[0]>0.)){
            if(rinf[1]>damage_ray_growth*rinf[0]) rflip=1;
            if((rflip)&&(damage_ray_shots<damage_ray_max)){
              rfull=1; rnpt=14; damage_ray_shots++;
            }
          }
        }

        printf("[DAMAGE RAY] inc=%" ITGFORMAT " iter=%" ITGFORMAT
               " time=%.12e dt=%.6e qam=%.6e"
               " R0=%.6e R1=%.6e ratio=%.3e flip=%s walked=%s%s",
               iinc,iit,theta**tper,dtime,qam[0],
               rinf[0],rinf[1],(rinf[0]>0.)?rinf[1]/rinf[0]:-1.,
               rflip?"YES":"no",rfull?"yes":"no",
               "\n");
        if(rfull){
          for(rii=0;rii<14;rii++){
            printf("   alpha=%.9f |R|inf=%.6e |R|2=%.6e phi=%.6e"
                   " |R|2/qam=%.6e pTR=%+.6e plastic_ip=%" ITGFORMAT
                   " ucomp_ip=%" ITGFORMAT "%s",
                   rayA[rii],rinf[rii],rn2[rii],0.5*rn2[rii]*rn2[rii],
                   rl2[rii],rptr[rii],rnpl[rii],rnuc[rii],
                   "\n");
          }
          /* Which measure would have accepted?  Armijo with c1=1e-4 on each
             of the three, against the alpha=0 reference. */
          for(rii=2;rii<12;rii++){
            printf("   ARMIJO alpha=%.9f  inf:%s  l2:%s  phi:%s%s",
                   rayA[rii],
                   (rinf[rii]<=(1.-1.e-4*rayA[rii])*rinf[0])?"PASS":"fail",
                   (rn2[rii] <=(1.-1.e-4*rayA[rii])*rn2[0]) ?"PASS":"fail",
                   (0.5*rn2[rii]*rn2[rii]<=
                    (1.-1.e-4*rayA[rii])*0.5*rn2[0]*rn2[0])?"PASS":"fail",
                   "\n");
          }
          /* componentwise defect against the linear model, and the branch
             census transitions.  This is the part that can actually tell a
             consistent tangent from a change of branch. */
          if(damage_ray_r0!=NULL){
            ITG dii,djj,dworst,dncat,dfirste,dfirstip,dfirsta,dfirstb;
            double de,dn2,di,r0n2,r0ni,dwv,dpred;
            r0n2=0.;r0ni=0.;
            for(djj=0;djj<neq[0];djj++){
              r0n2+=damage_ray_r0[djj]*damage_ray_r0[djj];
              if(fabs(damage_ray_r0[djj])>r0ni) r0ni=fabs(damage_ray_r0[djj]);
            }
            r0n2=sqrt(r0n2);
            if(r0n2<=0.) r0n2=1.e-300;
            if(r0ni<=0.) r0ni=1.e-300;
            for(dii=2;dii<14;dii++){
              SFREE(v);SFREE(stx);SFREE(fn);
              for(djj=0;djj<neq[1];djj++)
                b[djj]=rayA[dii]*damage_ray_p[djj];
              MNEW(v,double,mt**nk);
              isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
              NNEW(stx,double,6*mi[0]**ne);
              MNEW(fn,double,mt**nk);
              if(ne1d2d==1)NNEW(inum,ITG,*nk);
              results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
                  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
                  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
                  prestr,iprestr,filab,eme,emn,een,iperturb,
                  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
                  ndirboun,xbounact,nboun,ipompc,
                  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
                  accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
                  plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,matname,
                  mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,
                  ener,enern,emeini,xstaten,eei,enerini,cocon,ncocon,set,
                  nset,istartset,iendset,ialset,nprint,prlab,prset,qfx,qfn,
                  trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,ilmpc,istep,
                  &iinc,springarea,&reltime,&ne0,thicke,shcon,nshcon,
                  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,
                  pmastsurf,mortar,islavact,cdn,islavnode,nslavnode,ntie,
                  clearini,islavsurf,ielprop,prop,energyini,energy,&kscale,
                  iponoeln,inoeln,nener,orname,network,ipobody,xbodyact,
                  ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
                  t0g,t1g,islavquadel,aut,irowt,jqt,&mortartrafoflag,
                  &intscheme,physcon,dam,damn,iponoel);
              if(ne1d2d==1)SFREE(inum);
              calcresidual(nmethod,neq,damage_ray_res,fext,f,iexpl,nactdof,
                           aux2,vold,vini,&dtime,accold,nk,adb,aub,jq,irow,
                           nzl,alpha,fextini,fini,islavnode,nslavnode,mortar,
                           ntie,mi,nzs,&nasym,&idamping,veold,adc,auc,cvini,
                           cv,&alpham,&num_cpus);
              dn2=0.;di=0.;dworst=-1;dwv=0.;dpred=0.;
              for(djj=0;djj<neq[0];djj++){
                de=damage_ray_res[djj]-(1.-rayA[dii])*damage_ray_r0[djj];
                dn2+=de*de;
                if(fabs(de)>di){
                  di=fabs(de);dworst=djj;dwv=damage_ray_res[djj];
                  dpred=(1.-rayA[dii])*damage_ray_r0[djj];
                }
              }
              dn2=sqrt(dn2);
              dncat=damage_ray_census_diff(damage_ray_cat,xstate,xstateini,
                        dam,damdamageini,damage_damvisc,stx,ipkon,lakon,ne0,
                        mi[0],*nstate_,&dfirste,&dfirstip,&dfirsta,&dfirstb);
              printf("   DEFECT a=%.6f |E|2/|R0|2=%.6e |E|inf/|R0|inf=%.6e"
                     " worstdof=%" ITGFORMAT " R=%.6e lin=%.6e E=%.6e"
                     "  switched_ip=%" ITGFORMAT,
                     rayA[dii],dn2/r0n2,di/r0ni,dworst,dwv,dpred,dwv-dpred,
                     dncat);
              if(dncat>0)
                printf("  first: el=%" ITGFORMAT " ip=%" ITGFORMAT
                       " cat %" ITGFORMAT "->%" ITGFORMAT,
                       dfirste,dfirstip,dfirsta,dfirstb);
              printf("%s","\n");
            }
            fflush(stdout);
          }
          /* purity: index 10 and index 12 are the same alpha, twice */
          printf("   PURITY alpha=%.6f evaluated twice: |R|inf %.17e vs %.17e"
                 "  %s%s",rayA[10],rinf[10],rinf[12],
                 (rinf[10]==rinf[12])?"BITWISE IDENTICAL":"*** DIFFERS ***",
                 "\n");
        }
        fflush(stdout);

        /* restore the Newton direction; v/stx/fn already hold the alpha=1
           state because the ray ends there */
        isiz=neq[1];cpypardou(b,damage_ray_p,&isiz,&num_cpus);
      }

      /* ---- [DAMAGE BT] transactional backtracking in re-equilibration ----
         THIS IS A SOLVER CHANGE, not a diagnostic.  It decides which point
         the Newton iteration continues from.

         BK3 cannot serve here: it is locked out of idamagereeq by two
         independent gates, its floor is 0.10, and on failure it silently
         keeps the floor-length step whose residual is known to be worse.
         This block restores exactly and hands the increment back to the
         standard cutback when no step is acceptable.

         Contract, one probe at a time:
           - topology is whatever the current remastruct produced; no deletion
             scan and no active-set change happens inside the search;
           - dam / damvisc / xstate are restored from the snapshot taken
             before the search, so every probe starts from the same state;
           - a rejected probe leaves nothing behind: the loop always ends by
             re-evaluating the point it is going to keep, so the surviving
             state is that point's own and not a residue of the last trial;
           - if nothing is acceptable the FULL step is restored and the normal
             divergence / cutback machinery takes over unchanged.

         Acceptance is Armijo on the mechanical infinity norm:
             |R(alpha)|inf <= (1 - c1*alpha) * |R(0)|inf ,  c1 = 1e-4.
         Default off; unset the flag and not one array is allocated. */
      if(((damage_bt_mode==1)||(damage_rescue_bt_on==1))&&
         (idamagereeq==1)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
        static const double btA[]={0.,1.,0.5,0.25,0.125,0.0625,0.03125,
                                   0.015625};
        ITG bii,bjj,bacc,bnst,bevt;
        double ba,br,bref,br2;

        bnst=*nstate_;
        if(damage_ray_p==NULL){
          NNEW(damage_ray_p,double,neq[1]);
          NNEW(damage_ray_res,double,neq[1]);
        }
        if(damage_bt_dam==NULL){
          NNEW(damage_bt_dam,double,mi[0]**ne);
          NNEW(damage_bt_visc,double,mi[0]**ne);
          if(bnst>0) NNEW(damage_bt_xs,double,bnst*mi[0]**ne);
        }
        isiz=neq[1];cpypardou(damage_ray_p,b,&isiz,&num_cpus);
        isiz=mi[0]**ne;cpypardou(damage_bt_dam,dam,&isiz,&num_cpus);
        if(damage_damvisc!=NULL){
          isiz=mi[0]**ne;
          cpypardou(damage_bt_visc,damage_damvisc,&isiz,&num_cpus);
        }
        if((bnst>0)&&(damage_bt_xs!=NULL)){
          isiz=bnst*mi[0]**ne;cpypardou(damage_bt_xs,xstate,&isiz,&num_cpus);
        }

        bacc=-1; damage_bt_ntrial=0;
        /* non-monotone reference: the largest residual over the window.  The
           ring is reset at the first iteration of each same-load solve. */
        if(iit<=1){ damage_bt_nring=0; }
        for(bii=0;bii<8;bii++){
          ba=btA[bii];
          if((bii>=2)&&(ba<damage_bt_floor-1.e-12)) break;
          isiz=mi[0]**ne;cpypardou(dam,damage_bt_dam,&isiz,&num_cpus);
          if(damage_damvisc!=NULL){
            isiz=mi[0]**ne;
            cpypardou(damage_damvisc,damage_bt_visc,&isiz,&num_cpus);
          }
          if((bnst>0)&&(damage_bt_xs!=NULL)){
            isiz=bnst*mi[0]**ne;
            cpypardou(xstate,damage_bt_xs,&isiz,&num_cpus);
          }
          SFREE(v);SFREE(stx);SFREE(fn);
          for(bjj=0;bjj<neq[1];bjj++) b[bjj]=ba*damage_ray_p[bjj];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
              results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
                  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
                  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
                  prestr,iprestr,filab,eme,emn,een,iperturb,
                  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
                  ndirboun,xbounact,nboun,ipompc,
                  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
                  accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
                  plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,matname,
                  mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,
                  ener,enern,emeini,xstaten,eei,enerini,cocon,ncocon,set,
                  nset,istartset,iendset,ialset,nprint,prlab,prset,qfx,qfn,
                  trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,ilmpc,istep,
                  &iinc,springarea,&reltime,&ne0,thicke,shcon,nshcon,
                  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,
                  pmastsurf,mortar,islavact,cdn,islavnode,nslavnode,ntie,
                  clearini,islavsurf,ielprop,prop,energyini,energy,&kscale,
                  iponoeln,inoeln,nener,orname,network,ipobody,xbodyact,
                  ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
                  t0g,t1g,islavquadel,aut,irowt,jqt,&mortartrafoflag,
                  &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ray_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,nzs,
                       &nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          br=0.;br2=0.;
          for(bjj=0;bjj<neq[0];bjj++){
            if(fabs(damage_ray_res[bjj])>br) br=fabs(damage_ray_res[bjj]);
            br2+=damage_ray_res[bjj]*damage_ray_res[bjj];
          }
          br2=sqrt(br2);

          if((damage_evt_on==1)&&(bii==0)){
            if(damage_evt_sgn==NULL) NNEW(damage_evt_sgn,ITG,mi[0]*ne0);
            damage_evt_sign(stx,ipkon,lakon,ne0,mi[0],damage_evt_sgn);
            for(bjj=0;bjj<8;bjj++){
              damage_evt_nsw[bjj]=0;damage_evt_fe[bjj]=0;damage_evt_fp[bjj]=0;
            }
          }else if((damage_evt_on==1)&&(bii<8)){
            damage_evt_nsw[bii]=damage_evt_flips(stx,ipkon,lakon,ne0,mi[0],
                                                 damage_evt_sgn,
                                                 &damage_evt_fe[bii],
                                                 &damage_evt_fp[bii]);
          }
          if(bii==0){
            damage_bt_r0=br;
            bref=br;
            for(bjj=0;bjj<damage_bt_nring;bjj++)
              if(damage_bt_ring[bjj]>bref) bref=damage_bt_ring[bjj];
            continue;
          }
          if(bii==1){
            /* engagement gate: a full step that is not much worse than the
               current point is taken unchanged.  Newton is allowed to be
               temporarily worse; that is how the undamped run crosses. */
            if(br<=damage_bt_growth*damage_bt_r0){ bacc=1; break; }
          }
          damage_bt_ntrial++;
          printf("[DAMAGE BT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                 " trial=%" ITGFORMAT " alpha=%.6f R=%.6e R0=%.6e"
                 " Rref=%.6e R2=%.6e%s",
                 iinc,iit,damage_bt_ntrial,ba,br,damage_bt_r0,bref,br2,"\n");
          if(br<=(1.-1.e-4*ba)*bref){ bacc=bii; break; }
        }
        /* push the residual actually kept into the non-monotone ring */
        if(damage_bt_window>1){
          if(damage_bt_nring<damage_bt_window){
            damage_bt_ring[damage_bt_nring++]=damage_bt_r0;
          }else{
            for(bjj=1;bjj<damage_bt_window;bjj++)
              damage_bt_ring[bjj-1]=damage_bt_ring[bjj];
            damage_bt_ring[damage_bt_window-1]=damage_bt_r0;
          }
        }

        /* [DAMAGE EVT] Nothing satisfied Armijo.  Level one restores the
           full step here and hands the increment to the standard cutback.
           Level two first asks whether the obstruction is a UC6 crossing: if
           some ladder alpha changes the compression set, take the SMALLEST
           such alpha.  That leaves the facet definitely on its new branch, so
           the next assembly - built from vold, e_c3d_uc6.f:45 - carries
           ctan(1,1)=kn for it, which is the whole point.  If no alpha changes
           the set, nothing happens and level one's behaviour stands. */
        bevt=-1;
        if((bacc<0)&&(damage_evt_on==1)){
          for(bii=7;bii>=1;bii--){
            if(damage_evt_nsw[bii]>0){ bevt=bii; break; }
          }
          if(bevt>=0){
            bacc=bevt;
            damage_evt_nstep++;
            printf("[DAMAGE EVT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   " no alpha satisfied Armijo; the UC6 compression set first "
                   "changes at alpha=%.6f (%" ITGFORMAT " point(s) cross, "
                   "first el=%" ITGFORMAT " ip=%" ITGFORMAT
                   "); taking that EVENT STEP so the next tangent is built on "
                   "the new branch.  event steps so far: %" ITGFORMAT "%s",
                   iinc,iit,btA[bevt],damage_evt_nsw[bevt],
                   damage_evt_fe[bevt],damage_evt_fp[bevt],
                   damage_evt_nstep,"\n");
          }else{
            printf("[DAMAGE EVT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   " no alpha satisfied Armijo and NO ladder alpha changes the "
                   "UC6 compression set; the obstruction is not a crossing, so "
                   "level one behaviour stands (full step restored)%s",
                   iinc,iit,"\n");
          }
          fflush(stdout);
        }

        ba=(bacc>=0)?btA[bacc]:1.;
        isiz=mi[0]**ne;cpypardou(dam,damage_bt_dam,&isiz,&num_cpus);
        if(damage_damvisc!=NULL){
          isiz=mi[0]**ne;
          cpypardou(damage_damvisc,damage_bt_visc,&isiz,&num_cpus);
        }
        if((bnst>0)&&(damage_bt_xs!=NULL)){
          isiz=bnst*mi[0]**ne;cpypardou(xstate,damage_bt_xs,&isiz,&num_cpus);
        }
        SFREE(v);SFREE(stx);SFREE(fn);
        for(bjj=0;bjj<neq[1];bjj++) b[bjj]=ba*damage_ray_p[bjj];
        MNEW(v,double,mt**nk);
        isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
        NNEW(stx,double,6*mi[0]**ne);
        MNEW(fn,double,mt**nk);
        if(ne1d2d==1)NNEW(inum,ITG,*nk);
              results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
                  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
                  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
                  prestr,iprestr,filab,eme,emn,een,iperturb,
                  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
                  ndirboun,xbounact,nboun,ipompc,
                  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
                  accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
                  plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,matname,
                  mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,
                  ener,enern,emeini,xstaten,eei,enerini,cocon,ncocon,set,
                  nset,istartset,iendset,ialset,nprint,prlab,prset,qfx,qfn,
                  trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,ilmpc,istep,
                  &iinc,springarea,&reltime,&ne0,thicke,shcon,nshcon,
                  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,
                  pmastsurf,mortar,islavact,cdn,islavnode,nslavnode,ntie,
                  clearini,islavsurf,ielprop,prop,energyini,energy,&kscale,
                  iponoeln,inoeln,nener,orname,network,ipobody,xbodyact,
                  ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
                  t0g,t1g,islavquadel,aut,irowt,jqt,&mortartrafoflag,
                  &intscheme,physcon,dam,damn,iponoel);
        if(ne1d2d==1)SFREE(inum);
        printf("[DAMAGE BT] inc=%" ITGFORMAT " iter=%" ITGFORMAT
               " R0=%.6e trials=%" ITGFORMAT " %s%.6f%s",
               iinc,iit,damage_bt_r0,damage_bt_ntrial,
               (bevt>=0)?"EVENT-STEP alpha=":
               ((bacc>=0)?"ACCEPTED alpha=":
                          "NONE-ACCEPTED-full-step-restored alpha="),
               ba,"\n");
        fflush(stdout);
      }

      /* ---- [DAMAGE ABA] full-state purity test, fires once --------------- */
      damage_aba_hit=0;
      if((damage_aba_mode==1)&&(idamagereeq==1)&&(ncont==0)&&
         (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
        if(damage_aba_ninc==0){
          if(damage_aba_done==0) damage_aba_hit=1;
        }else{
          for(i=0;i<damage_aba_ninc;i++)
            if((damage_aba_inc[i]>0)&&(iinc==damage_aba_inc[i])){
              damage_aba_hit=1; damage_aba_inc[i]=-damage_aba_inc[i];
            }
        }
      }
      if(damage_aba_hit==1){
        ITG apass,abad=0,anst,nstx,nxst,neme;
        double aal;
        anst=*nstate_;
        nstx=6*mi[0]**ne;
        nxst=27*mi[0]**ne;
        neme=6*mi[0]**ne;
        if(damage_ray_p==NULL){
          NNEW(damage_ray_p,double,neq[1]);
          NNEW(damage_ray_res,double,neq[1]);
        }
        isiz=neq[1];cpypardou(damage_ray_p,b,&isiz,&num_cpus);
        NNEW(daba_res,double,neq[1]);
        NNEW(daba_v,double,mt**nk);
        NNEW(daba_stx,double,nstx);
        NNEW(daba_fn,double,mt**nk);
        NNEW(daba_f,double,neq[1]);
        NNEW(daba_dam,double,mi[0]**ne);
        NNEW(daba_visc,double,mi[0]**ne);
        if(anst>0) NNEW(daba_xs,double,anst*mi[0]**ne);
        NNEW(daba_eme,double,neme);
        NNEW(daba_stiff,double,nxst);
        NNEW(daba_qa,double,4);
        NNEW(daba_cam,double,5);

        for(apass=0;apass<3;apass++){
          aal=(apass==1)?(1.-damage_aba_a):damage_aba_a;
          SFREE(v);SFREE(stx);SFREE(fn);
          for(i=0;i<neq[1];i++) b[i]=aal*damage_ray_p[i];
          MNEW(v,double,mt**nk);
          isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
          NNEW(stx,double,6*mi[0]**ne);
          MNEW(fn,double,mt**nk);
          if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
              accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
              plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,matname,
              mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,
              ener,enern,emeini,xstaten,eei,enerini,cocon,ncocon,set,
              nset,istartset,iendset,ialset,nprint,prlab,prset,qfx,qfn,
              trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,ilmpc,istep,
              &iinc,springarea,&reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,
              pmastsurf,mortar,islavact,cdn,islavnode,nslavnode,ntie,
              clearini,islavsurf,ielprop,prop,energyini,energy,&kscale,
              iponoeln,inoeln,nener,orname,network,ipobody,xbodyact,
              ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
              t0g,t1g,islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
          if(ne1d2d==1)SFREE(inum);
          calcresidual(nmethod,neq,damage_ray_res,fext,f,iexpl,nactdof,aux2,
                       vold,vini,&dtime,accold,nk,adb,aub,jq,irow,nzl,alpha,
                       fextini,fini,islavnode,nslavnode,mortar,ntie,mi,nzs,
                       &nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
                       &num_cpus);
          if(apass==0){
            isiz=neq[1];cpypardou(daba_res,damage_ray_res,&isiz,&num_cpus);
            isiz=mt**nk;cpypardou(daba_v,v,&isiz,&num_cpus);
            isiz=nstx;cpypardou(daba_stx,stx,&isiz,&num_cpus);
            isiz=mt**nk;cpypardou(daba_fn,fn,&isiz,&num_cpus);
            isiz=neq[1];cpypardou(daba_f,f,&isiz,&num_cpus);
            isiz=mi[0]**ne;cpypardou(daba_dam,dam,&isiz,&num_cpus);
            if(damage_damvisc!=NULL){
              isiz=mi[0]**ne;
              cpypardou(daba_visc,damage_damvisc,&isiz,&num_cpus);
            }
            if(anst>0){
              isiz=anst*mi[0]**ne;cpypardou(daba_xs,xstate,&isiz,&num_cpus);
            }
            isiz=neme;cpypardou(daba_eme,eme,&isiz,&num_cpus);
            isiz=nxst;cpypardou(daba_stiff,xstiff,&isiz,&num_cpus);
            for(i=0;i<4;i++) daba_qa[i]=qa[i];
            for(i=0;i<5;i++) daba_cam[i]=cam[i];
          }else if(apass==2){
            printf("[DAMAGE ABA] inc=%" ITGFORMAT " iter=%" ITGFORMAT
                   " A=%.6f B=%.6f  full-state comparison of the two A"
                   " evaluations:%s",iinc,iit,damage_aba_a,
                   1.-damage_aba_a,"\n");
            damage_aba_cmp("residual",daba_res,damage_ray_res,neq[1],&abad);
            damage_aba_cmp("v",daba_v,v,mt**nk,&abad);
            damage_aba_cmp("stx",daba_stx,stx,nstx,&abad);
            damage_aba_cmp("fn",daba_fn,fn,mt**nk,&abad);
            damage_aba_cmp("f",daba_f,f,neq[1],&abad);
            damage_aba_cmp("dam",daba_dam,dam,mi[0]**ne,&abad);
            damage_aba_cmp("damvisc",daba_visc,damage_damvisc,
                           mi[0]**ne,&abad);
            damage_aba_cmp("xstate",daba_xs,xstate,anst*mi[0]**ne,&abad);
            damage_aba_cmp("eme",daba_eme,eme,neme,&abad);
            damage_aba_cmp("xstiff",daba_stiff,xstiff,nxst,&abad);
            damage_aba_cmp("qa",daba_qa,qa,4,&abad);
            damage_aba_cmp("cam",daba_cam,cam,5,&abad);
            printf("[DAMAGE ABA] VERDICT: %" ITGFORMAT
                   " array(s) differ -> %s%s",abad,
                   (abad==0)?"TRANSACTION IS CLEAN":
                             "STATE LEAK - a line search on this path compares"
                             " different physical states","\n");
            fflush(stdout);
          }
        }
        SFREE(daba_res);SFREE(daba_v);SFREE(daba_stx);SFREE(daba_fn);
        SFREE(daba_f);SFREE(daba_dam);SFREE(daba_visc);
        if(anst>0) SFREE(daba_xs);
        SFREE(daba_eme);SFREE(daba_stiff);SFREE(daba_qa);SFREE(daba_cam);
        damage_aba_done=1;

        /* leave the full step, exactly as the unprobed code would have */
        SFREE(v);SFREE(stx);SFREE(fn);
        isiz=neq[1];cpypardou(b,damage_ray_p,&isiz,&num_cpus);
        MNEW(v,double,mt**nk);
        isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
        NNEW(stx,double,6*mi[0]**ne);
        MNEW(fn,double,mt**nk);
        if(ne1d2d==1)NNEW(inum,ITG,*nk);
          results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
              elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
              ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
              prestr,iprestr,filab,eme,emn,een,iperturb,
              f,fn,nactdof,&iout,qa,vold,b,nodeboun,
              ndirboun,xbounact,nboun,ipompc,
              nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,
              accold,&bet,&gam,&dtime,&time,ttime,plicon,nplicon,
              plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,epn,matname,
              mi,&ielas,&icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,
              ener,enern,emeini,xstaten,eei,enerini,cocon,ncocon,set,
              nset,istartset,iendset,ialset,nprint,prlab,prset,qfx,qfn,
              trab,inotr,ntrans,fmpc,nelemload,nload,ikmpc,ilmpc,istep,
              &iinc,springarea,&reltime,&ne0,thicke,shcon,nshcon,
              sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,
              pmastsurf,mortar,islavact,cdn,islavnode,nslavnode,ntie,
              clearini,islavsurf,ielprop,prop,energyini,energy,&kscale,
              iponoeln,inoeln,nener,orname,network,ipobody,xbodyact,
              ibody,typeboun,itiefac,tieset,smscale,&mscalmethod,nbody,
              t0g,t1g,islavquadel,aut,irowt,jqt,&mortartrafoflag,
              &intscheme,physcon,dam,damn,iponoel);
        if(ne1d2d==1)SFREE(inum);
      }

      /* calculating the residual */

      // next line: change on 19072022
      if((*iexpl<=1)||(*nener==1)){
	calcresidual(nmethod,neq,b,fext,f,iexpl,nactdof,aux2,vold,vini,&dtime,
		     accold,nk,adb,aub,jq,irow,nzl,alpha,fextini,fini,
		     islavnode,nslavnode,mortar,ntie,mi,
		     nzs,&nasym,&idamping,veold,adc,auc,cvini,cv,&alpham,
		     &num_cpus);
      }

      if(damage_linesearch_applied){
        damage_linesearch_dampednorm=0.;
        for(i=0;i<neq[0];i++){
          if(fabs(b[i])>damage_linesearch_dampednorm)
            damage_linesearch_dampednorm=fabs(b[i]);
        }
        printf("[DAMAGE LINESEARCH BK3] inc=%" ITGFORMAT
               " attempt=%" ITGFORMAT " iter=%" ITGFORMAT
               " soft_elements=%" ITGFORMAT " max_dD=%.6e "
               "res_old=%.6e res_full=%.6e lambda=%.6f "
               "trials=%" ITGFORMAT " contracted=%" ITGFORMAT
               " res_damped=%.6e\n",iinc,icutb+1,iit,
               damage_linesearch_nsoft,damage_linesearch_maxdd,
               damage_linesearch_oldnorm,damage_linesearch_fullnorm,
               flinesearch,damage_linesearch_trial,
               damage_linesearch_contracted,damage_linesearch_dampednorm);
        fflush(stdout);
      }

      /* fix residuals for mortar contact, add contact forces */
      
      if(*mortar>1){
	for(k=0;k<neq[1];k++){
	  b[k]=b[k]-f_cs[k]-f_cm[k];}
      }	 

      isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);
      if(*ithermal!=2){
	// next line: change on 19072022
	//if((*ithermal!=2)&&((*iexpl<=1)||(*nener==1))){
	for(k=0;k<6*mi[0]*ne0;++k){
	  sti[k]=stx[k];
	}
      
	/* calculating the ratio of the smallest to largest pressure
	   for face-to-face contact
	   only done at the end of a step */

	if((*mortar==1)&&(1.-theta-dtheta<=1.e-6)){
	  FORTRAN(negativepressure,(&ne0,ne,mi,stx,&pressureratio));
	}else{pressureratio=0.;}
      }

      SFREE(v);SFREE(stx);SFREE(fn);
      
      if((idamping==1)&&(*iexpl<=1)){SFREE(adc);SFREE(auc);}

      if(*iexpl<=1){
	  
	/* store the residual forces for the next iteration */

	if(*ithermal!=2){
	  if(cam[0]>uam[0]){
	    uam[0]=cam[0];}
	  if(qau<1.e-10){
	    if(qa[0]>ea*qam[0]){
	      qam[0]=(qamold[0]*jnz+qa[0])/(jnz+1);}
	    else {
	      qam[0]=qamold[0];}
	  }
	  /* CCX_DAMAGE_QAM_FLOOR - DIAGNOSTIC, default off.

	     The force criterion is RELATIVE: checkconvergence.c tests
	     ram[0] against c1[0]*qam[0], and qam[0] is a running average of
	     the internal force over increments.  On a specimen that is
	     unloading as it breaks, qam[0] follows the load down and the
	     absolute tolerance collapses with it.  Measured on
	     run_s0_fine_grad14 at its wall: qa=0.001552, qam=0.001604,
	     residual 0.005474 - a miss by 3.4x - while the specimen was
	     still carrying 2.33 N, 22.7% of its peak, and the residual had
	     just fallen 14x in one iteration.  The run died on a vanishing
	     reference, not on a growing residual.

	     Flooring qam at a fraction of the largest value it ever reached
	     keeps the tolerance tied to the load the specimen ONCE carried.

	     THIS CHANGES THE CONVERGENCE CRITERION AND THEREFORE THE ANSWER.
	     It is a diagnostic for whether a given wall is criterial or
	     physical.  A run that goes further with it is NOT thereby a
	     success - that has to be shown on the physics (E-108) - and
	     adopting it needs the full verify + ladder gate. */
	  if(damage_qam_floor>0.){
	    if(qam[0]>damage_qam_peak){damage_qam_peak=qam[0];}
	    if(qam[0]<damage_qam_floor*damage_qam_peak){
	      qam[0]=damage_qam_floor*damage_qam_peak;}
	  }
	}
	if(*ithermal>1){
	  if(cam[1]>uam[1]){
	    uam[1]=cam[1];}      
	  if(qau<1.e-10){
	    if(qa[1]>ea*qam[1]){
	      qam[1]=(qamold[1]*jnz+qa[1])/(jnz+1);}
	    else {
	      qam[1]=qamold[1];}
	  }
	}
      
	/* calculating the maximum residual */

	for(k=0;k<2;++k){
	  ram2[k]=ram1[k];
	  ram1[k]=ram[k];
	  ram[k]=0.;
	}
	if(*ithermal!=2){
	  damage_spc_fmax=0.;damage_spc_fnode=0;damage_spc_fcount=0;
	  for(k=0;k<neq[0];++k){
	    err=fabs(b[k]);
	    if((damage_spc_force!=0)&&(damage_spc_mask!=NULL)&&
	       (nactdofinv!=NULL)){
	      ITG spcnd=nactdofinv[k]/mt;
	      if((spcnd>=0)&&(spcnd<damage_spc_nk)&&
	         (damage_spc_mask[spcnd]!=0)){
	        damage_spc_fcount++;
	        if(err>damage_spc_fmax){
	          damage_spc_fmax=err;damage_spc_fnode=spcnd+1;}
	        continue;
	      }
	    }
	    if(err>ram[0]){
	      ram[0]=err;
	      ram[2]=k+0.5;}
	  }
	}
	if(*ithermal>1){
	  for(k=neq[0];k<neq[1];++k){
	    err=fabs(b[k]);
	    if(err>ram[1]){
	      ram[1]=err;
	      ram[3]=k+0.5;}
	  }
	}
	  
	/*   Divergence criteria for face-to-face penalty is different */
	  
	if(*mortar==1){
	  for(k=4;k<6;++k){
	    ram2[k]=ram1[k];
	    ram1[k]=ram[k];
	  } 
	  ram[4]=ram[0]+ram1[0];
	  ram[5]=(*ne-ne0)-(neold-ne0)+0.5;
	}
	  
	/* next line is inserted to cope with stress-less
	   temperature calculations */
	  
	if(*ithermal!=2){
	  if(ram[0]<1.e-6){
	    ram[0]=0.;} 
	  printf(" average force= %f\n",qa[0]);
	  printf(" time avg. forc= %f\n",qam[0]);
	  if((damage_spc_force!=0)&&(damage_spc_fcount>0)){
	    /* Auditable by construction: the excluded peak is printed next to
	       the criterion it was excluded from, so a masked residual that
	       starts to grow is visible in the same place ram[0] is read. */
	    printf("[DAMAGE AUTOSPC-FORCE] excluded %" ITGFORMAT " dof(s) on "
	           "AUTOSPC-masked nodes from ram[0]; largest excluded "
	           "residual %.6e at node %" ITGFORMAT " (ram[0]=%.6e, "
	           "tolerance %.6e = %.4f x qam)%s",
	           damage_spc_fcount,damage_spc_fmax,damage_spc_fnode,
	           ram[0],ctrl[18]*qam[0],
	           (qam[0]>0.)?damage_spc_fmax/qam[0]:0.,"\n");
	    /* [CONVSTATE] the EXCLUDED node, read as a length.

	       Without this the probe below measures the peak of ram[0], and
	       ram[0] is what is left AFTER this exclusion - so the one node
	       the exclusion exists for was the one node never measured.  Cost
	       one s3rad run to notice: at the stall it reported node 8305 at
	       need_du/grip = 2.3e-05 while node 1246, the subject, was not in
	       the norm at all. */
	    if((damage_dstate.nk>0)&&(damage_spc_fnode>=1)&&
	       (damage_spc_fnode<=damage_dstate.nk)){
	      convstate_report(iinc,damage_spc_fnode,damage_spc_fmax,
	                       damage_dstate.diag[damage_spc_fnode-1],
	                       convstate_griplength(vold,mt,damage_lp.nodesb,
	                                            damage_lp.nb,
	                                            damage_lp.idir),
	                       ctrl[18]*qam[0]);
	    }
	  }
	  if((ITG)((double)nactdofinv[(ITG)ram[2]]/mt)+1==0){
	    printf(" largest residual force= %f\n",
		   ram[0]);
	  }else{
	    inode=(ITG)((double)nactdofinv[(ITG)ram[2]]/mt)+1;
	    idir=nactdofinv[(ITG)ram[2]]-mt*(inode-1);
	    printf(" largest residual force= %f in node %" ITGFORMAT
		   " and dof %" ITGFORMAT "\n",
		   ram[0],inode,idir);
	    /* [CONVSTATE] the same peak, read as a LENGTH.  Measurement only:
	       it decides nothing and changes no criterion - see convstate.c
	       for why the measurement has to come before the judgement. */
	    if((damage_dstate.nk>0)&&(inode>=1)&&(inode<=damage_dstate.nk)){
	      convstate_report(iinc,inode,ram[0],
	                       damage_dstate.diag[inode-1],
	                       convstate_griplength(vold,mt,damage_lp.nodesb,
	                                            damage_lp.nb,
	                                            damage_lp.idir),
	                       ctrl[18]*qam[0]);
	    }
	  }
	  printf(" largest increment of disp= %e\n",uam[0]);
	  if((ITG)cam[3]==0){
	    printf(" largest correction to disp= %e\n\n",
		   cam[0]);
	  }else{
	    inode=(ITG)((double)nactdofinv[(ITG)cam[3]]/mt)+1;
	    idir=nactdofinv[(ITG)cam[3]]-mt*(inode-1);
	    printf(" largest correction to disp= %e in node %" ITGFORMAT
		   " and dof %" ITGFORMAT "\n\n",cam[0],inode,idir);
	  }
	}
	if(*ithermal>1){
	  if(ram[1]<1.e-6){
	    ram[1]=0.;}      
	  printf(" average flux= %f\n",qa[1]);
	  printf(" time avg. flux= %f\n",qam[1]);
	  if((ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1==0){
	    printf(" largest residual flux= %f\n",
		   ram[1]);
	  }else{
	    inode=(ITG)((double)nactdofinv[(ITG)ram[3]]/mt)+1;
	    idir=nactdofinv[(ITG)ram[3]]-mt*(inode-1);
	    printf(" largest residual flux= %f in node %" ITGFORMAT
		   " and dof %" ITGFORMAT "\n",ram[1],inode,idir);
	  }
	  printf(" largest increment of temp= %e\n",uam[1]);
	  if((ITG)cam[4]==0){
	    printf(" largest correction to temp= %e\n\n",
		   cam[1]);
	  }else{
	    inode=(ITG)((double)nactdofinv[(ITG)cam[4]]/mt)+1;
	    idir=nactdofinv[(ITG)cam[4]]-mt*(inode-1);
	    printf(" largest correction to temp= %e in node %" ITGFORMAT
		   " and dof %" ITGFORMAT "\n\n",cam[1],inode,idir);
	  }
	}
	fflush(stdout);
	  
	FORTRAN(writecvg,(istep,&iinc,&icutb,&iit,ne,&ne0,ram,qam,cam,uam,
			  ithermal));

        /* NC1 damage-aware slow-convergence gate.

           Keep checkconvergence() itself stock.  ctrl[3] (ic) is raised only
           for this call when progressive damage is actually changing, or
           during the same-load equilibrium solve of a DE1.3 terminal
           deletion, and both the residual and correction have contracted
           twice.  If a
           previously extended solve loses that trend after stock ic, setting
           ic to the current iteration makes the existing "iit==ic" branch
           perform the normal dc cutback instead of risking an unbounded run.
           ctrl[3] is restored immediately after the call. */

        ctrl[3]=icref;
        damage_slow_active=0;
        damage_slow_allow=0;
        damage_slow_nsoft=0;
        damage_slow_maxdd=0.;
        damage_slow_estres=DAMAGE_SLOW_NEWTON_INVALID_EST;
        damage_slow_estcorr=DAMAGE_SLOW_NEWTON_INVALID_EST;
        damage_slow_esttotal=DAMAGE_SLOW_NEWTON_INVALID_EST;
        damage_slow_rratio=0.;
        damage_slow_cratio=0.;

        if((damage_de12_enabled)&&(damdamageini!=NULL)&&
           (*iexpl<=1)&&(*nmethod!=4)&&(*nmethod!=5)&&
           (*ithermal<2)&&(*idrct==0)&&(ncont==0)&&
           (iit>=(ITG)irref)){
          /* During a terminal-deletion transaction the surviving elements
             need not acquire additional D: the topology change itself is the
             damage event being equilibrated.  Treat that same-load solve as
             damage-active, while retaining all monotonicity and cost gates
             below.  Legacy A3 transactions do not set
             damage_de13_transaction and therefore remain stock. */
          if((idamagereeq==1)&&(damage_de13_transaction==1)){
            damage_slow_active=1;
            damage_slow_nsoft=damage_tent_count;
          }else{
            damage_slow_active=damage_de12_trial_softening(
                dam,damdamageini,ipkon,lakon,ielmat,mi[2],ndmcon,dmcon,
                *ndmat_,*ntmat_,ne0,mi[0],&damage_slow_nsoft,
                &damage_slow_maxdd);
          }
        }

        if(damage_slow_active){
          damage_slow_allow=damage_slow_newton_allow(
              iit,ram,ram1,ram2,cam,uam,damage_slow_camprev1,
              damage_slow_camprev2,qa,qam,ctrl,damage_slow_maxiters,
              &damage_slow_estres,&damage_slow_estcorr,
              &damage_slow_esttotal,&damage_slow_rratio,
              &damage_slow_cratio);
        }

        if((damage_slow_allow)&&
           ((damage_slow_esttotal>(ITG)icref)||damage_slow_extended)){
          ctrl[3]=(double)damage_slow_maxiters+0.5;
          if((damage_slow_extended==0)||(iit==(ITG)icref)||
             ((iit%5)==0)){
            printf("[DAMAGE NEWTON EXTEND] inc=%" ITGFORMAT
                   " attempt=%" ITGFORMAT " iter=%" ITGFORMAT
                   " soft_elements=%" ITGFORMAT " max_dD=%.6e "
                   "est_res=%" ITGFORMAT " est_corr=%" ITGFORMAT
                   " est_total=%" ITGFORMAT " cap=%" ITGFORMAT
                   " r_ratio=%.6f c_ratio=%.6f\n",
                   iinc,icutb+1,iit,damage_slow_nsoft,damage_slow_maxdd,
                   damage_slow_estres,damage_slow_estcorr,
                   damage_slow_esttotal,damage_slow_maxiters,
                   damage_slow_rratio,damage_slow_cratio);
            fflush(stdout);
          }
          damage_slow_extended=1;
        }else if((damage_slow_extended)&&(iit>(ITG)icref)){
          /* Avoid checkconvergence's fatal iit>ic guard: ic==iit reaches its
             ordinary too-slow cutback path on this iteration. */
          ctrl[3]=(double)iit+0.5;
          printf("[DAMAGE NEWTON STOP] inc=%" ITGFORMAT
                 " attempt=%" ITGFORMAT " iter=%" ITGFORMAT
                 " reason=lost-contraction-or-cost-bound; stock cutback\n",
                 iinc,icutb+1,iit);
          fflush(stdout);
        }

        damage_slow_camprev2=damage_slow_camprev1;
        damage_slow_camprev1=cam[0];

        /* NC2: a terminal topology solve has zero external load increment.
           Stock therefore normalizes its correction by the tiny displacement
           caused by deletion alone.  This can demand many factorizations even
           after the force residual satisfies the stock tolerance.  For this
           nested solve only, retain the displacement scale of the already
           converged physical increment.  No tolerance is relaxed: the same
           stock relative correction and force-residual tests are applied to
           the combined physical-increment/topology correction. */
        /* checkconvergence advances theta by dtheta on acceptance and
           sizes the next increment from dthetaref, which is untouched,
           so overriding dtheta here lands theta exactly on the load
           factor the constraint produced without disturbing the stock
           step controller. */

        /* While a descent is walking, lambda is no longer tied to
           theta, so the dissipation controller's handover would
           overwrite the step the descent just restored - which is
           why the descent only ever got three attempts. */
        if((damage_diss_ctrl>=1)&&(damage_path_desc==0)&&(damage_arc==0)){
          dtheta=damage_diss_lamcur-theta;
          if(dtheta<*tmin) dtheta=*tmin;
        }

        damage_reeq_uam_actual[0]=uam[0];
        damage_reeq_uam_actual[1]=uam[1];
        if((damage_reeq_scale_mode==1)&&(idamagereeq==1)&&
           (damage_de13_transaction==1)){
          if(uam[0]<damage_reeq_uam_ref[0])
            uam[0]=damage_reeq_uam_ref[0];
          if(uam[1]<damage_reeq_uam_ref[1])
            uam[1]=damage_reeq_uam_ref[1];
          if(iit==1){
            printf("[DAMAGE REEQ NC2] inc=%" ITGFORMAT
                   " actual_uam=%.6e physical_ref=%.6e\n",
                   iinc,damage_reeq_uam_actual[0],damage_reeq_uam_ref[0]);
            fflush(stdout);
          }
        }

        /* [DAMAGE RESCUE] the step the failing attempt actually used is
           the last admissible one; checkconvergence is about to shrink it. */
        damage_rescue_dtheta_last=dtheta;
        damage_rescue_dthetaref_last=dthetaref;
        ccx_rescue_req=0;
        ccx_rescue_arm=0;
        if((damage_rescue_mode==1)&&(damage_rec_disarmed==0)&&
           (damage_rescue_used<damage_rescue_maxlevel)&&(ncont==0)&&
           (*nmethod!=4)&&(*nmethod!=5)&&(*ithermal<2)&&(*idrct==0)){
          ccx_rescue_arm=1;
          /* [DAMAGE CT] level 4 is granted only if the continuation could
             actually arm.  The candidate scan needs committed ring data
             only, so it is evaluated HERE, before checkconvergence defers
             the stop - a refusal then leaves the stock path byte-identical
             to the Rescue2+dogleg control. */
          if((damage_ct_mode==1)&&(damage_rescue_used>=3)&&
             (damage_ct_on==0)){
            ITG cse,csi,csn,csr,csok=0;
            double csm[3],csk,csd,cst;
            if((damage_ct_alloc==1)&&(damage_ct_nring>=6)&&
               (damage_ct_refused==0)){
              csok=damage_ct_select(damage_ct_ring,damage_ct_fl,damage_ct_dt,
                                    ipkon,lakon,ielprop,prop,ne0,mi[0],
                                    damage_ct_head,&cse,&csi,csm,&csk,&csd,
                                    &cst,&csn,&csr,damage_ct_kaptol,damage_ct_bl,damage_ct_nbl);
            }else{csn=0;csr=1;}
            if(csok==0){
              printf("[DAMAGE CT] level 4 NOT granted at inc=%" ITGFORMAT
                     " (%s; candidates=%" ITGFORMAT ", ring=%" ITGFORMAT
                     "/6).  No continuation attempt is taken, the stock path "
                     "is untouched and this wall goes to the ORIGINAL stock "
                     "stop.%s",iinc,
                     (damage_ct_alloc==0)?"ring not allocated":
                     ((damage_ct_nring<6)?"ring incomplete":
                      ((damage_ct_refused!=0)?"already refused once":
                       ((csr==2)?"kappa unstable over the five intervals":
                        "no candidate with five admissible intervals"))),
                     csn,damage_ct_nring,"\n");
              fflush(stdout);
              damage_ct_refused=1;
              ccx_rescue_arm=0;
            }
          }
        }

	damage_arc_theta0=theta;
	checkconvergence(co,nk,kon,ipkon,lakon,ne,stn,nmethod, 
			 kode,filab,een,t1act,&time,epn,ielmat,matname,enern, 
			 xstaten,nstate_,istep,&iinc,iperturb,ener,mi,output,
			 ithermal,qfn,&mode,&noddiam,trab,inotr,ntrans,orab,
			 ielorien,norien,description,sti,&icutb,&iit,&dtime,qa,
			 vold,qam,ram1,ram2,ram,cam,uam,&ntg,ttime,&icntrl,
			 &theta,&dtheta,veold,vini,idrct,tper,&istab,tmax, 
			 nactdof,b,tmin,ctrl,amta,namta,itpamp,inext,&dthetaref,
			 &itp,&jprint,jout,&uncoupled,t1,&iitterm,nelemload,
			 nload,nodeboun,nboun,itg,ndirboun,&deltmx,&iflagact,
			 set,nset,istartset,iendset,ialset,emn,thicke,jobnamec,
			 mortar,nmat,ielprop,prop,&ialeatoric,&kscale,
			 energy,&allwk,&energyref,&emax,&r_abs,&enetoll,
			 energyini,
			 &allwkini,&temax,&sizemaxinc,&ne0,&neini,&dampwk,
			 &dampwkini,energystartstep);

        /* [DAMAGE RESCUE] checkconvergence deferred the stop.  It left
           exactly as an ordinary cutback leaves, so the standard rollback at
           `if(icutb!=0)` further down restores the start of this increment
           from the baselines saved when it began - no new snapshot.  Undo
           only its step reduction and arm BT for this one attempt. */
        if(ccx_rescue_req==1){
          ccx_rescue_req=0;
          ccx_rescue_arm=0;
          if(damage_corr_on==1){
            /* Wall inside the corridor.  Counted and judged HERE, not
               after some later converged increment, so a chain of walls
               is bounded even if nothing ever converges again. */
            ITG cwclose=0;
            damage_corr_nint++;
            damage_corr_nwalltot++;
            if(damage_corr_nwalltot>damage_corr_maxwall){
              printf("[DAMAGE CORR] CIRCUIT BREAKER (walls): %" ITGFORMAT
                     " walls inside the corridor%s",
                     damage_corr_nwalltot,"\n");
              cwclose=1;
            }else if(damage_corr_trial==1){
              /* the probe failed - fall back to the proven lambda.  This
                 is not a repeat: lambda changes. */
              damage_corr_trial=0;
              damage_corr_lam=damage_corr_lamstable;
              damage_corr_nsince=0;
              printf("[DAMAGE CORR] the probe did not hold at inc=%"
                     ITGFORMAT "; back to the proven lambda=%.3e "
                     "(intervention %" ITGFORMAT ")%s",
                     iinc,damage_corr_lam,damage_corr_nint,"\n");
            }else{
              /* Wall on the HELD lambda.  Repeating it would recompute
                 the identical attempt, so escalate one ladder step -
                 bounded - or close. */
              ITG ei;double lnext=0.;
              damage_corr_nwallstab++;
              for(ei=0;ei<damage_reg_nlam;ei++){
                if(damage_reg_lam[ei]>damage_corr_lam*1.0000001){
                  lnext=damage_reg_lam[ei];break;
                }
              }
              if((damage_corr_nwallstab>damage_corr_maxesc)||(lnext<=0.)){
                printf("[DAMAGE CORR] CIRCUIT BREAKER (escalation): wall "
                       "on the held lambda=%.3e at inc=%" ITGFORMAT ", %"
                       ITGFORMAT " consecutive, ladder %s%s",
                       damage_corr_lam,iinc,damage_corr_nwallstab,
                       (lnext<=0.)?"exhausted":"still open","\n");
                cwclose=1;
              }else{
                printf("[DAMAGE CORR] wall on the held lambda at inc=%"
                       ITGFORMAT ": an identical retry is refused; lambda "
                       "escalated %.3e -> %.3e (escalation %" ITGFORMAT
                       " of %" ITGFORMAT ")%s",
                       iinc,damage_corr_lam,lnext,damage_corr_nwallstab,
                       damage_corr_maxesc,"\n");
                damage_corr_lam=lnext;
                damage_corr_lamstable=lnext;
                damage_corr_nsince=0;
              }
            }
            if(cwclose==1){
              printf("[DAMAGE CORR] corridor CLOSES: %" ITGFORMAT
                     " increments, %" ITGFORMAT " intervention(s), %"
                     ITGFORMAT " regularized factorisation(s), step time "
                     "gained %.6e.  The mechanism DISARMS; the next wall "
                     "goes to the original stock stop%s",
                     damage_corr_ninc,damage_corr_nint,damage_corr_nfact,
                     (theta-damage_corr_theta0)**tper,"\n");
              damage_corr_on=0;damage_corr_lam=0.;
              damage_rec_disarmed=1;
            }
            damage_reg_lambda=damage_corr_lam;
            damage_reg_on=(damage_corr_lam>0.)?1:0;
            /* Bank BEFORE zeroing.  Dropping this counts only the
               factorisations of attempts that succeeded, so the price
               of the corridor would be understated by exactly the cost
               of its failures. */
            damage_corr_nfact+=damage_reg_napply;
            damage_reg_napply=0;
            damage_rescue_used=0;
            damage_rescue_bt_on=0;
            damage_evt_on=0;
            damage_rec_used_in_inc=1;
            fflush(stdout);
          }else{
          /* Recovery accounting.  A rescue that fires before WINDOW
             untouched increments have gone by has not recovered
             anything; MAXUNREC of those in a row and the mechanism
             stops pretending and disarms. */
          if(damage_rescue_used==0){
            if(damage_rec_healthy<damage_rec_window){
              damage_rec_unrec++;
            }else{
              damage_rec_unrec=0;
            }
            damage_rec_healthy=0;
            if(damage_rec_unrec>damage_rec_maxunrec){
              damage_rec_disarmed=1;
              printf("[DAMAGE RESCUE] recovery window BLOWN: %" ITGFORMAT
                     " consecutive rescues without %" ITGFORMAT
                     " clean increments in between.  The solver is being "
                     "carried, not recovering, so the mechanism DISARMS "
                     "and this wall goes to the original stock stop%s",
                     damage_rec_unrec,damage_rec_window,"\n");
              fflush(stdout);
              ccx_rescue_arm=0;
              ccx_rescue_req=0;
            }
          }
          if(damage_rec_disarmed==1){
            damage_rescue_used=damage_rescue_maxlevel;
          }else{
          damage_rescue_used++;
          damage_rec_used_in_inc=1;
          damage_rescue_bt_on=1;
          /* A wall reached with idamagereeq==0 never enters a same-load
             solve, so BT and the event step - both gated on
             idamagereeq==1 - cannot act on it, and levels 1 and 2 would
             recompute the identical attempt.  Measured: s3rad inc=569
             carries four [DAMAGE RESCUE] lines and not one [DAMAGE BT],
             and both retries there failed identically. */
          if((idamagereeq==0)&&(damage_reg_nlam>0)&&
             (damage_rescue_used<2)){
            printf("[DAMAGE RESCUE] this wall has idamagereeq=0: no "
                   "same-load solve, so levels 1 and 2 cannot act on it; "
                   "skipping straight to the regularized level%s","\n");
            damage_rescue_used=2;
          }
          /* [DAMAGE CT] level 4: Rescue2 levels 1-2 and the dogleg have
             all failed on this wall.  The standard cutback rollback has
             already restored the increment-start state, so everything below
             is measured on the last COMMITTED state.  Any refusal leaves
             that state untouched and hands the wall to the stock stop. */
          if((damage_ct_mode==1)&&(damage_rescue_used>=4)&&
             (damage_ct_on==0)&&(damage_ct_refused==0)){
            ITG cse,csi,csn,csr;
            double csm[3],csk,csd,cst;
            if((damage_ct_alloc==1)&&(damage_ct_nring>=6)&&
               (damage_ct_select(damage_ct_ring,damage_ct_fl,damage_ct_dt,
                                 ipkon,lakon,ielprop,prop,ne0,mi[0],
                                 damage_ct_head,&cse,&csi,csm,&csk,&csd,
                                 &cst,&csn,&csr,damage_ct_kaptol,damage_ct_bl,damage_ct_nbl)==1)){
              damage_ct_arm=1;   /* already gated above; this only reports */
              printf("[DAMAGE CT] level 4 pre-check PASSED at inc=%" ITGFORMAT
                     ": %" ITGFORMAT " candidate(s) with five admissible "
                     "committed intervals; best element %" ITGFORMAT " ip %"
                     ITGFORMAT ", kappa=%.6e (max/min=%.4f), ds0=%.6e.  The "
                     "increment gets one continuation attempt.%s",
                     iinc,csn,cse+1,csi+1,csk,cst,csd,"\n");
            }else{
              printf("[DAMAGE CT] level 4 REFUSED at inc=%" ITGFORMAT
                     " (%s; candidates=%" ITGFORMAT ", ring=%" ITGFORMAT
                     "/6).  No continuation attempt is granted, so the stock "
                     "path is untouched and this wall goes to the ORIGINAL "
                     "stock stop.%s",iinc,
                     (damage_ct_alloc==0)?"ring not allocated":
                     ((damage_ct_nring<6)?"ring incomplete":
                      ((csr==2)?"kappa unstable over the five intervals":
                       "no candidate with five admissible intervals")),
                     (damage_ct_alloc==1)?csn:0,damage_ct_nring,"\n");
              damage_ct_refused=1;
              damage_rescue_used=damage_rescue_maxlevel;
              ccx_rescue_arm=0;
            }
            fflush(stdout);
          }
          damage_evt_on=(damage_rescue_used==2)?1:0;
          damage_reg_on=0;
          /* [DAMAGE TR] guard: without damage_reg_nlam>0 this block indexes
             damage_reg_lam[-1] and switches the diagonal shift on with a
             garbage lambda.  Unreachable while only RESCUE3/CORRIDOR could
             raise the level count; reachable the moment any other mechanism
             claims level 3, which the dogleg does. */
          if((damage_rescue_used>=3)&&(damage_reg_nlam>0)){
            damage_reg_level=damage_rescue_used-3;
            if(damage_reg_level>=damage_reg_nlam)
              damage_reg_level=damage_reg_nlam-1;
            damage_reg_lambda=damage_reg_lam[damage_reg_level];
            damage_reg_on=1;
            damage_reg_napply=0;
          }
          }
          }
          /* [DAMAGE TR] level 3.  Levels 1 and 2 have already run and
             failed on THIS wall - unchanged Rescue2 - so the trajectory up
             to here is the Rescue2 trajectory.  Only now is the increment
             given one more attempt, with the trust region choosing every
             correction inside it. */
          damage_dl_lasthelp=iinc;
          damage_dl_selfrec=0;
          damage_dl_on=0;
          if((damage_dl_mode==1)&&(damage_rescue_used>=3)&&
             (damage_rec_disarmed==0)){
            if(damage_dl_have<0){
              printf("[DAMAGE TR] not re-arming: the transpose check has "
                     "already failed once in this run.  The wall goes to "
                     "the original stock stop.%s","\n");
              fflush(stdout);
              damage_rescue_used=damage_rescue_maxlevel;
            }else if(damage_dl_narm>=damage_dl_maxarm){
              printf("[DAMAGE TR] ATTEMPT BUDGET EXHAUSTED: %" ITGFORMAT
                     " armed attempts of %" ITGFORMAT " used.  The mechanism "
                     "DISARMS, the state is left as the standard cutback "
                     "rollback restored it, and this wall goes to the "
                     "ORIGINAL stock stop.%s",
                     damage_dl_narm,damage_dl_maxarm,"\n");
              fflush(stdout);
              damage_rec_disarmed=1;
              damage_rescue_used=damage_rescue_maxlevel;
            }else{
              if(damage_dl_narm>0){
                printf("[DAMAGE TR] inc=%" ITGFORMAT " running total "
                       "BEFORE this attempt: "
                       "armed %" ITGFORMAT ", accepted steps %" ITGFORMAT
                       " (Newton %" ITGFORMAT ", Cauchy %" ITGFORMAT
                       ", dogleg %" ITGFORMAT "), iterations that accepted "
                       "nothing %" ITGFORMAT ", rejected trials %" ITGFORMAT
                       ", residual evaluations %" ITGFORMAT
                       ", armed factorisations %" ITGFORMAT "%s",
                       iinc,damage_dl_narm,damage_dl_nacc,damage_dl_nnewt,
                       damage_dl_ncau,damage_dl_ndog,damage_dl_nfail,
                       damage_dl_nrej,damage_dl_neval,damage_dl_nfact,"\n");
                fflush(stdout);
              }
              damage_dl_narm++;
              damage_dl_on=1;
              damage_dl_have=0;
              damage_dl_delta=0.;
              damage_dl_banner=0;
              damage_dl_incarm=iinc;
              printf("[DAMAGE TR] LEVEL 3 armed for inc=%" ITGFORMAT
                     " (attempt %" ITGFORMAT " of %" ITGFORMAT ").  Levels 1 "
                     "and 2 both failed on this wall; idamagereeq=%" ITGFORMAT
                     ".  The increment is retried once more at the last "
                     "admissible dtheta, and inside it BK3 is replaced by a "
                     "dogleg trust region on phi=1/2|R|^2.  Nothing else "
                     "changes: no matrix entry, no material constant, no "
                     "deletion rule, and convergence is still decided by "
                     "checkconvergence on the unmodified residual.%s",
                     iinc,damage_dl_narm,damage_dl_maxarm,idamagereeq,"\n");
              fflush(stdout);
            }
          }
          damage_rescue_nfired++;
          dtheta=damage_rescue_dtheta_last;
          dthetaref=damage_rescue_dthetaref_last;
          printf("[DAMAGE RESCUE] FIRED #%" ITGFORMAT " inc=%" ITGFORMAT
                 " iter=%" ITGFORMAT " time=%.12e dtime=%.12e; dtheta restored "
                 "to the last admissible %.12e; transactional BT armed for "
                 "THIS attempt only (level %" ITGFORMAT ")%s",
                 damage_rescue_nfired,iinc,iit,time,dtime,dtheta,
                 damage_rescue_used,"\n");
          fflush(stdout);
        }

        /* checkconvergence advances theta itself, so an ACCEPTED increment is
           exactly one where theta moved.  Commit lambda there and nowhere
           else: the first version committed it in the dissipation-report
           block, gated on icutb==0 && idamagereeq==0, so every cutback and
           every re-equilibration advanced theta while leaving lambda behind.
           The drift accumulated - lambda 0.2733 against theta 0.2706 - and the
           run stalled at 98.6% of peak with 38 deletions. */
        if((damage_arc==1)&&(theta>damage_arc_theta0)){
          /* Until the constraint engages, lambda has no equation of its own
             and must track theta EXACTLY.  It cannot be taken from
             damage_diss_lamcur: that was built at the top of the increment
             from the TRIAL dtheta, while checkconvergence advances theta by
             whatever dtheta it settles on.  The two drift, the dG formula
             below then pairs a current theta with a stale lambda, dG comes
             out wrong and the engagement threshold never fires - measured,
             the arc arm never engaged while the identical run without it
             engaged and reached 70.9% of peak (E-96). */
          if(damage_diss_engaged==1){
            damage_arc_lam=damage_diss_lamcur;
          }else{
            damage_arc_lam+=theta-damage_arc_theta0;
          }
          damage_diss_lamcur=damage_arc_lam;
        }

        ctrl[3]=icref;
        uam[0]=damage_reeq_uam_actual[0];
        uam[1]=damage_reeq_uam_actual[1];

	if(*mortar>1){
	  SFREE(f_cs);SFREE(f_cm);
	} 
	  
      }else{

	/* explicit dynamics */

	icntrl=1;
	icutb=0;   

	theta=theta+dtheta;  
	if(dtheta>=1.-theta){
	  if(dtheta>1.-theta){
	    printf(" the increment size exceeds the remainder of the step and is decreased to %e\n\n",
		   dtheta**tper);
	  }
	  dtheta=1.-theta;
	  dthetaref=dtheta;
	}
	iflagact=0;
      }
    }

    if((*mortar==-1)&&(masslesslinear==0)&&(ncont!=0))
      {SFREE(auw);SFREE(jqw);SFREE(iroww);}

    if(*nmethod!=4)SFREE(resold);

    /*********************************************************/
    /*   end of the iteration loop                          */
    /*********************************************************/

    /* icutb=0 means that the iterations in the increment converged,
       icutb!=0 indicates that the increment has to be reiterated with
       another increment size (dtheta) */

    if(*mortar>1){
      SFREE(aubd);SFREE(jqbd);SFREE(irowbd);
      SFREE(aubdtil);SFREE(jqbdtil);SFREE(irowbdtil);
      SFREE(aubdtil2);SFREE(jqbdtil2);SFREE(irowbdtil2);
      SFREE(audd);SFREE(jqdd);SFREE(irowdd);
      SFREE(auddinv);SFREE(jqddinv);SFREE(irowddinv);
      SFREE(auddtil);SFREE(jqddtil);SFREE(irowddtil);
      SFREE(auddtil2);SFREE(jqddtil2);SFREE(irowddtil2);
      SFREE(bhat);	  
      SFREE(islavactdof);
    }

    if((icutb==0)&&(idamagereeq==1)){

      /* The same-load Newton solve on the current eroded topology has
         converged.  Recompute damage for every surviving element from the
         immutable physical-increment baseline.  This closes the coupling

              equilibrium -> damage -> topology -> equilibrium

         at one load level and, because dambase=damdamageini, avoids
         double-counting xstate-xstateini on repeated active-set passes. */

      damage_mode=1;
      idamage=0;
      damage_alphaevent=2.;
      damagebase=damdamageini;

      isiz=mi[0]*ne0;
      cpypardou(damde1prev,dam,&isiz,&num_cpus);

      FORTRAN(calcdamagebase,(ipkon,lakon,kon,co,mi,thicke,
                          ielmat,ielprop,prop,&ne0,ndmat_,ntmat_,
                          ndmcon,dmcon,dam,damagebase,&dtime,sti,
                          ithermal,t1,xstate,xstateini,nstate_,vold,
                          &idamage,&damage_mode,&damage_alphaevent));

      damage_de1_stats(dam,damde1prev,ipkon,lakon,ne0,mi[0],
                       &damage_de1_nactive,&damage_de1_gt01,
                       &damage_de1_gt05,&damage_de1_gt09,
                       &damage_de1_nfull,&damage_de1_nchanged,
                       &damage_de1_dmax,&damage_de1_maxdelta);

      /* DE1.3 terminal active-set extension.  The just-converged same-load
         equilibrium may have driven additional surviving DE1.2 elements
         to terminal degradation.  Mark only a bounded batch; the next
         same-load solve will expose any further terminal candidates. */
      damage_de13_new=0;
      damage_de13_batch_dmax=0.;
      if(damage_de12_enabled){
        damage_de13_new=damage_de13_mark_terminal(
            dam,ipkon,lakon,ielmat,mi[2],ndmcon,dmcon,*ndmat_,*ntmat_,
            ne0,mi[0],damage_de13_delete_d,DAMAGE_DE13_BATCH_MAX,
            &damage_de13_batch_dmax,damage_de13_trigger_value,
            damage_de13_trigger_ip,matname,damage_delete_filter,
            damage_damvisc,damage_delete_visc,&damage_de13_batch_vmin);
        damage_de13_term_only=damage_de13_new;
        if((damage_deadall_g>0.)&&(damage_de13_new<DAMAGE_DE13_BATCH_MAX)){
          ITG nda=damage_de13_mark_deadall(
              dam,damage_damvisc,damage_delete_visc,ipkon,lakon,kon,*nk,
              *ne,ne0,mi[0],damage_deadall_g,
              DAMAGE_DE13_BATCH_MAX-damage_de13_new,&damage_deadall_nodes);
          if(nda>0){
            damage_de13_new+=nda;
            damage_deadall_total+=nda;
            printf("[DAMAGE DEADALL] inc=%" ITGFORMAT " time=%.12e "
                   "nodes=%" ITGFORMAT " deleted=%" ITGFORMAT " total=%"
                   ITGFORMAT "\n",iinc,theta**tper,damage_deadall_nodes,
                   nda,damage_deadall_total);
            fflush(stdout);
          }
        }
        if((damage_deadsole_g>0.)&&(damage_de13_new<DAMAGE_DE13_BATCH_MAX)){
          ITG nds=damage_de13_mark_deadsole(
              dam,damage_damvisc,damage_delete_visc,ipkon,lakon,kon,*nk,
              ne0,mi[0],damage_deadsole_g,
              DAMAGE_DE13_BATCH_MAX-damage_de13_new);
          if(nds>0){
            damage_de13_new+=nds;
            damage_deadsole_total+=nds;
            printf("[DAMAGE DEADSOLE] inc=%" ITGFORMAT " time=%.12e "
                   "deleted=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                   iinc,theta**tper,nds,damage_deadsole_total);
            fflush(stdout);
          }
        }
        if(damage_de13_new>0){
          damage_de13_transaction=1;
          idamage+=damage_de13_new;
          printf("[DAMAGE DE1.3 EXTEND] inc=%" ITGFORMAT
                 " pass=%" ITGFORMAT " time=%.12e new_terminal=%"
                 ITGFORMAT " batch_Dmax=%.6e batch_Dvis=%.6e\n",
                 iinc,damage_active_pass+1,theta**tper,damage_de13_new,
                 damage_de13_batch_dmax,damage_de13_batch_vmin);
          fflush(stdout);
        }
      }

      if(idamage>0){

        /* Redistribution on the current topology has caused additional
           elements to reach the damage limit at the SAME load level.
           These are not committed yet.  Extend the tentative topology,
           rebuild the equation structure and equilibrate once more. */

        damage_active_pass++;

        /* Detached-island sweep.  BK4 only asks whether the nodes of the
           elements just deleted lost all their elements; a region that is
           internally connected but no longer attached to any support
           passes that test and leaves the operator singular in its
           rigid-body modes.  Marking those elements terminal here, before
           the tentative batch is collected, lets them travel through the
           unchanged transactional path. */

        if(damage_de12_enabled){
          damage_float_new=0;
          damage_float_reach=0;
          damage_float_coh=0;
          damage_float_isl=0;
          FORTRAN(damfloat,(ipkon,kon,lakon,ne,&ne0,nk,nodeboun,nboun,
                            ipompc,nodempc,nmpc,&damage_float_batch,
                            &damage_float_new,&damage_float_reach,
                            &damage_float_coh,&damage_float_isl));
          if(damage_float_new>0){
            damage_float_total+=damage_float_new;
            /* islands_marked and total keep their old meaning - the sum of
               both removals - so stored logs and every parser of them stay
               readable.  facets/islands are ADDED, never substituted: a
               cohesive facet leaves at gmin*Kn, a bulk island leaves at
               whatever it was carrying, and one number for both hid that. */
            printf("[DAMAGE DETACHED] inc=%" ITGFORMAT " time=%.12e "
                   "islands_marked=%" ITGFORMAT " nodes_with_load_path=%"
                   ITGFORMAT " total=%" ITGFORMAT
                   " facets=%" ITGFORMAT " islands=%" ITGFORMAT "\n",
                   iinc,theta**tper,damage_float_new,damage_float_reach,
                   damage_float_total,damage_float_coh,damage_float_isl);
            fflush(stdout);
          }

          /* Fully failed cohesive facets, removed on the same tentative
             topology and inside the same batch, so they re-equilibrate
             together with the bulk deletions instead of in a pass of their
             own.  Default OFF; see CCX_DAMAGE_FACET_DELETE above. */
          if(damage_facetdel){
            damage_facetdel_new=0;
            FORTRAN(damfaildead,(ipkon,lakon,ne,xstate,nstate_,mi,
                                 &damage_float_batch,&damage_facetdel_new));
            if(damage_facetdel_new>0){
              damage_facetdel_total+=damage_facetdel_new;
              printf("[DAMAGE FACET DEAD] inc=%" ITGFORMAT " time=%.12e "
                     "facets_removed=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_facetdel_new,
                     damage_facetdel_total);
              fflush(stdout);
            }
          }

          /* How far the assembled nodal stiffness has fallen.  Reported
             before anything acts on it: the point is to find out whether
             the node the solver actually throws is separable by this
             measure, which mesh topology and geometry both failed to
             do. */
          if((damage_stiff_probe>0)&&(damage_addiag!=NULL)){
            damage_stiff_n1=0;damage_stiff_n2=0;damage_stiff_n3=0;
            damage_stiff_nneg=0;
            damage_stiff_worst=-1;damage_addmin=2.;
            for(i=0;i<*nk;i++){
              if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
              if(damage_addiag0[i]<=0.) continue;
              if(damage_addiag[i]<=0.){damage_stiff_nneg++;continue;}
              damage_addrat=damage_addiag[i]/damage_addiag0[i];
              if(damage_addrat<1.e-3) damage_stiff_n1++;
              if(damage_addrat<1.e-2) damage_stiff_n2++;
              if(damage_addrat<1.e-1) damage_stiff_n3++;
              if(damage_addrat<damage_addmin){
                damage_addmin=damage_addrat;damage_stiff_worst=i+1;
              }
            }
            printf("[DAMAGE STIFFNESS] inc=%" ITGFORMAT " time=%.12e "
                   "below_1e-3=%" ITGFORMAT " below_1e-2=%" ITGFORMAT
                   " below_1e-1=%" ITGFORMAT " nonpositive=%" ITGFORMAT
                   " worst=node_%" ITGFORMAT "_at_%.4e\n",
                   iinc,theta**tper,damage_stiff_n1,damage_stiff_n2,
                   damage_stiff_n3,damage_stiff_nneg,
                   damage_stiff_worst,damage_addmin);
            fflush(stdout);
          }

          /* Material left hanging on a single element.  BK4 tests for
             zero elements and the sweep above tests for a lost load
             path; a node that falls from 22 elements to 1 passes both,
             and the next solve throws it most of a specimen length.
             See damdangle.f for the measurement this comes from. */
          if(damage_bare>damage_bare_rep){
            damage_bare_rep=damage_bare;
            printf("[DAMAGE BARE NODE] inc=%" ITGFORMAT " time=%.12e "
                   "%" ITGFORMAT " node(s) hold only cohesive facets, "
                   "every bulk element gone (new maximum)\n",
                   iinc,theta**tper,damage_bare);
            fflush(stdout);
          }
          /* A node with no bulk left is only in trouble if the facets
             still holding it are themselves debonded.  Measured: DHC1
             carries 30 bare nodes from t=0.209 and runs 300 more
             increments, and the UC6 disk completes with 63 of them, so
             bare alone is not the discriminator.  g = max(gmin,1-dvisc)
             from xstate slot 2 says whether a facet is a real support. */
          if(damage_free_probe==1){
            /* H4 probe (E-24).  For every node that has lost all its bulk,
               assemble the local support operator from its LIVE cohesive
               facets

                 K = sum_ip (A/3) * shape_i(ip)^2 * g_ip * Kn *
                     ( n (x) n + beta*(t1 (x) t1 + t2 (x) t2) )

               and report lambda_min(K)/lambda_min(K_intact), where K_intact
               is the same sum with g=1.

               This is area- and direction-weighted and additive, which the
               previous version of this probe was not.  Max over a node's
               facets is the WRONG statistic (E-15 correction 2: one
               surviving integration point at 0.698 masked seventeen dead
               ones), and a plain sum cannot see rank collapse (H3): with one
               point bonded out of eighteen the operator has one stiff
               direction and two at gmin, which only an eigenvalue exposes.

               Reported as a distribution, because no threshold may be
               proposed until PASS and FAIL separate in one. */
            ITG h4i,h4j,h4k,h4m,h4n,h4np,h4a,h4b,h4p;
            double h4e1[3],h4e2[3],h4cr[3],h4nr[3],h4t1[3],h4t2[3],h4T[9];
            double h4nc,h4n1,h4area,h4kn,h4tn,h4ts,h4gm,h4bta,h4g,h4w,h4sh;
            double h4A[9],h4V[9],h4ev[3],h4r,h4off,h4th,h4c,h4s,h4t,h4tau;
            double *h4loc=NULL,*h4int=NULL,*h4rat=NULL;
            ITG h4cnt,h4tot;

            NNEW(damage_free_nb,ITG,*nk);
            NNEW(h4loc,double,9**nk);
            NNEW(h4int,double,9**nk);
            NNEW(h4rat,double,*nk);
            for(h4i=0;h4i<*nk;h4i++){damage_free_nb[h4i]=0;h4rat[h4i]=-1.;}
            for(h4i=0;h4i<9**nk;h4i++){h4loc[h4i]=0.;h4int[h4i]=0.;}

            /* live bulk per node */
            for(h4i=0;h4i<*ne;h4i++){
              if(ipkon[h4i]<0) continue;
              if(lakon[8*h4i]!='C') continue;
              h4np=4;
              if(lakon[8*h4i+3]=='6') h4np=6;
              else if(lakon[8*h4i+3]=='8') h4np=8;
              else if(lakon[8*h4i+3]=='1') h4np=10;
              for(h4j=0;h4j<h4np;h4j++){
                h4k=kon[ipkon[h4i]+h4j]-1;
                if((h4k>=0)&&(h4k<*nk)) damage_free_nb[h4k]++;
              }
            }

            /* local support operator from live cohesive facets */
            for(h4i=0;h4i<*ne;h4i++){
              if(ipkon[h4i]<0) continue;
              if(lakon[8*h4i]!='U') continue;
              h4a=kon[ipkon[h4i]]-1;
              h4b=kon[ipkon[h4i]+1]-1;
              h4n=kon[ipkon[h4i]+2]-1;
              if((h4a<0)||(h4b<0)||(h4n<0)) continue;
              for(h4j=0;h4j<3;h4j++){
                h4e1[h4j]=co[3*h4b+h4j]-co[3*h4a+h4j];
                h4e2[h4j]=co[3*h4n+h4j]-co[3*h4a+h4j];
              }
              h4cr[0]=h4e1[1]*h4e2[2]-h4e1[2]*h4e2[1];
              h4cr[1]=h4e1[2]*h4e2[0]-h4e1[0]*h4e2[2];
              h4cr[2]=h4e1[0]*h4e2[1]-h4e1[1]*h4e2[0];
              h4nc=sqrt(h4cr[0]*h4cr[0]+h4cr[1]*h4cr[1]+h4cr[2]*h4cr[2]);
              h4n1=sqrt(h4e1[0]*h4e1[0]+h4e1[1]*h4e1[1]+h4e1[2]*h4e1[2]);
              if((h4nc<1.e-30)||(h4n1<1.e-30)) continue;
              h4area=0.5*h4nc;
              for(h4j=0;h4j<3;h4j++){
                h4nr[h4j]=h4cr[h4j]/h4nc;
                h4t1[h4j]=h4e1[h4j]/h4n1;
              }
              h4t2[0]=h4nr[1]*h4t1[2]-h4nr[2]*h4t1[1];
              h4t2[1]=h4nr[2]*h4t1[0]-h4nr[0]*h4t1[2];
              h4t2[2]=h4nr[0]*h4t1[1]-h4nr[1]*h4t1[0];
              h4p=ielprop[h4i];
              h4kn=prop[h4p];h4tn=prop[h4p+1];h4ts=prop[h4p+2];h4gm=prop[h4p+4];
              if((h4kn<=0.)||(h4tn<=0.)) continue;
              h4bta=(h4ts/h4tn)*(h4ts/h4tn);
              for(h4j=0;h4j<3;h4j++)
                for(h4k=0;h4k<3;h4k++)
                  h4T[3*h4j+h4k]=h4nr[h4j]*h4nr[h4k]
                    +h4bta*(h4t1[h4j]*h4t1[h4k]+h4t2[h4j]*h4t2[h4k]);
              for(h4m=0;h4m<3;h4m++){
                h4g=1.-xstate[*nstate_*(mi[0]*h4i+h4m)+1];
                if(h4g<h4gm) h4g=h4gm;
                for(h4n=0;h4n<6;h4n++){
                  h4sh=(h4n%3==h4m)?(2./3.):(1./6.);
                  h4k=kon[ipkon[h4i]+h4n]-1;
                  if((h4k<0)||(h4k>=*nk)) continue;
                  h4w=h4area/3.*h4sh*h4sh*h4kn;
                  for(h4j=0;h4j<9;h4j++){
                    h4loc[9*h4k+h4j]+=h4w*h4g*h4T[h4j];
                    h4int[9*h4k+h4j]+=h4w*h4T[h4j];
                  }
                }
              }
            }

            /* lambda_min of each 3x3 by cyclic Jacobi, for bare nodes only */
            h4tot=0;
            for(h4i=0;h4i<*nk;h4i++){
              if(damage_free_nb[h4i]>0) continue;
              if(h4int[9*h4i]+h4int[9*h4i+4]+h4int[9*h4i+8]<=0.) continue;
              h4tot++;
              for(h4b=0;h4b<2;h4b++){
                for(h4j=0;h4j<9;h4j++)
                  h4A[h4j]=(h4b==0)?h4loc[9*h4i+h4j]:h4int[9*h4i+h4j];
                for(h4n=0;h4n<20;h4n++){
                  h4off=0.;
                  for(h4j=0;h4j<3;h4j++)
                    for(h4k=0;h4k<3;h4k++)
                      if(h4j!=h4k) h4off+=h4A[3*h4j+h4k]*h4A[3*h4j+h4k];
                  if(h4off<1.e-30) break;
                  for(h4j=0;h4j<2;h4j++){
                    for(h4k=h4j+1;h4k<3;h4k++){
                      if(fabs(h4A[3*h4j+h4k])<1.e-300) continue;
                      h4th=(h4A[3*h4k+h4k]-h4A[3*h4j+h4j])
                           /(2.*h4A[3*h4j+h4k]);
                      h4t=(h4th>=0.?1.:-1.)/(fabs(h4th)+sqrt(h4th*h4th+1.));
                      h4c=1./sqrt(h4t*h4t+1.);h4s=h4t*h4c;
                      for(h4m=0;h4m<3;h4m++){
                        h4tau=h4A[3*h4j+h4m];
                        h4A[3*h4j+h4m]=h4c*h4tau-h4s*h4A[3*h4k+h4m];
                        h4A[3*h4k+h4m]=h4s*h4tau+h4c*h4A[3*h4k+h4m];
                      }
                      for(h4m=0;h4m<3;h4m++){
                        h4tau=h4A[3*h4m+h4j];
                        h4A[3*h4m+h4j]=h4c*h4tau-h4s*h4A[3*h4m+h4k];
                        h4A[3*h4m+h4k]=h4s*h4tau+h4c*h4A[3*h4m+h4k];
                      }
                    }
                  }
                }
                h4ev[h4b]=h4A[0];
                if(h4A[4]<h4ev[h4b]) h4ev[h4b]=h4A[4];
                if(h4A[8]<h4ev[h4b]) h4ev[h4b]=h4A[8];
              }
              if(h4ev[1]>0.) h4rat[h4i]=h4ev[0]/h4ev[1];
            }

            /* distribution, not a threshold */
            if(h4tot>0){
              ITG h4bin[7];double h4edge[6];
              h4edge[0]=1.e-6;h4edge[1]=1.e-5;h4edge[2]=1.e-4;
              h4edge[3]=1.e-3;h4edge[4]=1.e-2;h4edge[5]=1.e-1;
              for(h4j=0;h4j<7;h4j++) h4bin[h4j]=0;
              h4r=2.;h4cnt=-1;
              for(h4i=0;h4i<*nk;h4i++){
                if(h4rat[h4i]<0.) continue;
                for(h4j=0;h4j<6;h4j++) if(h4rat[h4i]<h4edge[h4j]) break;
                h4bin[h4j]++;
                if(h4rat[h4i]<h4r){h4r=h4rat[h4i];h4cnt=h4i;}
              }
              printf("[H4 SUPPORT] inc=%" ITGFORMAT " time=%.6e bare_nodes=%"
                     ITGFORMAT "\n",iinc,theta**tper,h4tot);
              printf("[H4 SUPPORT]   lambda_min(K)/lambda_min(K_intact) bins:"
                     " <1e-6:%" ITGFORMAT " <1e-5:%" ITGFORMAT
                     " <1e-4:%" ITGFORMAT " <1e-3:%" ITGFORMAT
                     " <1e-2:%" ITGFORMAT " <1e-1:%" ITGFORMAT
                     " >=1e-1:%" ITGFORMAT "\n",
                     h4bin[0],h4bin[1],h4bin[2],h4bin[3],h4bin[4],h4bin[5],
                     h4bin[6]);
              if(h4cnt>=0)
                printf("[H4 SUPPORT]   worst node %" ITGFORMAT
                       " at (%.3f, %.3f, %.3f) ratio=%.4e\n",
                       h4cnt+1,co[3*h4cnt],co[3*h4cnt+1],co[3*h4cnt+2],h4r);
              fflush(stdout);
            }
            SFREE(h4rat);SFREE(h4int);SFREE(h4loc);SFREE(damage_free_nb);
          }
          if((damage_dangle_max>0)||(damage_stiff_min>0.)){
            damage_dangle_new=0;
            damage_dangle_weak=0;
            if(damage_stiff_haz==NULL) NNEW(damage_stiff_haz,ITG,*nk);
            for(i=0;i<*nk;i++) damage_stiff_haz[i]=0;
            if((damage_stiff_min>0.)&&(damage_addiag!=NULL)){
              for(i=0;i<*nk;i++){
                if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
                if(damage_addiag0[i]<=0.) continue;
                if(damage_addiag[i]<=0.) continue;
                if(damage_addiag[i]<damage_stiff_min*damage_addiag0[i])
                  damage_stiff_haz[i]=1;
              }
            }
            FORTRAN(damdangle,(ipkon,kon,lakon,ne,nk,
                               &damage_float_batch,&damage_dangle_max,
                               damage_stiff_haz,
                               &damage_dangle_new,&damage_dangle_weak,
                               &damage_bare));
            if(damage_dangle_new>0){
              damage_dangle_total+=damage_dangle_new;
              printf("[DAMAGE DANGLING] inc=%" ITGFORMAT " time=%.12e "
                     "elements_marked=%" ITGFORMAT " nodes_stranded=%"
                     ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_dangle_new,
                     damage_dangle_weak,damage_dangle_total);
              fflush(stdout);
            }
          }
        }

        damage_scan_count=0;
        for(i=0;i<ne0;i++){
          if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
            damage_scan_count++;
          }
        }

        /* DE1 changes the constitutive stiffness but deliberately keeps
           the element topology and sparse pattern fixed.  If calcdamage
           requests another equilibrium and no topology changed, simply
           repeat Newton at the SAME load level. */
        if(damage_scan_count==0){
          damage_soft_reeq=1;

          printf("[DAMAGE DE1] pass=%" ITGFORMAT
                 " inc=%" ITGFORMAT " time=%.12e changed=%" ITGFORMAT
                 " action=same-load-reequilibrate\n",
                 damage_active_pass,iinc,theta**tper,idamage);
          printf("[DAMAGE DE1 STATUS] pass=%" ITGFORMAT
                 " active=%" ITGFORMAT " D>0.1=%" ITGFORMAT
                 " D>0.5=%" ITGFORMAT " D>0.9=%" ITGFORMAT
                 " Dfull=%" ITGFORMAT " Dmax=%.6e max_dD=%.6e "
                 "changed_tol=%" ITGFORMAT "\n",
                 damage_active_pass,damage_de1_nactive,damage_de1_gt01,
                 damage_de1_gt05,damage_de1_gt09,damage_de1_nfull,
                 damage_de1_dmax,damage_de1_maxdelta,
                 damage_de1_nchanged);
          fflush(stdout);

          /* DE1.1 closes the staggered fixed point on the largest scalar
             damage correction, not on the number of elements whose last
             bits changed.  The latter can stay O(1e5) for many passes even
             when every individual correction is negligible. */
          if(damage_de1_maxdelta<=DAMAGE_DE1_FP_TOL){
            printf("[DAMAGE DE1 FP CONVERGED] inc=%" ITGFORMAT
                   " pass=%" ITGFORMAT " max_dD=%.6e tol=%.6e\n",
                   iinc,damage_active_pass,damage_de1_maxdelta,
                   DAMAGE_DE1_FP_TOL);
            fflush(stdout);
            goto damage_active_set_closed;
          }

          if(damage_active_pass>=DAMAGE_DE1_MAX_PASSES){

            if(damage_de1_maxdelta<=DAMAGE_DE1_FP_RELAX_TOL){
              printf("[DAMAGE DE1 FP ACCEPT] inc=%" ITGFORMAT
                     " pass=%" ITGFORMAT
                     " max_dD=%.6e relaxed_tol=%.6e\n",
                     iinc,damage_active_pass,damage_de1_maxdelta,
                     DAMAGE_DE1_FP_RELAX_TOL);
              fflush(stdout);
              goto damage_active_set_closed;
            }

            /* Do not spend an unbounded number of PARDISO factorizations
               at one load level.  Reject this physical increment and let
               the existing transactional rollback retry a smaller step. */
            theta=thetadamage;
            dtheta=DAMAGE_DE1_CUTBACK_FACTOR*dthetadamage;
            if(dtheta<*tmin) dtheta=*tmin;
            dthetaref=dtheta;
            istab=0;
            icutb=1;

            printf("[DAMAGE DE1 CUTBACK] inc=%" ITGFORMAT
                   " pass=%" ITGFORMAT " max_dD=%.6e > %.6e "
                   "old_dt=%.6e retry_dt=%.6e\n",
                   iinc,damage_active_pass,damage_de1_maxdelta,
                   DAMAGE_DE1_FP_RELAX_TOL,dthetadamage**tper,
                   dtheta**tper);
            fflush(stdout);
            goto damage_controller_done;
          }

          theta=thetadamage;
          dtheta=dthetadamage;
          dthetaref=dthetarefdamage;
          idiscon=1;

          /* idamagereeq stays set: after convergence calcdamagebase will
             recompute DE1 from the immutable physical-increment baseline
             and close the constitutive fixed point. */
          continue;
        }

        damage_soft_reeq=0;

        /* Rebuild the tentative history from the complete difference
           between the physical-increment baseline topology and the current
           active-set topology.  Earlier deleted elements retain the damage
           value stored when they were removed because calcdamage skips
           ipkon<0 elements. */

        if(damage_tent_elem!=NULL){
          SFREE(damage_tent_elem); damage_tent_elem=NULL;
          SFREE(damage_tent_mat); damage_tent_mat=NULL;
          SFREE(damage_tent_ip); damage_tent_ip=NULL;
          SFREE(damage_tent_value); damage_tent_value=NULL;
          damage_tent_count=0;
        }

        if(damage_scan_count>0){
          NNEW(damage_tent_elem,ITG,damage_scan_count);
          NNEW(damage_tent_mat,ITG,damage_scan_count);
          NNEW(damage_tent_ip,ITG,damage_scan_count);
          NNEW(damage_tent_value,double,damage_scan_count);

          damage_tent_step=*istep;
          damage_tent_increment=iinc;
          damage_tent_step_time=theta**tper;
          damage_tent_total_time=*ttime+damage_tent_step_time;
          damage_tent_count=0;

          for(i=0;i<ne0;i++){
            if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
              damage_tent_elem[damage_tent_count]=i+1;
              damage_tent_mat[damage_tent_count]=ielmat[mi[2]*i];

              damage_nip_local=damage_history_nip(&lakon[8*i],mi[0]);
              if(damage_nip_local<1) damage_nip_local=1;
              if(damage_nip_local>mi[0]) damage_nip_local=mi[0];

              damage_dmax=dam[mi[0]*i];
              damage_tent_ip[damage_tent_count]=1;
              for(j=1;j<damage_nip_local;j++){
                if(dam[mi[0]*i+j]>damage_dmax){
                  damage_dmax=dam[mi[0]*i+j];
                  damage_tent_ip[damage_tent_count]=j+1;
                }
              }
              imat=damage_tent_mat[damage_tent_count];
              if((imat>0)&&
                 damage_progressive_material(imat,ndmcon,dmcon,*ndmat_,*ntmat_)&&
                 (damage_dmax>1.)) damage_dmax-=1.;
              if((damage_de13_transaction)&&(damage_de13_trigger_value!=NULL)&&
                 (damage_de13_trigger_value[i]>=0.)){
                damage_tent_value[damage_tent_count]=
                    damage_de13_trigger_value[i];
                damage_tent_ip[damage_tent_count]=
                    damage_de13_trigger_ip[i];
              }else{
                damage_tent_value[damage_tent_count]=damage_dmax;
              }
              damage_tent_count++;
            }
          }
        }

        printf("[DAMAGE ACTIVESET] pass=%" ITGFORMAT
               " inc=%" ITGFORMAT " time=%.12e new_deleted=%" ITGFORMAT
               " total_tentative=%" ITGFORMAT "\n",
               damage_active_pass,iinc,theta**tper,idamage,
               damage_tent_count);
        fflush(stdout);

        /* Which elements are actually in the batch.  Four hypotheses for
           the stall at lambda=0.3579 died for want of exactly this: the
           batch is what decides whether the surviving material still
           carries load, and the batch is the one thing the log never
           said. */
        if(damage_batch_list==1){
          printf("[DAMAGE BATCH] pass=%" ITGFORMAT " inc=%" ITGFORMAT
                 " elements:",damage_active_pass,iinc);
          for(k=0;k<damage_tent_count;k++){
            printf(" %" ITGFORMAT,damage_tent_elem[k]);
          }
          printf("\n");
          fflush(stdout);
        }

        /* The batch as a SET: sorted, hashed, and printed in full.  Two
           batches with the same hash are the same set of elements
           whatever order the scan produced them in. */

        if((td_armed!=0)&&(td_trace!=0)){
          td_batch=topodiag_hash_batch(damage_tent_elem,damage_tent_count,
                                       td_sort,topodiag_hash_seed());
          printf("[BATCHTRACE] batch inc=%" ITGFORMAT " icutb=%" ITGFORMAT
                 " pass=%" ITGFORMAT " n=%" ITGFORMAT
                 " committed_state=%016llx batch=%016llx sorted:",
                 iinc,icutb,damage_active_pass,damage_tent_count,
                 td_state,td_batch);
          for(k=0;k<damage_tent_count;k++)
            printf(" %" ITGFORMAT,td_sort[k]);
          printf("\n");
          fflush(stdout);
        }

        /* BK4 deferred sparse-structure compaction.  A terminal element is
           already absent from assembly because ipkon<0.  Keeping the old
           sparse graph is algebraically harmless while every previously
           active node still belongs to at least one active element: the
           obsolete entries are simply assembled as zero.  Rebuild as soon
           as an orphan node appears.  Detached but internally connected
           components remain the responsibility of the planned connectivity
           manager; a failed same-load solve still follows the full rollback
           path. */

        damage_topology_rebuild=1;
        damage_topology_orphans=0;
        if((damage_topology_deferred_mode==1)&&
           (damage_de13_transaction==1)){
          NNEW(damage_iponoel_trial,ITG,*nk);
          ITGMEMSET(damage_iponoel_trial,0,*nk,0);
          FORTRAN(nodebelongstoel,(damage_iponoel_trial,lakon,ipkon,kon,ne));
          NNEW(damage_orphan_seen,ITG,*nk);
          for(k=0;k<damage_tent_count;k++){
            i=damage_tent_elem[k]-1;
            if((i<0)||(i>=ne0)||(ipkon[i]>=0)) continue;
            if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
            damage_indexe=-ipkon[i]-2;
            for(j=0;j<4;j++){
              inode=kon[damage_indexe+j]-1;
              if((inode<0)||(inode>=*nk)) continue;
              if((damage_iponoel_trial[inode]==0)&&
                 (damage_orphan_seen[inode]==0)){
                for(idir=1;idir<=3;idir++){
                  if(nactdof[mt*inode+idir]>0){
                    damage_topology_orphans++;
                    damage_orphan_seen[inode]=1;
                    break;
                  }
                }
              }
            }
          }
          if(damage_topology_orphans==0){
            isiz=*nk;
            cpyparitg(iponoel,damage_iponoel_trial,&isiz,&num_cpus);
            damage_topology_rebuild=0;
          }
          SFREE(damage_orphan_seen);
          SFREE(damage_iponoel_trial);
        }

        printf("[DAMAGE TOPOLOGY BK4] inc=%" ITGFORMAT
               " pass=%" ITGFORMAT " tentative=%" ITGFORMAT
               " orphan_nodes=%" ITGFORMAT " action=%s\n",
               iinc,damage_active_pass,damage_tent_count,
               damage_topology_orphans,
               damage_topology_rebuild?"remastruct":"reuse-sparse-graph");
        fflush(stdout);

        /* [DAMAGE RELEASE] arm.  f still holds f_int(u*) of the converged
           state - the elements were marked above but no results() has run
           since - and nactdof is still the pre-remastruct numbering.  Both
           are mapped into node space here so they survive remastruct. */
        if(damage_release_probe){
          ITG ai,aj,ak;
          if(damage_frel==NULL){
            NNEW(damage_frel,double,mt**nk);
            NNEW(damage_ract,ITG,mt**nk);
          }
          for(ai=0;ai<*nk;ai++){
            for(aj=0;aj<mt;aj++){
              ak=nactdof[mt*ai+aj];
              damage_frel[mt*ai+aj]=(ak>0)?f[ak-1]:0.;
              damage_ract[mt*ai+aj]=(ak>0)?1:0;
            }
          }
          damage_release_armed=1;
          damage_release_pass++;
          damage_release_qa=qa[0];
          damage_release_qam=qam[0];
          damage_release_dt=dtime;
          damage_release_rebuild=damage_topology_rebuild;
          damage_release_iforbou=iforbou;
          damage_release_nterm=damage_de13_term_only;
          damage_release_nother=damage_de13_new-damage_de13_term_only;
          damage_release_nisl=damage_float_isl;
          damage_release_ncoh=damage_float_coh;

          /* level 2: one line per element in this batch.  The degradation it
             carried when it left is printed as a REFERENCE CHARACTERISTIC -
             it is a dimensionless multiplier, not a force, and it is NOT
             summed into any estimate of the released force.  The measured
             force is the surv/removed split above.

             It does classify the source for free, though: the terminal
             trigger cannot fire below its own threshold, so an element that
             left at Dvis < delete_d did NOT come from it - it came from
             damfloat, which reads no damage variable at all and therefore
             removes at whatever the element was carrying. */
          if(damage_release_probe>=2){
            ITG pe,pj,pnip,pel;
            double pd,pdv,pdmax,pdvmax,ptrig;
            for(pe=0;pe<damage_tent_count;pe++){
              pel=damage_tent_elem[pe]-1;
              if((pel<0)||(pel>=ne0)) continue;
              pnip=damage_history_nip(&lakon[8*pel],mi[0]);
              if(pnip<1) pnip=1;
              if(pnip>mi[0]) pnip=mi[0];
              pdmax=0.; pdvmax=0.;
              for(pj=0;pj<pnip;pj++){
                pd=dam[mi[0]*pel+pj]-1.;
                if(pd<0.) pd=0.;
                if(pd>1.) pd=1.;
                if(pd>pdmax) pdmax=pd;
                if(damage_damvisc!=NULL){
                  pdv=damage_damvisc[mi[0]*pel+pj];
                  if(pdv<0.) pdv=0.;
                  if(pdv>1.) pdv=1.;
                  if(pdv>pdvmax) pdvmax=pdv;
                }
              }
              ptrig=((damage_delete_visc==1)&&(damage_damvisc!=NULL))
                    ?pdvmax:pdmax;
              /* The bulk damage variable exists ONLY for bulk elements.
                 calcdamage.f:135 skips every lakon(1:1) != 'C', so dam is
                 identically zero for a UC6 facet - printing it would read
                 as "left fully intact" when in truth the facet was
                 conducting gmin*Kn and its state lives in xstate, not here.
                 Both decks in play carry UC6, so this is not hypothetical:
                 damfloatcoh removals land in the same tentative batch. */
              if(lakon[8*pel]!='C'){
                printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                       " pass=%" ITGFORMAT " el=%" ITGFORMAT
                       " mat=%" ITGFORMAT " type=%.8s"
                       " D=n/a Dvis=n/a g_ref=n/a"
                       " note=non-bulk-element-damage-lives-in-xstate%s",
                       iinc,damage_release_pass,damage_tent_elem[pe],
                       damage_tent_mat[pe],&lakon[8*pel],
                       "\n");
                continue;
              }
              printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                     " pass=%" ITGFORMAT " el=%" ITGFORMAT
                     " mat=%" ITGFORMAT " type=%.8s"
                     " D=%.6f Dvis=%.6f g_ref=%.4e"
                     " below_delete_d=%s%s",
                     iinc,damage_release_pass,damage_tent_elem[pe],
                     damage_tent_mat[pe],&lakon[8*pel],pdmax,pdvmax,1.-ptrig,
                     (ptrig<damage_de13_delete_d)?"YES-not-terminal":"no",
                     "\n");
            }
            fflush(stdout);
          }
        }

        if(damage_topology_rebuild){
          iitsav=iit;
          iit=-2;

          remastruct(ipompc,&coefmpc,&nodempc,nmpc,
                     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
                     ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
                     &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
                     &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
                     &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
                     mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
                     itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
                     tieset,ntie,&num_cpus,ielmat,matname);

          iit=iitsav;

          SFREE(nactdofinv);
          NNEW(nactdofinv,ITG,mt**nk);
          MNEW(nodorig,ITG,*nk);

          FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
                                 ipkon,lakon,kon,ne));

          SFREE(nodorig);

          ITGMEMSET(iponoel,0,*nk,0);
          FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));
        }

        /* checkconvergence() advanced theta when the just-finished
           re-equilibration converged.  Rewind to the saved beginning so
           the next Newton solve reaches exactly the same physical load. */

        theta=thetadamage;
        dtheta=dthetadamage;
        dthetaref=dthetarefdamage;
        idiscon=1;

        /* idamagereeq deliberately remains 1: the active set is still
           open and the next converged Newton solve must be checked again. */
        continue;
      }

damage_active_set_closed:

      /* No additional element failed after recomputing damage from the
         final post-deletion equilibrium.  The damage/topology active set
         is closed and the whole tentative deletion set may now be
         committed atomically to jobname.damage. */

      if(damage_tent_count>0){
        damage_batch++;

        for(i=0;i<damage_tent_count;i++){
          fprintf(fdamage,
                  "%" ITGFORMAT " %" ITGFORMAT " %" ITGFORMAT
                  " %.15e %.15e %" ITGFORMAT " %.15e %" ITGFORMAT
                  " %" ITGFORMAT "\n",
                  damage_tent_elem[i],damage_tent_step,
                  damage_tent_increment,damage_tent_step_time,
                  damage_tent_total_time,damage_tent_mat[i],
                  damage_tent_value[i],damage_tent_ip[i],damage_batch);
        }

        fflush(fdamage);

        printf("[DAMAGE COMMIT] batch=%" ITGFORMAT
               " inc=%" ITGFORMAT " time=%.12e deleted=%" ITGFORMAT
               " active_passes=%" ITGFORMAT "\n",
               damage_batch,damage_tent_increment,
               damage_tent_step_time,damage_tent_count,
               damage_active_pass);
        if(damage_de13_transaction){
          printf("[DAMAGE DE1.3 COMMIT] inc=%" ITGFORMAT
                 " time=%.12e terminal_deleted=%" ITGFORMAT
                 " active_passes=%" ITGFORMAT "\n",
                 damage_tent_increment,damage_tent_step_time,
                 damage_tent_count,damage_active_pass);
        }
        fflush(stdout);

        /* [LOADPATH] The connectivity test that used to live here has
           moved to the per-increment census in loadpath.c.

           It was measurably weaker in the same three ways every time:
           it ran only inside this block, so it needed a bulk deletion in
           the very increment the specimen came apart; it needed
           CCX_FRACTURE_TERMINATION to name two node sets by hand; and it
           needed CCX_FRACTURE_DEADFACET on top of that before a fully
           failed facet stopped counting as a load path.  Same kernel, same
           facet judgement, three conditions instead of none.  Running both
           would be a seventh overlapping mechanism, so this one is gone
           rather than kept as insurance. */


        SFREE(damage_tent_elem); damage_tent_elem=NULL;
        SFREE(damage_tent_mat); damage_tent_mat=NULL;
        SFREE(damage_tent_ip); damage_tent_ip=NULL;
        SFREE(damage_tent_value); damage_tent_value=NULL;
        damage_tent_count=0;
      }

      if(damage_soft_reeq==1){
        printf("[DAMAGE DE1 CONVERGED] inc=%" ITGFORMAT
               " time=%.12e passes=%" ITGFORMAT
               " active=%" ITGFORMAT " Dmax=%.6e max_dD=%.6e\n",
               iinc,theta**tper,damage_active_pass,
               damage_de1_nactive,damage_de1_dmax,
               damage_de1_maxdelta);

        damage_de1_append_stats(jobnamec,*istep,iinc,theta**tper,
                                *ttime+theta**tper,damage_active_pass,
                                damage_de1_nactive,damage_de1_gt01,
                                damage_de1_gt05,damage_de1_gt09,
                                damage_de1_nfull,damage_de1_dmax,
                                damage_de1_maxdelta);

        damage_de1_write_vtk(jobnamec,co,vold,*nk,mt,kon,ipkon,lakon,
                             ielmat,mi[2],dam,mi[0],ne0,*istep,iinc,
                             theta**tper);

        printf("[DAMAGE DE1 OUTPUT] exact cell snapshot: %s.de1.vtk; "
               "history: %s.de1stats\n",jobnamec,jobnamec);
        fflush(stdout);
      }

      damage_soft_reeq=0;
      damage_de13_transaction=0;

      /* If a bounded fast event had to cut the physical increment,
         do not regrow dtheta through many 1.5x increments.  The local
         controller below will TRY the whole remaining distance to
         theta_goal in one Newton solve.  Failure is harmless: ordinary
         CalculiX cutback remains in charge. */
      if((damage_fast_used==1)&&(ilocalsubstep==1))
        damage_fast_recover=1;

      idamagereeq=0;

    }else if((icutb==0)&&(*ndmat_>0)){
      damage_event_cut=0;
      damage_mode=0;
      damage_predict_count=0;
      damage_alphaevent=2.;

      /* Use the immutable damage state at the beginning of the physical
         increment as the baseline.  In explicit dynamics no rollback
         baseline is allocated, so the current committed dam is used. */

      if((*iexpl<=1)&&(damdamageini!=NULL)){
        damagebase=damdamageini;
      }else{
        damagebase=dam;
      }

      /* First pass: predict whether the accumulated damage reaches the
         deletion limit inside this already-converged physical increment.
         This pass does not modify dam or ipkon. */

      FORTRAN(calcdamagebase,(ipkon,lakon,kon,co,mi,thicke,
			  ielmat,ielprop,prop,&ne0,ndmat_,ntmat_,
			  ndmcon,dmcon,dam,damagebase,&dtime,sti,ithermal,t1,xstate,
			  xstateini,nstate_,vold,&damage_predict_count,
			  &damage_mode,&damage_alphaevent));

      /* Adaptive-fast event controller.

         A2 localized essentially every first crossing, which is robust but
         can collapse dtheta by factors of 50-100 when many C3D4 elements
         approach the limit one after another.

         A3 treats event placement as a TRIAL only:
           1) if the already-converged coarse increment predicts a modest
              batch and the first event is not too early, accept the coarse
              endpoint and let the transactional active-set test the batch;
           2) otherwise localize, but never below a fixed fraction of the
              current physical increment on the first fast attempt;
           3) if the post-deletion Newton solve fails, the existing rollback
              restores topology/state and damage_fast_retry forces the next
              attempt back to the exact A2 locator (no floor, no batching).

         Thus aggressive batching can cost at most one failed trial.  It can
         never be committed unless same-load equilibrium and active-set
         closure both converge. */

      /* damage_fast_used is intentionally NOT cleared here: when a
         bounded-localization trial is rolled back and re-solved, the flag
         must survive until the resulting deletion transaction commits so
         the fast recovery-to-theta_goal logic can be armed. */

      if((damage_predict_count>0)&&(*iexpl<=1)&&
         (damage_alphaevent>1.e-8)&&(damage_alphaevent<0.95)&&
         (dthetadamage>1.01**tmin)){

        damage_event_raw=1.02*damage_alphaevent*dthetadamage;
        damage_event_dtheta=damage_event_raw;

        if(damage_fast_retry==0){

          /* Fast path A: use the already-converged coarse endpoint for a
             bounded predicted batch.  Actual deletion is still based on the
             re-solved/current damage state, not on predictor membership. */

          if((damage_predict_count<=DAMAGE_FAST_BATCH_MAX)&&
             (damage_alphaevent>=DAMAGE_FAST_DIRECT_ALPHA)){

            damage_fast_used=1;

            printf("[DAMAGE FAST] inc=%" ITGFORMAT
                   " action=coarse-batch alpha=%.6e dt=%.6e "
                   "predicted=%" ITGFORMAT "\n",
                   iinc,damage_alphaevent,dthetadamage**tper,
                   damage_predict_count);
            fflush(stdout);

          }else{

            /* Fast path B: bounded localization.  This prevents a single
               early crossing from forcing an extremely small physical step.
               Any excessive topology jump is protected by transactional
               rollback and the exact-locator retry below. */

            damage_event_floor=DAMAGE_FAST_MIN_FRACTION*dthetadamage;
            if(damage_event_dtheta<damage_event_floor)
              damage_event_dtheta=damage_event_floor;

            if(damage_event_dtheta<*tmin) damage_event_dtheta=*tmin;
            if(damage_event_dtheta>0.95*dthetadamage)
              damage_event_dtheta=0.95*dthetadamage;

            if(damage_event_dtheta<0.999*dthetadamage){
              damage_fast_used=1;

              printf("[DAMAGE FAST] inc=%" ITGFORMAT
                     " action=bounded-localize alpha=%.6e old_dt=%.6e "
                     "raw_dt=%.6e trial_dt=%.6e predicted=%" ITGFORMAT
                     "\n",
                     iinc,damage_alphaevent,dthetadamage**tper,
                     damage_event_raw**tper,damage_event_dtheta**tper,
                     damage_predict_count);
              fflush(stdout);

              theta=thetadamage;
              dtheta=damage_event_dtheta;
              dthetaref=dtheta;
              istab=0;
              icutb=1;
              damage_event_cut=1;
            }
          }

        }else{

          /* Safety path after a failed fast post-deletion equilibrium:
             restore the exact A2 locator.  No artificial lower bound is
             applied. */

          if(damage_event_dtheta<*tmin) damage_event_dtheta=*tmin;
          if(damage_event_dtheta>0.95*dthetadamage)
            damage_event_dtheta=0.95*dthetadamage;

          if(damage_event_dtheta<0.999*dthetadamage){
            printf("[DAMAGE SAFE RETRY] inc=%" ITGFORMAT
                   " alpha=%.6e old_dt=%.6e exact_dt=%.6e "
                   "predicted=%" ITGFORMAT " fast_failures=%" ITGFORMAT
                   "\n",
                   iinc,damage_alphaevent,dthetadamage**tper,
                   damage_event_dtheta**tper,damage_predict_count,
                   damage_fast_failures);
            fflush(stdout);

            theta=thetadamage;
            dtheta=damage_event_dtheta;
            dthetaref=dtheta;
            istab=0;
            icutb=1;
            damage_event_cut=1;
          }
        }
      }

      if(damage_event_cut==0){

        /* Second pass: apply the damage increment and retain the original
           hard-deletion law. Because strongly interior crossings were
           localized above, the simultaneous deletion batch should now be
           much smaller and closer to the first physical damage event. */

        damage_mode=1;
        idamage=0;
        damage_alphaevent=2.;

        isiz=mi[0]*ne0;
        cpypardou(damde1prev,dam,&isiz,&num_cpus);

        FORTRAN(calcdamagebase,(ipkon,lakon,kon,co,mi,thicke,
			    ielmat,ielprop,prop,&ne0,ndmat_,ntmat_,
			    ndmcon,dmcon,dam,damagebase,&dtime,sti,ithermal,t1,xstate,
			    xstateini,nstate_,vold,&idamage,
			    &damage_mode,&damage_alphaevent));

        damage_de1_stats(dam,damde1prev,ipkon,lakon,ne0,mi[0],
                         &damage_de1_nactive,&damage_de1_gt01,
                         &damage_de1_gt05,&damage_de1_gt09,
                         &damage_de1_nfull,&damage_de1_nchanged,
                         &damage_de1_dmax,&damage_de1_maxdelta);

      if((idamage>0)&&(*iexpl<=1)){

        /* Build a tentative deletion batch by comparing the topology
           snapshot at the beginning of the physical increment with the
           topology immediately after calcdamage().  Nothing is written
           to disk here: the batch is committed only after successful
           same-load re-equilibration. */

        /* Detached-island sweep.  BK4 only asks whether the nodes of the
           elements just deleted lost all their elements; a region that is
           internally connected but no longer attached to any support
           passes that test and leaves the operator singular in its
           rigid-body modes.  Marking those elements terminal here, before
           the tentative batch is collected, lets them travel through the
           unchanged transactional path. */

        if(damage_de12_enabled){
          damage_float_new=0;
          damage_float_reach=0;
          damage_float_coh=0;
          damage_float_isl=0;
          FORTRAN(damfloat,(ipkon,kon,lakon,ne,&ne0,nk,nodeboun,nboun,
                            ipompc,nodempc,nmpc,&damage_float_batch,
                            &damage_float_new,&damage_float_reach,
                            &damage_float_coh,&damage_float_isl));
          if(damage_float_new>0){
            damage_float_total+=damage_float_new;
            /* islands_marked and total keep their old meaning - the sum of
               both removals - so stored logs and every parser of them stay
               readable.  facets/islands are ADDED, never substituted: a
               cohesive facet leaves at gmin*Kn, a bulk island leaves at
               whatever it was carrying, and one number for both hid that. */
            printf("[DAMAGE DETACHED] inc=%" ITGFORMAT " time=%.12e "
                   "islands_marked=%" ITGFORMAT " nodes_with_load_path=%"
                   ITGFORMAT " total=%" ITGFORMAT
                   " facets=%" ITGFORMAT " islands=%" ITGFORMAT "\n",
                   iinc,theta**tper,damage_float_new,damage_float_reach,
                   damage_float_total,damage_float_coh,damage_float_isl);
            fflush(stdout);
          }

          /* Fully failed cohesive facets, removed on the same tentative
             topology and inside the same batch, so they re-equilibrate
             together with the bulk deletions instead of in a pass of their
             own.  Default OFF; see CCX_DAMAGE_FACET_DELETE above. */
          if(damage_facetdel){
            damage_facetdel_new=0;
            FORTRAN(damfaildead,(ipkon,lakon,ne,xstate,nstate_,mi,
                                 &damage_float_batch,&damage_facetdel_new));
            if(damage_facetdel_new>0){
              damage_facetdel_total+=damage_facetdel_new;
              printf("[DAMAGE FACET DEAD] inc=%" ITGFORMAT " time=%.12e "
                     "facets_removed=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_facetdel_new,
                     damage_facetdel_total);
              fflush(stdout);
            }
          }

          /* How far the assembled nodal stiffness has fallen.  Reported
             before anything acts on it: the point is to find out whether
             the node the solver actually throws is separable by this
             measure, which mesh topology and geometry both failed to
             do. */
          if((damage_stiff_probe>0)&&(damage_addiag!=NULL)){
            damage_stiff_n1=0;damage_stiff_n2=0;damage_stiff_n3=0;
            damage_stiff_nneg=0;
            damage_stiff_worst=-1;damage_addmin=2.;
            for(i=0;i<*nk;i++){
              if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
              if(damage_addiag0[i]<=0.) continue;
              if(damage_addiag[i]<=0.){damage_stiff_nneg++;continue;}
              damage_addrat=damage_addiag[i]/damage_addiag0[i];
              if(damage_addrat<1.e-3) damage_stiff_n1++;
              if(damage_addrat<1.e-2) damage_stiff_n2++;
              if(damage_addrat<1.e-1) damage_stiff_n3++;
              if(damage_addrat<damage_addmin){
                damage_addmin=damage_addrat;damage_stiff_worst=i+1;
              }
            }
            printf("[DAMAGE STIFFNESS] inc=%" ITGFORMAT " time=%.12e "
                   "below_1e-3=%" ITGFORMAT " below_1e-2=%" ITGFORMAT
                   " below_1e-1=%" ITGFORMAT " nonpositive=%" ITGFORMAT
                   " worst=node_%" ITGFORMAT "_at_%.4e\n",
                   iinc,theta**tper,damage_stiff_n1,damage_stiff_n2,
                   damage_stiff_n3,damage_stiff_nneg,
                   damage_stiff_worst,damage_addmin);
            fflush(stdout);
          }

          /* Material left hanging on a single element.  BK4 tests for
             zero elements and the sweep above tests for a lost load
             path; a node that falls from 22 elements to 1 passes both,
             and the next solve throws it most of a specimen length.
             See damdangle.f for the measurement this comes from. */
          if(damage_bare>damage_bare_rep){
            damage_bare_rep=damage_bare;
            printf("[DAMAGE BARE NODE] inc=%" ITGFORMAT " time=%.12e "
                   "%" ITGFORMAT " node(s) hold only cohesive facets, "
                   "every bulk element gone (new maximum)\n",
                   iinc,theta**tper,damage_bare);
            fflush(stdout);
          }
          /* A node with no bulk left is only in trouble if the facets
             still holding it are themselves debonded.  Measured: DHC1
             carries 30 bare nodes from t=0.209 and runs 300 more
             increments, and the UC6 disk completes with 63 of them, so
             bare alone is not the discriminator.  g = max(gmin,1-dvisc)
             from xstate slot 2 says whether a facet is a real support. */
          if(damage_free_probe==1){
            NNEW(damage_free_nb,ITG,*nk);
            NNEW(damage_free_g,double,*nk);
            for(i=0;i<*nk;i++){damage_free_nb[i]=0;damage_free_g[i]=-1.;}
            for(i=0;i<*ne;i++){
              if(ipkon[i]<0) continue;
              if(lakon[8*i]=='C'){
                for(j=0;j<4;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k>=0)&&(k<*nk)) damage_free_nb[k]++;
                }
              }else if(lakon[8*i]=='U'){
                damage_free_gm=0.;
                for(j=0;j<3;j++){
                  damage_free_dv=xstate[*nstate_*(mi[0]*i+j)+1];
                  if(1.-damage_free_dv>damage_free_gm)
                    damage_free_gm=1.-damage_free_dv;
                }
                for(j=0;j<6;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k<0)||(k>=*nk)) continue;
                  if(damage_free_gm>damage_free_g[k])
                    damage_free_g[k]=damage_free_gm;
                }
              }
            }
            damage_free_cnt=0;
            for(i=0;i<*nk;i++){
              if(damage_free_nb[i]>0) continue;
              if(damage_free_g[i]<0.) continue;
              if(damage_free_g[i]<1.e-2) damage_free_cnt++;
            }
            if(damage_free_cnt>damage_free_rep){
              damage_free_rep=damage_free_cnt;
              printf("[DAMAGE FREE NODE] inc=%" ITGFORMAT " time=%.12e "
                     "%" ITGFORMAT " node(s) with no bulk AND every "
                     "remaining facet debonded (new maximum)\n",
                     iinc,theta**tper,damage_free_cnt);
              fflush(stdout);
            }
            SFREE(damage_free_nb);SFREE(damage_free_g);
          }
          if((damage_dangle_max>0)||(damage_stiff_min>0.)){
            damage_dangle_new=0;
            damage_dangle_weak=0;
            if(damage_stiff_haz==NULL) NNEW(damage_stiff_haz,ITG,*nk);
            for(i=0;i<*nk;i++) damage_stiff_haz[i]=0;
            if((damage_stiff_min>0.)&&(damage_addiag!=NULL)){
              for(i=0;i<*nk;i++){
                if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
                if(damage_addiag0[i]<=0.) continue;
                if(damage_addiag[i]<=0.) continue;
                if(damage_addiag[i]<damage_stiff_min*damage_addiag0[i])
                  damage_stiff_haz[i]=1;
              }
            }
            FORTRAN(damdangle,(ipkon,kon,lakon,ne,nk,
                               &damage_float_batch,&damage_dangle_max,
                               damage_stiff_haz,
                               &damage_dangle_new,&damage_dangle_weak,
                               &damage_bare));
            if(damage_dangle_new>0){
              damage_dangle_total+=damage_dangle_new;
              printf("[DAMAGE DANGLING] inc=%" ITGFORMAT " time=%.12e "
                     "elements_marked=%" ITGFORMAT " nodes_stranded=%"
                     ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_dangle_new,
                     damage_dangle_weak,damage_dangle_total);
              fflush(stdout);
            }
          }
        }

        damage_scan_count=0;
        for(i=0;i<ne0;i++){
          if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
            damage_scan_count++;
          }
        }

        /* DE1: constitutive damage changed, but no element was physically
           deleted.  Keep the matrix graph intact and re-equilibrate the
           degraded material at the same physical load level. */
        if(damage_scan_count==0){
          damage_active_pass=1;
          damage_soft_reeq=1;
          damage_tent_increment=iinc;
          damage_tent_step_time=theta**tper;
          damage_tent_total_time=*ttime+damage_tent_step_time;

          printf("[DAMAGE DE1] pass=1 inc=%" ITGFORMAT
                 " time=%.12e changed=%" ITGFORMAT
                 " action=same-load-reequilibrate\n",
                 iinc,theta**tper,idamage);
          printf("[DAMAGE DE1 STATUS] pass=1 active=%" ITGFORMAT
                 " D>0.1=%" ITGFORMAT " D>0.5=%" ITGFORMAT
                 " D>0.9=%" ITGFORMAT " Dfull=%" ITGFORMAT
                 " Dmax=%.6e max_dD=%.6e changed_tol=%" ITGFORMAT "\n",
                 damage_de1_nactive,damage_de1_gt01,damage_de1_gt05,
                 damage_de1_gt09,damage_de1_nfull,damage_de1_dmax,
                 damage_de1_maxdelta,damage_de1_nchanged);
          fflush(stdout);

          theta=thetadamage;
          dtheta=dthetadamage;
          dthetaref=dthetarefdamage;
          idiscon=1;
          idamagereeq=1;
          continue;
        }

        damage_soft_reeq=0;

        /* A non-empty old buffer would indicate an internal logic error.
           Discard it rather than allowing trial data to leak into a later
           commit. */
        if(damage_tent_elem!=NULL){
          SFREE(damage_tent_elem); damage_tent_elem=NULL;
          SFREE(damage_tent_mat); damage_tent_mat=NULL;
          SFREE(damage_tent_ip); damage_tent_ip=NULL;
          SFREE(damage_tent_value); damage_tent_value=NULL;
          damage_tent_count=0;
        }

        if(damage_scan_count>0){
          NNEW(damage_tent_elem,ITG,damage_scan_count);
          NNEW(damage_tent_mat,ITG,damage_scan_count);
          NNEW(damage_tent_ip,ITG,damage_scan_count);
          NNEW(damage_tent_value,double,damage_scan_count);

          damage_tent_step=*istep;
          damage_tent_increment=iinc;
          damage_tent_step_time=theta**tper;
          damage_tent_total_time=*ttime+damage_tent_step_time;
          damage_tent_count=0;

          for(i=0;i<ne0;i++){
            if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){

              damage_tent_elem[damage_tent_count]=i+1;

              /* For the present non-composite Zr/ZrH C3D4 model this is
                 exactly the same mapping as imat=ielmat(1,i) in
                 calcdamage.f. */
              damage_tent_mat[damage_tent_count]=ielmat[mi[2]*i];

              damage_nip_local=damage_history_nip(&lakon[8*i],mi[0]);
              if(damage_nip_local<1) damage_nip_local=1;
              if(damage_nip_local>mi[0]) damage_nip_local=mi[0];

              damage_dmax=dam[mi[0]*i];
              damage_tent_ip[damage_tent_count]=1;
              for(j=1;j<damage_nip_local;j++){
                if(dam[mi[0]*i+j]>damage_dmax){
                  damage_dmax=dam[mi[0]*i+j];
                  damage_tent_ip[damage_tent_count]=j+1;
                }
              }
              imat=damage_tent_mat[damage_tent_count];
              if((imat>0)&&
                 damage_progressive_material(imat,ndmcon,dmcon,*ndmat_,*ntmat_)&&
                 (damage_dmax>1.)) damage_dmax-=1.;
              if((damage_de13_transaction)&&(damage_de13_trigger_value!=NULL)&&
                 (damage_de13_trigger_value[i]>=0.)){
                damage_tent_value[damage_tent_count]=
                    damage_de13_trigger_value[i];
                damage_tent_ip[damage_tent_count]=
                    damage_de13_trigger_ip[i];
              }else{
                damage_tent_value[damage_tent_count]=damage_dmax;
              }
              damage_tent_count++;
            }
          }
        }

        /* First topology change at this physical load level. */
        damage_active_pass=1;

	iitsav=iit;
	iit=-2;

	remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		   &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
		   ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
		   &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
		   &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
		   &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
		   mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
		   itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
		   tieset,ntie,&num_cpus,ielmat,matname);

	iit=iitsav;

	SFREE(nactdofinv);
	NNEW(nactdofinv,ITG,mt**nk);
	MNEW(nodorig,ITG,*nk);

	FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			       ipkon,lakon,kon,ne));

	SFREE(nodorig);

	ITGMEMSET(iponoel,0,*nk,0);

	FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));

	theta=thetadamage;
	dtheta=dthetadamage;
	dthetaref=dthetarefdamage;
	idiscon=1;
	idamagereeq=1;

	continue;
      }
      }
    }

damage_controller_done:

    /* DE1.3 terminal failure manager.

       DE1.2 has already converged the continuous damage state inside Newton.
       Only now, at a converged physical load level, are near-fully degraded
       elements allowed to change topology.  The batch is tentative: after
       remastruct the model must converge again at exactly the same load.
       Failure enters the existing idamagereeq rollback path and restores the
       topology/damage state from the start of the physical increment. */
    /* Dissipation measure for path following.

       The loading here is displacement controlled, u_p = lambda*u_hat, so
       the Gutierrez dissipation collapses to a scalar built from the load
       factor and the reaction conjugate to the prescribed pattern:

           work      = integral P dlambda ~ (P_n+P_n+1)(l_n+1-l_n)/2
           stored    = (P_n+1 l_n+1 - P_n l_n)/2
           dissipated= (P_n l_n+1 - P_n+1 l_n)/2

       P = sum over the prescribed dofs of fn*xboun, because
       xbounact = lambda*xboun and fn is the reaction there.  resultsini.c
       sets calcul_fn=1 for the NLGEOM Newton path, so fn is current on
       every iteration and this costs one pass over nboun.

       This is the control variable a dissipation-controlled step would
       prescribe.  Measuring it first, before making it the control, keeps
       the two questions separate: is the quantity well behaved, and does
       driving it fix the limit point. */

    if((icutb==0)&&(idamagereeq==0)&&(*nmethod!=4)&&(damage_diss_report==1)){
      if(damage_diss_init==1){
        damage_diss_lamnow=(damage_arc==1)?damage_arc_lam:theta;
        damage_diss_dg=0.5*(damage_diss_pprev*damage_diss_lamnow-
                            damage_diss_p*damage_diss_lprev);
        damage_diss_total+=damage_diss_dg;

        if((damage_diss_ctrl==2)&&(damage_diss_engaged==0)&&
           (damage_diss_dg>DAMAGE_DISS_ENGAGE*damage_diss_target)&&
           ((damage_diss_engage_t<0.)||(theta>=damage_diss_engage_t))){
          damage_diss_engaged=1;
          printf("[DISSIPATION CONTROL] engaged at inc=%" ITGFORMAT
                 " lambda=%.6f dG=%.6e\n",iinc,theta,damage_diss_dg);
          fflush(stdout);
        }

        /* Dissipation-based step control.

           The measured history says the controller walks at full dtmax
           right up to the failure and never pre-emptively cuts: it has no
           way of knowing the dissipation is climbing, because the
           iteration count stays low until the step that cannot converge
           at all.  Sizing the next increment so its predicted dissipation
           stays near a target makes the step shrink as the crack starts
           to run, which is exactly the information the stock controller
           lacks.

           This is NOT path following: lambda still only increases, so a
           genuine snap-back remains out of reach.  It is the cheap half
           of the idea, worth measuring before rebuilding the Newton loop
           around a bordered system. */

        if((damage_diss_target>0.)&&(damage_diss_dg>1.e-30)&&
           (dtheta>0.)&&(*idrct==0)&&(damage_diss_step==1)){
          damage_diss_scale=damage_diss_target/damage_diss_dg;
          if(damage_diss_scale>DAMAGE_DISS_GROW)
            damage_diss_scale=DAMAGE_DISS_GROW;
          if(damage_diss_scale<DAMAGE_DISS_SHRINK)
            damage_diss_scale=DAMAGE_DISS_SHRINK;
          damage_diss_dtheta=dtheta*damage_diss_scale;
          if(damage_diss_dtheta>dthetaref) damage_diss_dtheta=dthetaref;
          if(damage_diss_dtheta<(*tmin)) damage_diss_dtheta=*tmin;
          if(damage_diss_dtheta<0.98*dtheta){
            printf("[DISSIPATION STEP] inc=%" ITGFORMAT
                   " dG=%.6e target=%.6e dtheta %.6e -> %.6e\n",
                   iinc,damage_diss_dg,damage_diss_target,dtheta,
                   damage_diss_dtheta);
            fflush(stdout);
          }
          dtheta=damage_diss_dtheta;
        }

        printf("[DISSIPATION] inc=%" ITGFORMAT " lambda=%.6f P=%.6e "
               "dG=%.6e G=%.6e\n",iinc,theta,damage_diss_p,
               damage_diss_dg,damage_diss_total);
        fflush(stdout);
      }
      damage_diss_init=1;
      damage_diss_lprev=damage_diss_lamnow;
      damage_diss_pprev=damage_diss_p;
    }

    if((icutb==0)&&(idamagereeq==0)){
      damage_path_retry=0;
      damage_path_lamcom=damage_path_lam;
      if(damage_path_desc>0) damage_path_desc--;
    /* Arm only when the step controller has actually run out of
       room.  Arming on icutb>=3 fired at increment 385 on an
       ordinary cutback the stock logic recovers from, derailed the
       path, and never reached the real stall at all. */
    }else if((damage_path_on>=2)&&(dtheta<10.*(*tmin))&&
             (damage_path_desc==0)){
      damage_path_desc=damage_path_nstep;
      damage_path_used=1;
      printf("[DAMAGE PATH] inc=%" ITGFORMAT " arming a %"
             ITGFORMAT "-increment descent at lambda=%.6f\n",
             iinc,damage_path_nstep,damage_path_lamcom);
      fflush(stdout);
    }
    if((damage_de12_enabled)&&(icutb==0)&&(idamagereeq==0)&&
       (damdamageini!=NULL)){

      damage_de13_batch_dmax=0.;
      damage_de13_new=damage_de13_mark_terminal(
          dam,ipkon,lakon,ielmat,mi[2],ndmcon,dmcon,*ndmat_,*ntmat_,
            ne0,mi[0],damage_de13_delete_d,DAMAGE_DE13_BATCH_MAX,
          &damage_de13_batch_dmax,damage_de13_trigger_value,
          damage_de13_trigger_ip,matname,damage_delete_filter,
          damage_damvisc,damage_delete_visc,&damage_de13_batch_vmin);
      damage_de13_term_only=damage_de13_new;

      if((damage_deadall_g>0.)&&(damage_de13_new<DAMAGE_DE13_BATCH_MAX)){
        ITG nda=damage_de13_mark_deadall(
            dam,damage_damvisc,damage_delete_visc,ipkon,lakon,kon,*nk,
            *ne,ne0,mi[0],damage_deadall_g,
            DAMAGE_DE13_BATCH_MAX-damage_de13_new,&damage_deadall_nodes);
        if(nda>0){
          damage_de13_new+=nda;
          damage_deadall_total+=nda;
          printf("[DAMAGE DEADALL] inc=%" ITGFORMAT " time=%.12e "
                 "nodes=%" ITGFORMAT " deleted=%" ITGFORMAT " total=%"
                 ITGFORMAT "\n",iinc,theta**tper,damage_deadall_nodes,
                 nda,damage_deadall_total);
          fflush(stdout);
        }
      }
      if((damage_deadsole_g>0.)&&(damage_de13_new<DAMAGE_DE13_BATCH_MAX)){
        ITG nds=damage_de13_mark_deadsole(
            dam,damage_damvisc,damage_delete_visc,ipkon,lakon,kon,*nk,
            ne0,mi[0],damage_deadsole_g,
            DAMAGE_DE13_BATCH_MAX-damage_de13_new);
        if(nds>0){
          damage_de13_new+=nds;
          damage_deadsole_total+=nds;
          printf("[DAMAGE DEADSOLE] inc=%" ITGFORMAT " time=%.12e "
                 "deleted=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                 iinc,theta**tper,nds,damage_deadsole_total);
          fflush(stdout);
        }
      }

      if(damage_de13_new>0){
        damage_de13_transaction=1;
        damage_soft_reeq=0;
        damage_active_pass=1;
        damage_reeq_uam_ref[0]=uam[0];
        damage_reeq_uam_ref[1]=uam[1];
        /* J-09: keep the reference from collapsing with the step.  The peak
           is a running maximum over the whole step, never reset per
           increment - the same shape as damage_qam_peak. */
        if(damage_reeq_uam_floor>0.){
          if(damage_reeq_uam_ref[0]>damage_reeq_uam_peak[0])
            damage_reeq_uam_peak[0]=damage_reeq_uam_ref[0];
          if(damage_reeq_uam_ref[1]>damage_reeq_uam_peak[1])
            damage_reeq_uam_peak[1]=damage_reeq_uam_ref[1];
          if(damage_reeq_uam_ref[0]<damage_reeq_uam_floor*damage_reeq_uam_peak[0])
            damage_reeq_uam_ref[0]=damage_reeq_uam_floor*damage_reeq_uam_peak[0];
          if(damage_reeq_uam_ref[1]<damage_reeq_uam_floor*damage_reeq_uam_peak[1])
            damage_reeq_uam_ref[1]=damage_reeq_uam_floor*damage_reeq_uam_peak[1];
        }

        if(damage_tent_elem!=NULL){
          SFREE(damage_tent_elem); damage_tent_elem=NULL;
          SFREE(damage_tent_mat); damage_tent_mat=NULL;
          SFREE(damage_tent_ip); damage_tent_ip=NULL;
          SFREE(damage_tent_value); damage_tent_value=NULL;
          damage_tent_count=0;
        }

        /* Detached-island sweep.  BK4 only asks whether the nodes of the
           elements just deleted lost all their elements; a region that is
           internally connected but no longer attached to any support
           passes that test and leaves the operator singular in its
           rigid-body modes.  Marking those elements terminal here, before
           the tentative batch is collected, lets them travel through the
           unchanged transactional path. */

        if(damage_de12_enabled){
          damage_float_new=0;
          damage_float_reach=0;
          damage_float_coh=0;
          damage_float_isl=0;
          FORTRAN(damfloat,(ipkon,kon,lakon,ne,&ne0,nk,nodeboun,nboun,
                            ipompc,nodempc,nmpc,&damage_float_batch,
                            &damage_float_new,&damage_float_reach,
                            &damage_float_coh,&damage_float_isl));
          if(damage_float_new>0){
            damage_float_total+=damage_float_new;
            /* islands_marked and total keep their old meaning - the sum of
               both removals - so stored logs and every parser of them stay
               readable.  facets/islands are ADDED, never substituted: a
               cohesive facet leaves at gmin*Kn, a bulk island leaves at
               whatever it was carrying, and one number for both hid that. */
            printf("[DAMAGE DETACHED] inc=%" ITGFORMAT " time=%.12e "
                   "islands_marked=%" ITGFORMAT " nodes_with_load_path=%"
                   ITGFORMAT " total=%" ITGFORMAT
                   " facets=%" ITGFORMAT " islands=%" ITGFORMAT "\n",
                   iinc,theta**tper,damage_float_new,damage_float_reach,
                   damage_float_total,damage_float_coh,damage_float_isl);
            fflush(stdout);
          }

          /* Fully failed cohesive facets, removed on the same tentative
             topology and inside the same batch, so they re-equilibrate
             together with the bulk deletions instead of in a pass of their
             own.  Default OFF; see CCX_DAMAGE_FACET_DELETE above. */
          if(damage_facetdel){
            damage_facetdel_new=0;
            FORTRAN(damfaildead,(ipkon,lakon,ne,xstate,nstate_,mi,
                                 &damage_float_batch,&damage_facetdel_new));
            if(damage_facetdel_new>0){
              damage_facetdel_total+=damage_facetdel_new;
              printf("[DAMAGE FACET DEAD] inc=%" ITGFORMAT " time=%.12e "
                     "facets_removed=%" ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_facetdel_new,
                     damage_facetdel_total);
              fflush(stdout);
            }
          }

          /* How far the assembled nodal stiffness has fallen.  Reported
             before anything acts on it: the point is to find out whether
             the node the solver actually throws is separable by this
             measure, which mesh topology and geometry both failed to
             do. */
          if((damage_stiff_probe>0)&&(damage_addiag!=NULL)){
            damage_stiff_n1=0;damage_stiff_n2=0;damage_stiff_n3=0;
            damage_stiff_nneg=0;
            damage_stiff_worst=-1;damage_addmin=2.;
            for(i=0;i<*nk;i++){
              if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
              if(damage_addiag0[i]<=0.) continue;
              if(damage_addiag[i]<=0.){damage_stiff_nneg++;continue;}
              damage_addrat=damage_addiag[i]/damage_addiag0[i];
              if(damage_addrat<1.e-3) damage_stiff_n1++;
              if(damage_addrat<1.e-2) damage_stiff_n2++;
              if(damage_addrat<1.e-1) damage_stiff_n3++;
              if(damage_addrat<damage_addmin){
                damage_addmin=damage_addrat;damage_stiff_worst=i+1;
              }
            }
            printf("[DAMAGE STIFFNESS] inc=%" ITGFORMAT " time=%.12e "
                   "below_1e-3=%" ITGFORMAT " below_1e-2=%" ITGFORMAT
                   " below_1e-1=%" ITGFORMAT " nonpositive=%" ITGFORMAT
                   " worst=node_%" ITGFORMAT "_at_%.4e\n",
                   iinc,theta**tper,damage_stiff_n1,damage_stiff_n2,
                   damage_stiff_n3,damage_stiff_nneg,
                   damage_stiff_worst,damage_addmin);
            fflush(stdout);
          }

          /* Material left hanging on a single element.  BK4 tests for
             zero elements and the sweep above tests for a lost load
             path; a node that falls from 22 elements to 1 passes both,
             and the next solve throws it most of a specimen length.
             See damdangle.f for the measurement this comes from. */
          if(damage_bare>damage_bare_rep){
            damage_bare_rep=damage_bare;
            printf("[DAMAGE BARE NODE] inc=%" ITGFORMAT " time=%.12e "
                   "%" ITGFORMAT " node(s) hold only cohesive facets, "
                   "every bulk element gone (new maximum)\n",
                   iinc,theta**tper,damage_bare);
            fflush(stdout);
          }
          /* A node with no bulk left is only in trouble if the facets
             still holding it are themselves debonded.  Measured: DHC1
             carries 30 bare nodes from t=0.209 and runs 300 more
             increments, and the UC6 disk completes with 63 of them, so
             bare alone is not the discriminator.  g = max(gmin,1-dvisc)
             from xstate slot 2 says whether a facet is a real support. */
          if(damage_free_probe==1){
            NNEW(damage_free_nb,ITG,*nk);
            NNEW(damage_free_g,double,*nk);
            for(i=0;i<*nk;i++){damage_free_nb[i]=0;damage_free_g[i]=-1.;}
            for(i=0;i<*ne;i++){
              if(ipkon[i]<0) continue;
              if(lakon[8*i]=='C'){
                for(j=0;j<4;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k>=0)&&(k<*nk)) damage_free_nb[k]++;
                }
              }else if(lakon[8*i]=='U'){
                damage_free_gm=0.;
                for(j=0;j<3;j++){
                  damage_free_dv=xstate[*nstate_*(mi[0]*i+j)+1];
                  if(1.-damage_free_dv>damage_free_gm)
                    damage_free_gm=1.-damage_free_dv;
                }
                for(j=0;j<6;j++){
                  k=kon[ipkon[i]+j]-1;
                  if((k<0)||(k>=*nk)) continue;
                  if(damage_free_gm>damage_free_g[k])
                    damage_free_g[k]=damage_free_gm;
                }
              }
            }
            damage_free_cnt=0;
            for(i=0;i<*nk;i++){
              if(damage_free_nb[i]>0) continue;
              if(damage_free_g[i]<0.) continue;
              if(damage_free_g[i]<1.e-2) damage_free_cnt++;
            }
            if(damage_free_cnt>damage_free_rep){
              damage_free_rep=damage_free_cnt;
              printf("[DAMAGE FREE NODE] inc=%" ITGFORMAT " time=%.12e "
                     "%" ITGFORMAT " node(s) with no bulk AND every "
                     "remaining facet debonded (new maximum)\n",
                     iinc,theta**tper,damage_free_cnt);
              fflush(stdout);
            }
            SFREE(damage_free_nb);SFREE(damage_free_g);
          }
          if((damage_dangle_max>0)||(damage_stiff_min>0.)){
            damage_dangle_new=0;
            damage_dangle_weak=0;
            if(damage_stiff_haz==NULL) NNEW(damage_stiff_haz,ITG,*nk);
            for(i=0;i<*nk;i++) damage_stiff_haz[i]=0;
            if((damage_stiff_min>0.)&&(damage_addiag!=NULL)){
              for(i=0;i<*nk;i++){
                if((damage_addok!=NULL)&&(damage_addok[i]==0)) continue;
                if(damage_addiag0[i]<=0.) continue;
                if(damage_addiag[i]<=0.) continue;
                if(damage_addiag[i]<damage_stiff_min*damage_addiag0[i])
                  damage_stiff_haz[i]=1;
              }
            }
            FORTRAN(damdangle,(ipkon,kon,lakon,ne,nk,
                               &damage_float_batch,&damage_dangle_max,
                               damage_stiff_haz,
                               &damage_dangle_new,&damage_dangle_weak,
                               &damage_bare));
            if(damage_dangle_new>0){
              damage_dangle_total+=damage_dangle_new;
              printf("[DAMAGE DANGLING] inc=%" ITGFORMAT " time=%.12e "
                     "elements_marked=%" ITGFORMAT " nodes_stranded=%"
                     ITGFORMAT " total=%" ITGFORMAT "\n",
                     iinc,theta**tper,damage_dangle_new,
                     damage_dangle_weak,damage_dangle_total);
              fflush(stdout);
            }
          }
        }

        damage_scan_count=0;
        for(i=0;i<ne0;i++){
          if((ipkondamageini[i]>=0)&&(ipkon[i]<0)) damage_scan_count++;
        }

        if(damage_scan_count>0){
          NNEW(damage_tent_elem,ITG,damage_scan_count);
          NNEW(damage_tent_mat,ITG,damage_scan_count);
          NNEW(damage_tent_ip,ITG,damage_scan_count);
          NNEW(damage_tent_value,double,damage_scan_count);

          damage_tent_step=*istep;
          damage_tent_increment=iinc;
          damage_tent_step_time=theta**tper;
          damage_tent_total_time=*ttime+damage_tent_step_time;
          damage_tent_count=0;

          for(i=0;i<ne0;i++){
            if((ipkondamageini[i]>=0)&&(ipkon[i]<0)){
              damage_tent_elem[damage_tent_count]=i+1;
              damage_tent_mat[damage_tent_count]=ielmat[mi[2]*i];

              damage_nip_local=damage_history_nip(&lakon[8*i],mi[0]);
              if(damage_nip_local<1) damage_nip_local=1;
              if(damage_nip_local>mi[0]) damage_nip_local=mi[0];

              damage_dmax=dam[mi[0]*i];
              damage_tent_ip[damage_tent_count]=1;
              for(j=1;j<damage_nip_local;j++){
                if(dam[mi[0]*i+j]>damage_dmax){
                  damage_dmax=dam[mi[0]*i+j];
                  damage_tent_ip[damage_tent_count]=j+1;
                }
              }

              imat=damage_tent_mat[damage_tent_count];
              if((imat>0)&&
                 damage_progressive_material(imat,ndmcon,dmcon,*ndmat_,*ntmat_)&&
                 (damage_dmax>1.)) damage_dmax-=1.;
              if((damage_de13_transaction)&&(damage_de13_trigger_value!=NULL)&&
                 (damage_de13_trigger_value[i]>=0.)){
                damage_tent_value[damage_tent_count]=
                    damage_de13_trigger_value[i];
                damage_tent_ip[damage_tent_count]=
                    damage_de13_trigger_ip[i];
              }else{
                damage_tent_value[damage_tent_count]=damage_dmax;
              }
              damage_tent_count++;
            }
          }
        }

        /* The batch as a SET: sorted, hashed and printed in full, next to
           the committed-state hash of the attempt that produced it.  Two
           batches with the same hash are the same set of elements
           whatever order the scan produced them in, and the PAIR
           (committed state, batch) is what would have to repeat for the
           batch/rollback loop to be real. */

        if((td_armed!=0)&&(td_trace!=0)&&(damage_tent_count>0)){
          td_batch=topodiag_hash_batch(damage_tent_elem,damage_tent_count,
                                       td_sort,topodiag_hash_seed());
          printf("[BATCHTRACE] batch inc=%" ITGFORMAT " icutb=%" ITGFORMAT
                 " pass=%" ITGFORMAT " n=%" ITGFORMAT
                 " committed_state=%016llx batch=%016llx sorted:",
                 iinc,icutb,damage_active_pass,damage_tent_count,
                 td_state,td_batch);
          for(k=0;k<damage_tent_count;k++)
            printf(" %" ITGFORMAT,td_sort[k]);
          printf("\n");
          fflush(stdout);
        }

        damage_topology_rebuild=1;
        damage_topology_orphans=0;
        if(damage_topology_deferred_mode==1){
          NNEW(damage_iponoel_trial,ITG,*nk);
          ITGMEMSET(damage_iponoel_trial,0,*nk,0);
          FORTRAN(nodebelongstoel,(damage_iponoel_trial,lakon,ipkon,kon,ne));
          NNEW(damage_orphan_seen,ITG,*nk);
          for(k=0;k<damage_tent_count;k++){
            i=damage_tent_elem[k]-1;
            if((i<0)||(i>=ne0)||(ipkon[i]>=0)) continue;
            if(strncmp(&lakon[8*i],"C3D4",4)!=0) continue;
            damage_indexe=-ipkon[i]-2;
            for(j=0;j<4;j++){
              inode=kon[damage_indexe+j]-1;
              if((inode<0)||(inode>=*nk)) continue;
              if((damage_iponoel_trial[inode]==0)&&
                 (damage_orphan_seen[inode]==0)){
                for(idir=1;idir<=3;idir++){
                  if(nactdof[mt*inode+idir]>0){
                    damage_topology_orphans++;
                    damage_orphan_seen[inode]=1;
                    break;
                  }
                }
              }
            }
          }
          if(damage_topology_orphans==0){
            isiz=*nk;
            cpyparitg(iponoel,damage_iponoel_trial,&isiz,&num_cpus);
            damage_topology_rebuild=0;
          }
          SFREE(damage_orphan_seen);
          SFREE(damage_iponoel_trial);
        }

        printf("[DAMAGE DE1.3 TERMINAL] inc=%" ITGFORMAT
               " time=%.12e new_terminal=%" ITGFORMAT
               " tentative=%" ITGFORMAT " batch_Dmax=%.6e "
               "batch_Dvis=%.6e "
               "orphan_nodes=%" ITGFORMAT " action=%s+same-load-Newton\n",
               iinc,theta**tper,damage_de13_new,damage_tent_count,
               damage_de13_batch_dmax,damage_de13_batch_vmin,
               damage_topology_orphans,
               damage_topology_rebuild?"remastruct":"reuse-sparse-graph");
        fflush(stdout);

        /* [DAMAGE RELEASE] arm.  f still holds f_int(u*) of the converged
           state - the elements were marked above but no results() has run
           since - and nactdof is still the pre-remastruct numbering.  Both
           are mapped into node space here so they survive remastruct. */
        if(damage_release_probe){
          ITG ai,aj,ak;
          if(damage_frel==NULL){
            NNEW(damage_frel,double,mt**nk);
            NNEW(damage_ract,ITG,mt**nk);
          }
          for(ai=0;ai<*nk;ai++){
            for(aj=0;aj<mt;aj++){
              ak=nactdof[mt*ai+aj];
              damage_frel[mt*ai+aj]=(ak>0)?f[ak-1]:0.;
              damage_ract[mt*ai+aj]=(ak>0)?1:0;
            }
          }
          damage_release_armed=1;
          damage_release_pass++;
          damage_release_qa=qa[0];
          damage_release_qam=qam[0];
          damage_release_dt=dtime;
          damage_release_rebuild=damage_topology_rebuild;
          damage_release_iforbou=iforbou;
          damage_release_nterm=damage_de13_term_only;
          damage_release_nother=damage_de13_new-damage_de13_term_only;
          damage_release_nisl=damage_float_isl;
          damage_release_ncoh=damage_float_coh;

          /* level 2: one line per element in this batch.  The degradation it
             carried when it left is printed as a REFERENCE CHARACTERISTIC -
             it is a dimensionless multiplier, not a force, and it is NOT
             summed into any estimate of the released force.  The measured
             force is the surv/removed split above.

             It does classify the source for free, though: the terminal
             trigger cannot fire below its own threshold, so an element that
             left at Dvis < delete_d did NOT come from it - it came from
             damfloat, which reads no damage variable at all and therefore
             removes at whatever the element was carrying. */
          if(damage_release_probe>=2){
            ITG pe,pj,pnip,pel;
            double pd,pdv,pdmax,pdvmax,ptrig;
            for(pe=0;pe<damage_tent_count;pe++){
              pel=damage_tent_elem[pe]-1;
              if((pel<0)||(pel>=ne0)) continue;
              pnip=damage_history_nip(&lakon[8*pel],mi[0]);
              if(pnip<1) pnip=1;
              if(pnip>mi[0]) pnip=mi[0];
              pdmax=0.; pdvmax=0.;
              for(pj=0;pj<pnip;pj++){
                pd=dam[mi[0]*pel+pj]-1.;
                if(pd<0.) pd=0.;
                if(pd>1.) pd=1.;
                if(pd>pdmax) pdmax=pd;
                if(damage_damvisc!=NULL){
                  pdv=damage_damvisc[mi[0]*pel+pj];
                  if(pdv<0.) pdv=0.;
                  if(pdv>1.) pdv=1.;
                  if(pdv>pdvmax) pdvmax=pdv;
                }
              }
              ptrig=((damage_delete_visc==1)&&(damage_damvisc!=NULL))
                    ?pdvmax:pdmax;
              /* The bulk damage variable exists ONLY for bulk elements.
                 calcdamage.f:135 skips every lakon(1:1) != 'C', so dam is
                 identically zero for a UC6 facet - printing it would read
                 as "left fully intact" when in truth the facet was
                 conducting gmin*Kn and its state lives in xstate, not here.
                 Both decks in play carry UC6, so this is not hypothetical:
                 damfloatcoh removals land in the same tentative batch. */
              if(lakon[8*pel]!='C'){
                printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                       " pass=%" ITGFORMAT " el=%" ITGFORMAT
                       " mat=%" ITGFORMAT " type=%.8s"
                       " D=n/a Dvis=n/a g_ref=n/a"
                       " note=non-bulk-element-damage-lives-in-xstate%s",
                       iinc,damage_release_pass,damage_tent_elem[pe],
                       damage_tent_mat[pe],&lakon[8*pel],
                       "\n");
                continue;
              }
              printf("[DAMAGE RELEASE ELEM] inc=%" ITGFORMAT
                     " pass=%" ITGFORMAT " el=%" ITGFORMAT
                     " mat=%" ITGFORMAT " type=%.8s"
                     " D=%.6f Dvis=%.6f g_ref=%.4e"
                     " below_delete_d=%s%s",
                     iinc,damage_release_pass,damage_tent_elem[pe],
                     damage_tent_mat[pe],&lakon[8*pel],pdmax,pdvmax,1.-ptrig,
                     (ptrig<damage_de13_delete_d)?"YES-not-terminal":"no",
                     "\n");
            }
            fflush(stdout);
          }
        }

        if(damage_topology_rebuild){
          iitsav=iit;
          iit=-2;

          remastruct(ipompc,&coefmpc,&nodempc,nmpc,
                     &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
                     ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
                     &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
                     &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
                     &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
                     mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
                     itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
                     tieset,ntie,&num_cpus,ielmat,matname);

          iit=iitsav;

          SFREE(nactdofinv);
          NNEW(nactdofinv,ITG,mt**nk);
          MNEW(nodorig,ITG,*nk);

          FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
                                 ipkon,lakon,kon,ne));

          SFREE(nodorig);

          ITGMEMSET(iponoel,0,*nk,0);
          FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));
        }

        theta=thetadamage;
        dtheta=dthetadamage;
        dthetaref=dthetarefdamage;
        idiscon=1;
        idamagereeq=1;
        continue;
      }
    }

    /* DE1.2 has already produced the converged trial damage inside the last
       Newton iteration.  There is no outer same-load damage fixed point to
       close.  Once the physical increment is accepted, report/serialize the
       exact integration-point state directly. */
    if((damage_de12_enabled)&&(icutb==0)&&(idamagereeq==0)&&
       (damdamageini!=NULL)){
      damage_de1_stats(dam,damdamageini,ipkon,lakon,ne0,mi[0],
                       &damage_de1_nactive,&damage_de1_gt01,
                       &damage_de1_gt05,&damage_de1_gt09,
                       &damage_de1_nfull,&damage_de1_nchanged,
                       &damage_de1_dmax,&damage_de1_maxdelta);

      if(damage_de1_nactive>0){
        printf("[DAMAGE DE1.2 COMMIT] inc=%" ITGFORMAT
               " time=%.12e active=%" ITGFORMAT
               " D>0.1=%" ITGFORMAT " D>0.5=%" ITGFORMAT
               " D>0.9=%" ITGFORMAT " Dfull=%" ITGFORMAT
               " Dmax=%.6e max_dD_inc=%.6e\n",
               iinc,theta**tper,damage_de1_nactive,damage_de1_gt01,
               damage_de1_gt05,damage_de1_gt09,damage_de1_nfull,
               damage_de1_dmax,damage_de1_maxdelta);

        damage_de1_append_stats(jobnamec,*istep,iinc,theta**tper,
                                *ttime+theta**tper,1,damage_de1_nactive,
                                damage_de1_gt01,damage_de1_gt05,
                                damage_de1_gt09,damage_de1_nfull,
                                damage_de1_dmax,damage_de1_maxdelta);

        damage_de1_write_vtk(jobnamec,co,vold,*nk,mt,kon,ipkon,lakon,
                             ielmat,mi[2],dam,mi[0],ne0,*istep,iinc,
                             theta**tper);
        fflush(stdout);
      }
    }

    /* Local damage substepping.  checkconvergence() is left in charge
       of the actual cutback.  Here we only keep the original end point
       of the failed coarse increment and prevent the following successful
       reduced increments from stepping beyond it.  A physical substep is
       considered accepted only after any damage re-equilibration at its
       end has also converged. */

    if((icutb==0)&&(idamagereeq==0)&&(ilocalsubstep==1)){
      dtheta_remaining=theta_goal-theta;

      if((damage_fast_recover==1)&&(dtheta_remaining>1.e-12)&&
         (*itpamp==0)){

        /* Jump directly toward the saved coarse-increment goal instead of
           spending many accepted increments regrowing a tiny event step.
           This is only a proposed next step: Newton/checkconvergence can
           cut it back normally if the jump is too aggressive. */

        dtheta=dtheta_remaining;
        if(dtheta>dtheta_restore) dtheta=dtheta_restore;
        if(dtheta>*tmax) dtheta=*tmax;
        if(dtheta>1.-theta) dtheta=1.-theta;
        dthetaref=dtheta;
        istab=0;
        damage_fast_recover=0;

        printf("[DAMAGE FAST RECOVER] time=%e goal=%e next_dt=%e\n\n",
               theta**tper,theta_goal**tper,dtheta**tper);

      }else if(dtheta_remaining<=1.e-12){

        /* Eliminate round-off drift at the saved target. */

        theta=theta_goal;
        ilocalsubstep=0;

        /* Do not let the two-success CCX growth history accumulated by
           the internal substeps influence the first normal increment. */

        istab=0;

        printf(" Local damage substepping complete at step time %e.\n",
               theta**tper);

        /* TIME POINTS are handled inside checkconvergence().  Replacing
           its selected dtheta here could skip a requested time point, so
           only restore the pre-cutback coarse increment when no TIME
           POINTS sequence is active. */

        if((*itpamp==0)&&((1.-theta)>1.e-12)){
          dtheta=dtheta_restore;
          if(dtheta>*tmax) dtheta=*tmax;
          if(dtheta>1.-theta) dtheta=1.-theta;
          dthetaref=dtheta;
          printf(" Restoring normal increment size to %e.\n\n",
                 dtheta**tper);
        }else{
          printf(" Returning increment-size control to CalculiX.\n\n");
        }

      }else if(dtheta>dtheta_remaining){

        dtheta=dtheta_remaining;
        dthetaref=dtheta;
        printf(" Local damage substepping: next increment is clipped to "
               "%e to reach saved goal %e.\n\n",
               dtheta**tper,theta_goal**tper);
      }
    }

    /* printing the energies (only for dynamic calculations) */

    if((icutb==0)&&(*nmethod==4)&&(*ithermal<2)&&(jout[0]==jprint)&&
       (*nener==1)){

      printenergy(iexpl,ttime,&theta,tper,energy,ne,nslavs,ener,&energyref,
		  &allwk,&dampwk,&ea,&energym,&energymold,&jnz,&mscalmethod,
		  mortar,mi);

    }

    if(uncoupled){
      SFREE(iruc);
    }

    if(((qa[0]>ea*qam[0])||(qa[1]>ea*qam[1]))&&(icutb==0)){jnz++;}
    iit=0;

    if(icutb!=0){

      if(damage_rescue_bt_on==1){
        printf("[DAMAGE RESCUE] STANDARD cutback rollback executed: vold, "
               "xbounact, f, sti, eme, ener, xstate, dam, damvisc and ipkon "
               "restored from the increment-start baselines; no additional "
               "snapshot was taken%s","\n");
        fflush(stdout);
      }

      /* A cutback may come either from the ordinary physical Newton
         solve or from the same-load damage re-equilibration.  In both
         cases checkconvergence() has already reduced dtheta.  Arm local
         substepping once and keep the target established at the start of
         the failed coarse increment.  The existing state/topology rollback
         below remains unchanged. */

      if((ilocalsubstep==0)&&(*ndmat_>0)&&(*iexpl<=1)&&
         (*nmethod!=4)&&(*idrct==0)&&
         (theta_goal>theta+1.e-12)){
        ilocalsubstep=1;
        if(damage_event_cut==1){
          printf(" Local damage substepping activated by damage-event "
                 "localization.\n");
        }else if(idamagereeq==1){
          printf(" Local damage substepping activated after damage "
                 "re-equilibration divergence.\n");
        }else{
          printf(" Local damage substepping activated after physical "
                 "Newton divergence.\n");
        }
        printf(" Saved interval: %e -> %e; reduced retry increment: %e.\n\n",
               theta_local_start**tper,theta_goal**tper,dtheta**tper);
      }

      isiz=mt**nk;cpypardou(vold,vini,&isiz,&num_cpus);

      isiz=*nboun;cpypardou(xbounact,xbounini,&isiz,&num_cpus);
      if((*ithermal==1)||(*ithermal>=3)){
	isiz=*nk;cpypardou(t1act,t1ini,&isiz,&num_cpus);
      }
      isiz=neq[1];cpypardou(f,fini,&isiz,&num_cpus);
      if(*nmethod==4){
	isiz=mt**nk;
	cpypardou(veold,veini,&isiz,&num_cpus);
	cpypardou(accold,accini,&isiz,&num_cpus);
	isiz=neq[1];
	cpypardou(fext,fextini,&isiz,&num_cpus);
	cpypardou(cv,cvini,&isiz,&num_cpus);
	if(*ithermal<2){
	  allwk=allwkini;
	  if(idamping==1)dampwk=dampwkini;
	  for(k=0;k<4;k++){
	    energy[k]=energyini[k];
	  }
	}
      }
      if(*ithermal!=2){
	isiz=6*mi[0]*ne0;
	cpypardou(sti,stiini,&isiz,&num_cpus);
	cpypardou(eme,emeini,&isiz,&num_cpus);
      }
      if(*nener==1){
	isiz=2*mi[0]*ne0;
	cpypardou(ener,enerini,&isiz,&num_cpus);
      }

      isiz=*nstate_*mi[0]*(ne0+maxprevcontel);cpypardou(xstate,xstateini,
							&isiz,&num_cpus);

      qam[0]=qamold[0];
      qam[1]=qamold[1];


      /* the rollback restores qam from before the failed attempt, which can
         be under the floor again; re-apply it so the criterion the retry
         faces is the same one the attempt faced (CCX_DAMAGE_QAM_FLOOR) */
      if((damage_qam_floor>0.)&&(qam[0]<damage_qam_floor*damage_qam_peak)){
        qam[0]=damage_qam_floor*damage_qam_peak;
      }

      if(*mortar>1){
	for (i=0;i<*ntie;i++){
	  for(j=nslavnode[i];j<nslavnode[i+1];j++){
	    islavact[j]=islavactini[j];
	    bp[j]=bpini[j];
	    for(k=0;k<mt;k++){
	      cstress[mt*j+k]=cstressini[mt*j+k];
	    }
	  }    
	} 
      }
      /* if the failed attempt was a damage re-equilibration,
         restore the topology and damage state at the beginning
         of the physical increment before retrying with the smaller step */

      if(idamagereeq==1){

        if(damage_soft_reeq==1){
          printf("[DAMAGE DE1 ROLLBACK] inc=%" ITGFORMAT
                 " time=%.12e -> restore constitutive baseline\n",
                 iinc,theta**tper);
        }else if(damage_de13_transaction){
          printf("[DAMAGE DE1.3 ROLLBACK] inc=%" ITGFORMAT
                 " time=%.12e tentative_terminal=%" ITGFORMAT
                 " -> restore topology and constitutive baseline\n",
                 damage_tent_increment,damage_tent_step_time,
                 damage_tent_count);
        }else{
          printf("[DAMAGE ROLLBACK] inc=%" ITGFORMAT
                 " time=%.12e tentative=%" ITGFORMAT "\n",
                 damage_tent_increment,damage_tent_step_time,
                 damage_tent_count);
        }
        fflush(stdout);

        /* Trial deletions must never enter jobname.damage. */
        if(damage_tent_elem!=NULL){
          SFREE(damage_tent_elem); damage_tent_elem=NULL;
          SFREE(damage_tent_mat); damage_tent_mat=NULL;
          SFREE(damage_tent_ip); damage_tent_ip=NULL;
          SFREE(damage_tent_value); damage_tent_value=NULL;
        }
        damage_tent_count=0;
        if(damage_de13_transaction) damage_path_retry++;
        damage_de13_transaction=0;
        if(damage_de13_trigger_value!=NULL){
          for(i=0;i<ne0;i++){
            damage_de13_trigger_value[i]=-1.;
            damage_de13_trigger_ip[i]=0;
          }
        }

        /* A3 safety latch: if this failed transaction came from an
           aggressive fast trial, the retry uses the original exact A2
           event locator.  The latch survives the current increment retry
           and is reset only after that physical increment is accepted. */
        if(damage_fast_used==1){
          damage_fast_retry=1;
          damage_fast_used=0;
          damage_fast_recover=0;
          damage_fast_failures++;
          printf("[DAMAGE FAST BACKOFF] inc=%" ITGFORMAT
                 " failures=%" ITGFORMAT
                 " -> exact event localization on retry\n",
                 iinc,damage_fast_failures);
          fflush(stdout);
        }

	isiz=ne0;
	cpyparitg(ipkon,ipkondamageini,&isiz,&num_cpus);
	isiz=mi[0]*ne0;
	cpypardou(dam,damdamageini,&isiz,&num_cpus);

	if(*nmethod!=4){
	  isiz=mt**nk;
	  cpypardou(veold,veolddamageini,&isiz,&num_cpus);
	}

        if(damage_soft_reeq==0){
	iitsav=iit;
	iit=-2;

	remastruct(ipompc,&coefmpc,&nodempc,nmpc,
		   &mpcfree,nodeboun,ndirboun,nboun,ikmpc,ilmpc,
		   ikboun,ilboun,labmpc,nk,&memmpc_,&icascade,
		   &maxlenmpc,kon,ipkon,lakon,ne,nactdof,icol,jq,
		   &irow,isolver,neq,nzs,nmethod,&f,&fext,&b,&aux2,
		   &fini,&fextini,&adb,&aub,ithermal,iperturb,mass,
		   mi,iexpl,mortar,typeboun,&cv,&cvini,&iit,network,
		   itiefac,&ne0,&nkon0,nintpoint,islavsurf,pmastsurf,
		   tieset,ntie,&num_cpus,ielmat,matname);

	iit=iitsav;

	SFREE(nactdofinv);
	NNEW(nactdofinv,ITG,mt**nk);
	MNEW(nodorig,ITG,*nk);

	FORTRAN(gennactdofinv,(nactdof,nactdofinv,nk,mi,nodorig,
			       ipkon,lakon,kon,ne));

	SFREE(nodorig);

	ITGMEMSET(iponoel,0,*nk,0);

	FORTRAN(nodebelongstoel,(iponoel,lakon,ipkon,kon,ne));
        }

	theta=thetadamage;
	idiscon=0;
	idamagereeq=0;
        damage_soft_reeq=0;

      }

    }
    /* ---- [DAMAGE CT] CLEAN PARTIAL EXIT ---------------------------
       Placed at the COMMON post-rollback join: the point reached once per
       outer pass, immediately after the complete if(icutb!=0) rollback
       block closes.  icutb!=0 here means checkconvergence() has just
       decided a cutback (it sets *icutb=0 on convergence and (*icutb)++
       otherwise), so the rollback above has just restored vold, xbounact,
       f, sti, eme, ener, xstate, dam, damvisc, ipkon and veold from the
       increment-start baselines: the reported state IS the last committed
       one, on EVERY rollback path.
       The earlier placement inside if(idamagereeq==1) was wrong twice
       over: that branch is only one of the rollback paths, so an ordinary
       cutback after a refusal never reached it, and the run ended at the
       stock stop with no PARTIAL report at all.  Placing it in the commit
       path was wrong for the same reason: a refusal never converges.
       It reports every quantity the record needs, states plainly that this
       is NOT a completed step and NOT a restart point, and leaves without
       the normal end-of-step bookkeeping.  Chosen deliberately (logged):
       the exit uses the stock 201 code, which already means "did not
       complete" and therefore cannot be mistaken for COMPLETED; a
       dedicated code would need a change outside nonlingeo.c, which this
       MVP is not allowed to make. */
    if((damage_ct_partial==1)&&(icutb!=0)){
        ITG cq1=0,cq2=0,cq3=0,cqi,cqj,cqk;
        double cqd0,cqdl[3],cqrm[9],cqsh[3],cqop=0.;
        for(cqi=0;cqi<ne0;cqi++){
          if(ipkon[cqi]<0){cq2++;continue;}
          if(lakon[8*cqi]!='U') continue;
          if(ielprop[cqi]<0) continue;
          cqd0=prop[ielprop[cqi]+1];
          if(cqd0<=0.) continue;
          cqd0=cqd0/prop[ielprop[cqi]];
          for(cqj=0;cqj<3;cqj++){
            if(cqj>=mi[0]) break;
            cqk=mi[0]*cqi+cqj;
            if(xstate[*nstate_*cqk+3]>=0.5) cq1++;
            if(xstate[*nstate_*cqk]>cqd0) cq3++;
          }
        }
        if((damage_ct_elem>=0)&&(damage_ct_elem<ne0)&&
           (ipkon[damage_ct_elem]>=0)&&(damage_ct_w!=NULL)){
          damage_ct_kin(co,kon,ipkon[damage_ct_elem],vold,mt,damage_ct_ip,
                        cqdl,cqrm,cqsh);
          cqop=damage_ct_m[0]*(cqdl[0]-damage_ct_dc[0])
              +damage_ct_m[1]*(cqdl[1]-damage_ct_dc[1])
              +damage_ct_m[2]*(cqdl[2]-damage_ct_dc[2]);
        }
        printf("%s","\n");
        printf("[DAMAGE CT] ===== PARTIAL - bounded experimental "
               "continuation segment =====%s","\n");
        printf("[DAMAGE CT]   last COMMITTED state is preserved; this is NOT "
               "a completed CalculiX step and NOT a valid restart point.%s",
               "\n");
        printf("[DAMAGE CT]   theta (pseudo-time) = %.12e   time = %.12e   "
               "ttime = %.12e%s",theta,theta**tper,*ttime,"\n");
        printf("[DAMAGE CT]   lambda = %.12e   lambda_c = %.12e   "
               "ds = %.6e   kappa = %.6e%s",
               damage_ct_lam,damage_ct_lamc,damage_ct_ds,damage_ct_kappa,
               "\n");
        printf("[DAMAGE CT]   control point: element %" ITGFORMAT " ip %"
               ITGFORMAT "   m = (%.6f,%.6f,%.6f)   m.(delta-delta_c) = "
               "%.6e   |c| last = %.6e   tol_c = %.6e%s",
               damage_ct_elem+1,damage_ct_ip+1,damage_ct_m[0],
               damage_ct_m[1],damage_ct_m[2],cqop,fabs(damage_ct_cprev),
               damage_ct_tolc,"\n");
        printf("[DAMAGE CT]   observables: P1 failed facets %" ITGFORMAT
               ", P2 deleted elements %" ITGFORMAT ", P3 ip with dmax>d0 %"
               ITGFORMAT "%s",cq1,cq2,cq3,"\n");
        printf("[DAMAGE CT]   cost: %" ITGFORMAT " accepted continuation "
               "commits, %" ITGFORMAT " corrector iterations, %" ITGFORMAT
               " factorisations, %" ITGFORMAT " residual evaluations, %"
               ITGFORMAT " blacklisted control point(s)%s",
               damage_ct_ncommit,damage_ct_ncorr,damage_ct_nfact,
               damage_ct_neval,damage_ct_nbl,"\n");
        printf("[DAMAGE CT] ===== end PARTIAL report =====%s","\n");
        fflush(stdout);
        FORTRAN(stop,());
      }

    
    /* damage_event_cut is only a label for the rollback just consumed.
       Clear it before the next attempt so a later ordinary Newton
       cutback cannot be misclassified. */
    damage_event_cut=0;

    /* face-to-face penalty */

    if((*mortar==1)&&(icutb==0)&&(ncont!=0)){
	
      ntrimax=0;
      for(i=0;i<*ntie;i++){	    
	if(itietri[2*i+1]-itietri[2*i]+1>ntrimax)		
	  ntrimax=itietri[2*i+1]-itietri[2*i]+1;  	
      }
      MNEW(xo,double,ntrimax);	    
      MNEW(yo,double,ntrimax);	    
      MNEW(zo,double,ntrimax);	    
      MNEW(x,double,ntrimax);	    
      MNEW(y,double,ntrimax);	    
      MNEW(z,double,ntrimax);	   
      MNEW(nx,ITG,ntrimax);	   
      MNEW(ny,ITG,ntrimax);	    
      MNEW(nz,ITG,ntrimax);
      
      /*  Determination of active nodes (islavact) */
      
      FORTRAN(islavactive,(tieset,ntie,itietri,cg,straight,
			   co,vold,xo,yo,zo,x,y,z,nx,ny,nz,mi,
			   imastop,nslavnode,islavnode,islavact));

      SFREE(xo);SFREE(yo);SFREE(zo);SFREE(x);SFREE(y);SFREE(z);SFREE(nx);
      SFREE(ny);SFREE(nz);

      if(*ithermal!=2){
	if(negpres==0){
	  if((*mortar==1)&&(1.-theta-dtheta<=1.e-6)&&(itruecontact==1)){
	    printf(" pressure ratio (smallest/largest pressure over all contact areas) =%e\n\n",pressureratio);
	    	    if(pressureratio<-0.05){
	    //	    if((pressureratio<-0.05)||((*nmethod==1)&&(iperturb[1]==1))){
	      printf(" zero-size increment is appended\n\n");
	      negpres=1;theta=1.-1.e-6;dtheta=1.e-6;
	    }
	  }
	}else{negpres=0;}
      }

    }

    /* output */

    if((jout[0]==jprint)&&(icutb==0)){

      jprint=0;

      /* calculating the displacements and the stresses and storing */
      /* the results in frd format  */
	
      MNEW(v,double,mt**nk);
      MNEW(fn,double,mt**nk);
      NNEW(stn,double,6**nk);
      if(*ithermal>1) NNEW(qfn,double,3**nk);
      NNEW(inum,ITG,*nk);
      NNEW(stx,double,6*mi[0]**ne);
      
      if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
      if(strcmp1(&filab[435],"PEEQ")==0) NNEW(epn,double,*nk);
      if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
      if(strcmp1(&filab[609],"SDV ")==0) NNEW(xstaten,double,*nstate_**nk);
      if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);
      if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
      if(strcmp1(&filab[4785],"DUCT")==0) NNEW(damn,double,*nk);

      isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);

      if((*mortar==-1)&&(idispfrdonly==1)){

	/* nothing to do if massless explicit dynamics and all output
           consists of displacements */
	
	cpyparitg(inum,inumcp,nk,&num_cpus);

      }else{
      
	iout=2;
	icmd=3;
      
#ifdef COMPANY
	FORTRAN(uinit,());
#endif
	results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
		prestr,iprestr,filab,eme,emn,een,iperturb,
		f,fn,nactdof,&iout,qa,vold,b,nodeboun,
		ndirboun,xbounact,nboun,ipompc,
		nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		&bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
		xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
		ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
		&reltime,&ne0,thicke,shcon,nshcon,
		sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
		mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
		inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
		itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
		islavquadel,aut,irowt,jqt,&mortartrafoflag,
		&intscheme,physcon,dam,damn,iponoel);
      
	isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);

      }
      
      iout=0;
      if(*iexpl<=1) icmd=0;
      
      ++*kode;
      if(*mcs!=0){
	ptime=*ttime+time;

	if(*mortar>1){
	  mortar_prefrd(ne,nslavs,mi,nk,nkon,&stx,cdisp,fn,cfs,cfm);       
	}
	
	frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	       t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
               nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
               ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	       norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0,
               cdn,mortar,nmat,qfx,ielprop,prop,damn,&errn);

	if(*mortar>1){
	  mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
	}
#ifdef COMPANY
	FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
      }
      else{
	if(strcmp1(&filab[1044],"ZZS")==0){
	  NNEW(neigh,ITG,40**ne);
	  MNEW(ipneigh,ITG,*nk);
	}

	ptime=*ttime+time;

	if(*mortar>1){
	  mortar_prefrd(ne,nslavs,mi,nk,nkon, &stx,cdisp,fn,cfs,cfm);       
	}
	frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,
	    prop,sti,damn,&errn);
	if(*mortar>1){
	  mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
	}

	if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
#ifdef COMPANY
	FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
      }

      /* mesh refinement */
  
      if(strcmp1(&filab[4089],"RM")==0){
	refinemesh(nk,ne,co,ipkon,kon,v,veold,stn,een,emn,epn,enern,
		   qfn,errn,filab,mi,lakon,jobnamec,istartset,iendset,
		   ialset,set,nset,matname,ithermal,output,nmat,
		   nelemload,nload,sideload,nodeforc,
		   nforc,nodeboun,nboun,nodempc,ipompc,nmpc);

	/* free errn */
	
	if(((*nmethod!=5)||(mode==-1))&&
	   ((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2))) SFREE(errn);
      }
      
      SFREE(v);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(stx);
      if(*ithermal>1){SFREE(qfn);}
      
      if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
      if(strcmp1(&filab[435],"PEEQ")==0) SFREE(epn);
      if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
      if(strcmp1(&filab[609],"SDV ")==0) SFREE(xstaten);
      if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);
      if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
      if(strcmp1(&filab[4785],"DUCT")==0) SFREE(damn);
    }
    
  }

  /*********************************************************/
  /*   end of the increment loop                          */
  /*********************************************************/

  if(jprint!=0){
    
    /* printing the energies (only for dynamic calculations) */

    if((*nmethod==4)&&(*ithermal<2)&&(*nener==1)){

      printenergy(iexpl,ttime,&theta,tper,energy,ne,nslavs,ener,&energyref,
		  &allwk,&dampwk,&ea,&energym,&energymold,&jnz,&mscalmethod,
		  mortar,mi);
    }

    /* calculating the displacements and the stresses and storing  
       the results in frd format */
  
    MNEW(v,double,mt**nk);
    MNEW(fn,double,mt**nk);
    NNEW(stn,double,6**nk);
    if(*ithermal>1) NNEW(qfn,double,3**nk);
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
    if(strcmp1(&filab[435],"PEEQ")==0) NNEW(epn,double,*nk);
    if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
    if(strcmp1(&filab[609],"SDV ")==0) NNEW(xstaten,double,*nstate_**nk);
    if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
    if(strcmp1(&filab[4785],"DUCT")==0) NNEW(damn,double,*nk);
    
    isiz=mt**nk;cpypardou(v,vold,&isiz,&num_cpus);
    iout=2;
    icmd=3;

#ifdef COMPANY
    FORTRAN(uinit,());
#endif
    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	    prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xbounact,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
            &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,&icmd,
            ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
            xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
            ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
	    nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
            &reltime,&ne0,thicke,shcon,nshcon,
            sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
            mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoeln,
            inoeln,nener,orname,network,ipobody,xbodyact,ibody,typeboun,
	    itiefac,tieset,smscale,&mscalmethod,nbody,t0g,t1g,
	    islavquadel,aut,irowt,jqt,&mortartrafoflag,
	    &intscheme,physcon,dam,damn,iponoel);
    
    isiz=mt**nk;cpypardou(vold,v,&isiz,&num_cpus);

    iout=0;
    if(*iexpl<=1) icmd=0;
    
    ++*kode;
    if(*mcs>0){
      ptime=*ttime+time;
      if(*mortar>1){
	mortar_prefrd(ne,nslavs,mi,nk,nkon, &stx,cdisp,fn,cfs,cfm);       
      }
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,
	     t1act,fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
             nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,qfn,
             ialset,istartset,iendset,trab,inotr,ntrans,orab,ielorien,
	     norien,stx,veold,&noddiam,set,nset,emn,thicke,jobnamec,&ne0,
             cdn,mortar,nmat,qfx,ielprop,prop,damn,&errn);
      if(*mortar>1){
	mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
      }
#ifdef COMPANY
      FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif

    }else{
      if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	MNEW(ipneigh,ITG,*nk);
      }

      ptime=*ttime+time;
      if(*mortar>1){
	mortar_prefrd(ne,nslavs,mi,nk,nkon, &stx,cdisp,fn,cfs,cfm);       
      }
      frd(co,nk,kon,ipkon,lakon,&ne0,v,stn,inum,nmethod,
	  kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	  nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	  ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	  mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	  cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	  thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat,ielprop,
	  prop,sti,damn,&errn);
      if(*mortar>1){
	mortar_postfrd(ne,nslavs,mi,nk,nkon,fn,cfs,cfm);      
      }

      if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
#ifdef COMPANY
      FORTRAN(uout,(v,mi,ithermal,filab,kode,output,jobnamec));
#endif
    }

    /* mesh refinement */
  
    if(strcmp1(&filab[4089],"RM")==0){
      refinemesh(nk,ne,co,ipkon,kon,v,veold,stn,een,emn,epn,enern,
		 qfn,errn,filab,mi,lakon,jobnamec,istartset,iendset,
		 ialset,set,nset,matname,ithermal,output,nmat,
		 nelemload,nload,sideload,nodeforc,
		 nforc,nodeboun,nboun,nodempc,ipompc,nmpc);

      /* free errn */
	
      if(((*nmethod!=5)||(mode==-1))&&
	 ((strcmp1(&filab[1044],"ERR")==0)&&(*ithermal!=2))) SFREE(errn);
    }

    SFREE(v);SFREE(fn);SFREE(stn);SFREE(inum);SFREE(stx);
    if(*ithermal>1){SFREE(qfn);}
    
    if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
    if(strcmp1(&filab[435],"PEEQ")==0) SFREE(epn);
    if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
    if(strcmp1(&filab[609],"SDV ")==0) SFREE(xstaten);
    if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);
    if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
    if(strcmp1(&filab[4785],"DUCT")==0) SFREE(damn);

  }
    
  /* writing out the latest stiffness matrix for a subsequent
     sensitivity analysis */

  if(isensitivity){
      
    strcpy2(stiffmatrix,jobnamec,132);
    strcat(stiffmatrix,".stm");
      
    if((f1=fopen(stiffmatrix,"wb"))==NULL){
      printf(" *ERROR in nonlingeo: cannot open stiffness matrix file for writing...");
      exit(0);
    }
      
    /* storing the stiffness matrix */

    /* nzs,irow,jq and icol have to be stored too, since the static analysis
       can involve contact, whereas in the sensitivity analysis contact is not
       taken into account while determining the structure of the stiffness
       matrix (in mastruct.c)
    */
      
    if(fwrite(&nasym,sizeof(ITG),1,f1)!=1){
      printf(" *ERROR in nonlingeo saving the symmetry flag to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(nzs,sizeof(ITG),3,f1)!=3){
      printf(" *ERROR in nonlingeo saving the number of subdiagonal nonzeros to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(irow,sizeof(ITG),nzs[2],f1)!=nzs[2]){
      printf(" *ERROR in nonlingeo saving irow to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(jq,sizeof(ITG),neq[1]+1,f1)!=neq[1]+1){
      printf(" *ERROR in nonlingeo saving jq to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(icol,sizeof(ITG),neq[1],f1)!=neq[1]){
      printf(" *ERROR in nonlingeo saving icol to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(adcpy,sizeof(double),neq[1],f1)!=neq[1]){
      printf(" *ERROR in nonlingeo saving the diagonal of the stiffness matrix to the stiffness matrix file...");
      exit(0);
    }
    if(fwrite(aucpy,sizeof(double),(nasym+1)*nzs[2],f1)!=(nasym+1)*nzs[2]){
      printf(" *ERROR in nonlingeo saving the off-diagonal terms of the stiffness matrix to the stiffness matrix file...");
      exit(0);
    }
    fclose(f1);
    SFREE(adcpy);SFREE(aucpy);
  }
  
  /* restoring the distributed loading  */

  if((*ithermal==3)&&(ncont!=0)&&(*mortar==1)&&(*ncmat_>=11)){
    *nload=nloadref;
    RENEW(nelemload,ITG,2**nload);
    isiz=2**nload;cpyparitg(nelemload,nelemloadref,&isiz,&num_cpus);
    if(*nam>0){
      RENEW(iamload,ITG,2**nload);
      isiz=2**nload;cpyparitg(iamload,iamloadref,&isiz,&num_cpus);
    }
    RENEW(sideload,char,20**nload);memcpy(&sideload[0],&sideloadref[0],
					  sizeof(char)*20**nload);
      
    /* freeing the temporary fields */
      
    SFREE(nelemloadref);if(*nam>0){SFREE(iamloadref);};
    SFREE(sideloadref);
  }

  /* setting the velocity to zero at the end of a quasistatic or stationary
     step */

  if(abs(*nmethod)==1){
    for(k=0;k<mt**nk;++k){
      veold[k]=0.;}
  }

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  for(k=0;k<*nboun;++k){

    /* thermal boundary conditions are updated only if the
       step was thermal or thermomechanical */

    if(ndirboun[k]==0){
      if(*ithermal<2) continue;

      /* mechanical boundary conditions are updated only
	 if the step was not thermal or the node is a
	 network node */

    }else if((ndirboun[k]>0)&&(ndirboun[k]<4)){
      node=nodeboun[k];
      FORTRAN(nident,(itg,&node,&ntg,&id));
      networknode=0;
      if(id>0){
	if(itg[id-1]==node) networknode=1;
      }
      if((*ithermal==2)&&(networknode==0)) continue;
    }
    xbounold[k]=xbounact[k];
  }
  isiz=*nforc;cpypardou(xforcold,xforcact,&isiz,&num_cpus);
  isiz=2**nload;cpypardou(xloadold,xloadact,&isiz,&num_cpus);
  isiz=7**nbody;cpypardou(xbodyold,xbodyact,&isiz,&num_cpus);
  if(*ithermal==1){
    cpypardou(t1old,t1act,nk,&num_cpus);
    for(k=0;k<*nk;++k){
      vold[mt*k]=t1act[k];}
  }
  else if(*ithermal>1){
    for(k=0;k<*nk;++k){
      t1[k]=vold[mt*k];}
    if(*ithermal>=3){
      cpypardou(t1old,t1act,nk,&num_cpus);
    }
  }

  qaold[0]=qa[0];
  qaold[1]=qa[1];
  
  if(*iexpl>1){
    SFREE(smscale);

    if((mscalmethod==1)||(mscalmethod==3)||(*mortar==-1)){
      if(*isolver==0){
#ifdef SPOOLES
	spooles_cleanup();
#endif
      }
      else if(*isolver==4){
#ifdef SGI
	sgi_cleanup(token);
#endif
      }
      else if(*isolver==5){
#ifdef TAUCS
	tau_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	pardiso_cleanup(&neq[0],&symmetryflag,&inputformat);
#endif
      }
      else if(*isolver==8){
#ifdef PASTIX
#endif
      }
    }
  }
  
  SFREE(f);SFREE(b);
  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(xbodyact);

  if(*inewton==1){SFREE(cgr);}
  SFREE(fext);SFREE(ampli);SFREE(xbounini);SFREE(xstiff);
  if((*ithermal==1)||(*ithermal>=3)){SFREE(t1act);SFREE(t1ini);}

  if(*ithermal>1){
    SFREE(itg);SFREE(ieg);SFREE(kontri);SFREE(nloadtr);
    SFREE(nactdog);SFREE(nacteq);SFREE(ineighe);
    SFREE(tarea);SFREE(tenv);SFREE(fenv);SFREE(qfx);
    SFREE(erad);SFREE(ac);SFREE(bc);SFREE(ipiv);
    SFREE(bcr);SFREE(ipivr);SFREE(adview);SFREE(auview);SFREE(adrad);
    SFREE(aurad);SFREE(irowrad);SFREE(jqrad);SFREE(icolrad);
    if((*mcs>0)&&(ntr>0)){SFREE(inocs);}
    if((*network>0)||(ntg>0)){SFREE(iponoeln);SFREE(inoeln);}
    if(ntr>0){
    }
  }

  if(icfd==1){
  }else if(icfd==2){
    SFREE(sideface);SFREE(nelemface);SFREE(ifreestream);
    SFREE(isolidsurf);SFREE(neighsolidsurf);SFREE(iponoelf);SFREE(inoelf);
    SFREE(inomat);SFREE(ipface);
    if(*ithermal==1) SFREE(qfx);
	 
    SFREE(ipkonf);SFREE(lakonf);SFREE(ielmatf);SFREE(nelold);SFREE(nelnew);
    SFREE(cof);SFREE(voldf);SFREE(nkold);SFREE(nknew);SFREE(konf);
    SFREE(ipompcf);SFREE(nodempcf);SFREE(coefmpcf);SFREE(nodebounf);
    SFREE(ndirbounf);SFREE(xbounf);SFREE(nelemloadf);SFREE(xloadf);
    SFREE(sideloadf);SFREE(ikbounf);SFREE(ilbounf);
    SFREE(ikmpcf);SFREE(ilmpcf);SFREE(xbounoldf);SFREE(xbounactf);
    SFREE(xloadoldf);SFREE(xloadactf);SFREE(inotrf);
    if(*norien>0) SFREE(ielorienf);
    if(*nbody>0) SFREE(ipobodyf);
    if(*nam>0){SFREE(iambounf);SFREE(iamloadf);}
  }

  SFREE(fini);
  if(*nmethod==4){
    SFREE(aux2);SFREE(fextini);SFREE(veini);SFREE(accini);
    SFREE(adb);SFREE(aub);SFREE(cvini);SFREE(cv);SFREE(fnext);
    SFREE(fnextini);
  }
  SFREE(eei);SFREE(stiini);SFREE(emeini);
  if(*nener==1)SFREE(enerini);
  if(*nstate_!=0){SFREE(xstateini);}

  SFREE(aux);SFREE(iaux);SFREE(vini);

  if((*ndmat_>0)&&(*iexpl<=1)){
    /* The results context owns no storage.  Clear it before releasing the
       nonlingeo damage baselines so later steps/output calls cannot retain
       dangling pointers. */
    results_set_de12_context(0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0.);
    if(damage_tent_elem!=NULL){
      SFREE(damage_tent_elem); damage_tent_elem=NULL;
      SFREE(damage_tent_mat); damage_tent_mat=NULL;
      SFREE(damage_tent_ip); damage_tent_ip=NULL;
      SFREE(damage_tent_value); damage_tent_value=NULL;
    }
    damage_tent_count=0;

    if(fdamage!=NULL){
      fflush(fdamage);
      fclose(fdamage);
      fdamage=NULL;
    }

    SFREE(ipkondamageini);
    SFREE(damdamageini);
    SFREE(damde1prev);
    if(damage_damjac!=NULL) SFREE(damage_damjac);
    if(damage_fracture_seta!=NULL) free(damage_fracture_seta);
    if(damage_diss_fhat!=NULL) SFREE(damage_diss_fhat);
    if(damage_diss_uf!=NULL) SFREE(damage_diss_uf);
    if(damage_damvisc!=NULL) SFREE(damage_damvisc);
    if(damage_damviscini!=NULL) SFREE(damage_damviscini);
    if(damage_bt_dam!=NULL) SFREE(damage_bt_dam);
    if(damage_bt_visc!=NULL) SFREE(damage_bt_visc);
    if(damage_bt_xs!=NULL) SFREE(damage_bt_xs);
    if(damage_evt_sgn!=NULL) SFREE(damage_evt_sgn);
    if(damage_ray_r0!=NULL) SFREE(damage_ray_r0);
    if(damage_ray_cat!=NULL) SFREE(damage_ray_cat);
    if(damage_ray_p!=NULL) SFREE(damage_ray_p);
    if(damage_ray_res!=NULL) SFREE(damage_ray_res);
    if(damage_ct_sgn!=NULL) SFREE(damage_ct_sgn);
    if(damage_ct_ring!=NULL) SFREE(damage_ct_ring);
    if(damage_ct_fl!=NULL) SFREE(damage_ct_fl);
    if(damage_ct_bl!=NULL) SFREE(damage_ct_bl);
    if(damage_ct_w!=NULL) SFREE(damage_ct_w);
    if(damage_ct_r0!=NULL) SFREE(damage_ct_r0);
    if(damage_ct_beps!=NULL) SFREE(damage_ct_beps);
    if(damage_ct_y!=NULL) SFREE(damage_ct_y);
    if(damage_ct_z!=NULL) SFREE(damage_ct_z);
    if(damage_ct_qh!=NULL) SFREE(damage_ct_qh);
    if(damage_ct_dam!=NULL) SFREE(damage_ct_dam);
    if(damage_ct_visc!=NULL) SFREE(damage_ct_visc);
    if(damage_ct_xs!=NULL) SFREE(damage_ct_xs);
    if(damage_ct_jac!=NULL) SFREE(damage_ct_jac);
    if(damage_dl_r0!=NULL) SFREE(damage_dl_r0);
    if(damage_dl_d!=NULL) SFREE(damage_dl_d);
    if(damage_dl_w!=NULL) SFREE(damage_dl_w);
    if(damage_dl_pn!=NULL) SFREE(damage_dl_pn);
    if(damage_dl_res!=NULL) SFREE(damage_dl_res);
    if(damage_dl_dam!=NULL) SFREE(damage_dl_dam);
    if(damage_dl_visc!=NULL) SFREE(damage_dl_visc);
    if(damage_dl_xs!=NULL) SFREE(damage_dl_xs);
    if(damage_dl_mode==1){
      printf("[DAMAGE TR] SUMMARY: armed %" ITGFORMAT " time(s); accepted "
             "steps %" ITGFORMAT " (Newton %" ITGFORMAT ", Cauchy %" ITGFORMAT
             ", dogleg %" ITGFORMAT "); iterations that accepted nothing %"
             ITGFORMAT "; rejected trials %" ITGFORMAT "; residual "
             "evaluations %" ITGFORMAT "; armed factorisations %" ITGFORMAT
             "; last transpose check %.12e; operator asymmetry %.6e\n",
             damage_dl_narm,damage_dl_nacc,damage_dl_nnewt,damage_dl_ncau,
             damage_dl_ndog,damage_dl_nfail,damage_dl_nrej,damage_dl_neval,
             damage_dl_nfact,damage_dl_ident,damage_dl_asym);
      fflush(stdout);
    }
    if(damage_frel!=NULL) SFREE(damage_frel);
    if(damage_ract!=NULL) SFREE(damage_ract);
    if(damage_de13_trigger_value!=NULL) SFREE(damage_de13_trigger_value);
    if(damage_de13_trigger_ip!=NULL) SFREE(damage_de13_trigger_ip);
    if(*nmethod!=4) SFREE(veolddamageini);
  }

  if(icascade==2){
    memmpc_=memmpcref_;mpcfree=mpcfreeref;maxlenmpc=maxlenmpcref;
    RENEW(nodempc,ITG,3*memmpcref_);
    for(k=0;k<3*memmpcref_;k++){
      nodempc[k]=nodempcref[k];}
    RENEW(coefmpc,double,memmpcref_);
    for(k=0;k<memmpcref_;k++){
      coefmpc[k]=coefmpcref[k];}
    SFREE(nodempcref);SFREE(coefmpcref);
  }

  if(ncont!=0){
    *ne=ne0;*nkon=nkon0;
    if((*nener==1)&&(*mortar==1)){
      RENEW(ener,double,mi[0]**ne*2);
    }
    RENEW(ipkon,ITG,*ne);
    RENEW(lakon,char,8**ne);
    RENEW(kon,ITG,*nkon);
    if(*norien>0){
      RENEW(ielorien,ITG,mi[2]**ne);
    }
    RENEW(ielmat,ITG,mi[2]**ne);

    if(*mortar>1){
      
      /// needed for next step coloumb friction
      
      for (i=0;i<*ntie;i++){
	if(tieset[i*(81*3)+80]=='C'){
	  if(*nstate_*mi[0]>0){
	    for(j=nslavnode[i];j<nslavnode[i+1];j++){	  	     
	      for(k=0;k<3;k++){	    	       
		xstate[*nstate_*mi[0]*(*ne+j)+k]=cstress[mt*j+k]; 
	      } 	    	       
	      xstate[*nstate_*mi[0]*(*ne+j)+3]=islavact[j]+0.5;
	    }
		      
	  }
	}
      }
    }
      
    SFREE(cg);SFREE(straight);
    SFREE(imastop);SFREE(itiefac);SFREE(islavnode);
    SFREE(nslavnode);SFREE(iponoels);SFREE(inoels);SFREE(imastnode);
    SFREE(nmastnode);SFREE(itietri);SFREE(koncont);SFREE(xnoels);
    SFREE(springarea);SFREE(xmastnor);

    if(*mortar==-1){
      if(ncont!=0){SFREE(kslav);SFREE(lslav);SFREE(ktot);SFREE(ltot);
	SFREE(aloc);SFREE(alglob);SFREE(areaslav);SFREE(fric);}
      SFREE(adc);SFREE(auc);
      if(idispfrdonly==1){SFREE(inumcp);}
      if(masslesslinear>0){
	SFREE(ad);SFREE(au);
	if(ncont!=0){SFREE(auw);SFREE(jqw);SFREE(iroww);
	  SFREE(fullgmatrix);SFREE(fullr);}
	iclean=1;
        massless(kslav,lslav,ktot,ltot,au,ad,auc,adc,jq,irow,neq,nzs,auw,jqw,
		 iroww,&nzsw,islavnode,nslavnode,nslavs,imastnode,nmastnode,
		 ntie,nactdof,mi,vold,volddof,veold,nk,fext,isolver,
		 &masslesslinear,co,springarea,&neqtot,qb,b,&dtime,aloc,fric,
		 iexpl,nener,ener,ne,&jqbi,&aubi,&irowbi,&jqib,&auib,&irowib,
		 &iclean,&iinc,fullgmatrix,fullr,alglob,&num_cpus,&ncont);
      }
      if(masslesslinear==2){SFREE(fextload);}

    }else if(*mortar==0){
      SFREE(areaslav);
    }else if(*mortar==1){
      SFREE(pmastsurf);SFREE(ipe);SFREE(ime);
      SFREE(islavact);
    }else if(*mortar>1){
      SFREE(islavact);SFREE(gap);SFREE(slavnor);SFREE(slavtan);
      SFREE(cstress);SFREE(ipe);SFREE(ime);SFREE(cfs);SFREE(cfm);
      SFREE(cdisp);SFREE(bp);SFREE(islavtie);
      SFREE(nslavspc);SFREE(islavspc);SFREE(nslavmpc);SFREE(islavmpc);
      SFREE(nmastspc);SFREE(imastspc);SFREE(nmastmpc);SFREE(imastmpc);
      SFREE(pslavdual);
      SFREE(cstressini);SFREE(bpini);SFREE(islavactini);
      SFREE(aut);SFREE(irowt);SFREE(jqt);
      SFREE(autinv);SFREE(irowtinv);SFREE(jqtinv);
      SFREE(Bd);SFREE(irowb);SFREE(jqb);
      SFREE(Bdhelp);SFREE(irowbhelp);SFREE(jqbhelp);
      SFREE(Dd);SFREE(irowd);SFREE(jqd);
      SFREE(Ddtil);SFREE(irowdtil);SFREE(jqdtil);
      SFREE(Bdtil);SFREE(irowbtil);SFREE(jqbtil);
      SFREE(islavnodeinv);SFREE(islavquadel);
    }
  }

  /* reset icascade */

  if(icascade==1){icascade=0;}

  mpcinfo[0]=memmpc_;mpcinfo[1]=mpcfree;mpcinfo[2]=icascade;
  mpcinfo[3]=maxlenmpc;

  if(iglob==1){SFREE(integerglob);SFREE(doubleglob);}

  *icolp=icol;*irowp=irow;*cop=co;*voldp=vold;

  *ipompcp=ipompc;*labmpcp=labmpc;*ikmpcp=ikmpc;*ilmpcp=ilmpc;
  *fmpcp=fmpc;*nodempcp=nodempc;*coefmpcp=coefmpc;*nelemloadp=nelemload;
  *iamloadp=iamload;*sideloadp=sideload;

  *ipkonp=ipkon;*lakonp=lakon;*konp=kon;*ielorienp=ielorien;
  *ielmatp=ielmat;*enerp=ener;*xstatep=xstate;

  *islavsurfp=islavsurf;*pslavsurfp=pslavsurf;*clearinip=clearini;

  (*tmin)*=(*tper);
  (*tmax)*=(*tper);

  SFREE(nactdofinv);
  // MPADD start
  if((*nmethod==4)&&(*ithermal!=2)&&(*iexpl<=1)&&(icfd==0)){ SFREE(adblump);}
  // MPADD end
  
  (*ttime)+=(*tper);

  SFREE(iponoel);
  

  return;
}
