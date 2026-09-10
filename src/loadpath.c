/*     CalculiX - damage/fracture extension                              */
/*     loadpath.c: THE judgement about whether the specimen is still     */
/*                 one piece.                                            */

/* Why this module exists
   ----------------------
   `damstate.c` owns the judgement about a single NODE losing its load path.
   This module owns the judgement one level up: whether a load path still
   links the grips at all - that is, whether what is being solved is still a
   specimen.

   Before this module that judgement existed but was reachable only through a
   chain of three independent conditions, and every one of them was off by
   default.  Measured on the fast plain deck, 1 thread, PARDISO:

     - the specimen loses grip-to-grip connectivity at increment 65,
       theta = 0.1175;
     - the run continues to increment 507, theta = 1.0, deleting nothing
       after increment 65;
     - the grip reaction over that stretch is a straight line through the
       origin - Fx/theta constant to 0.04% over theta 0.19..1.0 - because
       the failed facets are a linear spring of stiffness gmin*Kn joining
       two halves that are no longer attached.

   So the specimen breaks, and the model then reports a load that RISES
   linearly, from 1.2% of peak back to 6.5% of peak.  That is the defect: not
   a wrong number, a plausible one.

   The three gaps, each measured on this branch, were:

     1. the judgement was asked only when CCX_FRACTURE_TERMINATION named two
        node sets by hand, so it was off for any deck nobody had configured;
     2. when asked, a fully failed cohesive facet still counted as
        conducting unless a SECOND switch, CCX_FRACTURE_DEADFACET, was also
        set - and test/s3rad/run_s3rad.sh explicitly unset it;
     3. the call site sat inside `if(damage_tent_count>0)`, the bulk-deletion
        transaction, so on an interface-dominated fracture - where facets
        fail and no bulk erodes - it could never fire at all.

   Gap 3 is the one that makes the other two hard to see: the guard could be
   configured correctly and still never run.

   Contract
   --------
     - one census per converged increment, asked of this module, never
       recomputed by a consumer;
     - an element conducts iff it is not deleted AND it is not a fully failed
       cohesive facet.  The facet half of that is asked of
       damstate_facet_dead - this module composes the judgement, it does not
       re-derive it;
     - the graph walk is the existing, in-use Fortran kernel (damconnect /
       damconnectface).  There is one implementation of the walk;
     - the endpoints come from the deck's own *BOUNDARY cards when nobody
       names them, so the judgement does not depend on configuration;
     - severance LATCHES.  A specimen that has come apart does not become a
       specimen again, and every increment after it is stamped;
     - the module refuses to arm if its self test fails, the discipline
       lsladder.c set and damstate.c followed.

   It changes no equation and takes no decision of its own.  Whether a run
   stops at severance is a policy its consumer applies; this module only
   answers, and says so loudly enough that a post-severance number cannot be
   quoted by accident.                                                    */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"

/* ------------------------------------------------------------ endpoints */

/* Derive the grips from the boundary conditions.

   The load path that matters is the one between the degrees of freedom that
   are DRIVEN and the ones that REACT in the same direction.  A dof held at
   zero in some other direction is a rigid-body suppressor, not a grip: on
   the target deck FIXPOINTA and FIXPOINTB hold directions 2 and 3 and must
   not be mistaken for the reacting face.

   So: pick the translational direction that carries the largest prescribed
   magnitude; the driven set is every node prescribed nonzero in it, the held
   set every node prescribed zero in it.  On every deck in this tree that
   reproduces the hand-named pair exactly - FACE_XL_NSET against
   FACE_X0_NSET - which is the point: the judgement stops depending on
   somebody having configured it.

   Returns the driven direction, or 0 when the deck does not present a
   grip pair and the module must not arm. */
ITG loadpath_derive(const ITG *nodeboun,const ITG *ndirboun,const double *xboun,
                    ITG nboun,ITG nk,ITG **held,ITG *nheld,
                    ITG **driven,ITG *ndriven)
{
  ITG k,d,idir=0,na=0,nb=0,*a=NULL,*b=NULL;
  double mag,best=0.;
  char *seena=NULL,*seenb=NULL;

  *held=NULL; *driven=NULL; *nheld=0; *ndriven=0;
  if((nboun<=0)||(nk<=0)) return 0;

  /* which direction is driven */
  for(k=0;k<nboun;k++){
    d=ndirboun[k];
    if((d<1)||(d>3)) continue;
    mag=xboun[k]; if(mag<0.) mag=-mag;
    if(mag>best){best=mag; idir=d;}
  }
  if((idir==0)||(best<=0.)) return 0;

  /* a node may appear more than once; count each at most once per side */
  NNEW(seena,char,nk+1);
  NNEW(seenb,char,nk+1);
  for(k=0;k<nboun;k++){
    ITG node=nodeboun[k];
    if(ndirboun[k]!=idir) continue;
    if((node<1)||(node>nk)) continue;
    mag=xboun[k]; if(mag<0.) mag=-mag;
    if(mag>0.){ if(!seenb[node]){seenb[node]=1;nb++;} }
    else      { if(!seena[node]){seena[node]=1;na++;} }
  }
  /* A node both driven and held in the same direction is over-specified;
     the driven side wins, and it is removed from the held side so the two
     lists cannot share a node and read as trivially connected. */
  for(k=1;k<=nk;k++) if(seenb[k]&&seena[k]){seena[k]=0;na--;}

  if((na<=0)||(nb<=0)){SFREE(seena);SFREE(seenb);return 0;}

  NNEW(a,ITG,na); NNEW(b,ITG,nb);
  na=0; nb=0;
  for(k=1;k<=nk;k++){
    if(seena[k]) a[na++]=k;
    if(seenb[k]) b[nb++]=k;
  }
  SFREE(seena); SFREE(seenb);
  *held=a; *nheld=na; *driven=b; *ndriven=nb;
  return idir;
}

/* ----------------------------------------------------------------- arm */

void loadpath_init(loadpath *lp)
{
  memset(lp,0,sizeof(*lp));
  lp->connected=1;
  lp->sev_inc=-1;
  /* Three consecutive converged increments.  One is what a flickering
     judgement can produce; three is still nothing against the 178
     increments the unguarded defect used to spend. */
  lp->nconfirm=3;
}

void loadpath_free(loadpath *lp)
{
  if(lp->nodesa!=NULL) SFREE(lp->nodesa);
  if(lp->nodesb!=NULL) SFREE(lp->nodesb);
  lp->nodesa=NULL; lp->nodesb=NULL;
  lp->armed=0; lp->na=0; lp->nb=0;
}

/* Resolve the endpoints and arm.  A named pair (CCX_FRACTURE_TERMINATION)
   wins over the derivation, so an existing configuration keeps its exact
   meaning; the derivation is what happens when nobody configured anything.
   Returns 1 when armed. */
ITG loadpath_arm(loadpath *lp,ITG nk,
                 char *set,ITG *nset,ITG *istartset,ITG *iendset,ITG *ialset,
                 const char *namea,const char *nameb,
                 const ITG *nodeboun,const ITG *ndirboun,const double *xboun,
                 ITG nboun,ITG imode)
{
  ITG na=0,nb=0,*a=NULL,*b=NULL,idir,dummy[1];
  char pa[81],pb[81];

  /* re-arming on a later step must not leak the previous endpoint lists */
  if(lp->nodesa!=NULL){SFREE(lp->nodesa); lp->nodesa=NULL;}
  if(lp->nodesb!=NULL){SFREE(lp->nodesb); lp->nodesb=NULL;}
  lp->na=0; lp->nb=0;
  lp->nk=nk; lp->imode=imode;

  if((namea!=NULL)&&(nameb!=NULL)){
    /* blank padded to 81: the FORTRAN macro passes no hidden string
       length, so an assumed-length dummy there would read garbage */
    memset(pa,' ',81); memset(pb,' ',81);
    memcpy(pa,namea,(strlen(namea)<80)?strlen(namea):80);
    memcpy(pb,nameb,(strlen(nameb)<80)?strlen(nameb):80);
    dummy[0]=0; na=0; nb=0;
    FORTRAN(damsetnodes,(set,nset,istartset,iendset,ialset,pa,dummy,&na));
    FORTRAN(damsetnodes,(set,nset,istartset,iendset,ialset,pb,dummy,&nb));
    if((na>0)&&(nb>0)){
      NNEW(a,ITG,na); NNEW(b,ITG,nb);
      FORTRAN(damsetnodes,(set,nset,istartset,iendset,ialset,pa,a,&na));
      FORTRAN(damsetnodes,(set,nset,istartset,iendset,ialset,pb,b,&nb));
      lp->nodesa=a; lp->na=na; lp->nodesb=b; lp->nb=nb;
      lp->idir=0; lp->armed=1;
      snprintf(lp->origin,sizeof(lp->origin),"named sets %s:%s",namea,nameb);
      return 1;
    }
    printf("[LOADPATH] *WARNING: %s or %s is not a node set in this deck; "
           "falling back to the deck's own *BOUNDARY cards\n",namea,nameb);
  }

  idir=loadpath_derive(nodeboun,ndirboun,xboun,nboun,nk,&a,&na,&b,&nb);
  if(idir==0){
    printf("[LOADPATH] no grip pair in this deck: no direction carries both "
           "a driven and a reacting boundary condition.  The specimen "
           "cannot be told when it has come apart, so nothing is judged.\n");
    return 0;
  }
  lp->nodesa=a; lp->na=na; lp->nodesb=b; lp->nb=nb;
  lp->idir=idir; lp->armed=1;
  snprintf(lp->origin,sizeof(lp->origin),
           "derived from *BOUNDARY, direction %d: %d reacting node(s) "
           "against %d driven",(int)idir,(int)na,(int)nb);
  return 1;
}

/* -------------------------------------------------------------- census */

/* Walk the surviving topology and answer.  Cheap: one pass to mark the dead
   facets, then the existing kernel, which is O(nkon+nk).  Against a PARDISO
   solve this is free, which is the whole reason it can run every increment
   instead of being a python script somebody remembers to run afterwards. */
/* Is a facet OPEN?

   A fully failed facet still carries load while its faces are pressed
   together - that is how a real crack transmits compression, and
   cohesive_uc6.f says so explicitly: for deltal(1) < 0 the normal stiffness
   is the FULL kn, not g*kn.  So "failed" is not the same as "not a load
   path", and treating it as such was a modelling decision dressed up as a
   topological fact.

   The normal traction is the sign that decides it, and the law already
   writes it to stx slot 1 for every integration point:

     deltal(1) >= 0   traction(1) = g*kn*deltal(1)                   >= 0
     deltal(1) <  0   traction(1) = kn*deltal(1)                     <  0
     inside the CCX_UC6_CONTACT_SMOOTH band, deltal(1) < 0:
                      traction(1) = g*kn*deltal(1) - (1-g)*kn*psism  <  0

   so traction(1) < 0 means compression in every branch, smoothed or sharp,
   and needs no threshold.  ONE point in compression is enough: the facet
   transmits through that point.

   With stx unavailable the facet is reported open, which is the judgement
   this module made before it could read the traction at all. */
ITG loadpath_facet_open(const double *stx,ITG mi0,ITG elem,ITG nip)
{
  ITG j;
  if(stx==NULL) return 1;
  if(nip>mi0) nip=mi0;
  for(j=0;j<nip;j++){
    if(stx[6*(mi0*elem+j)]<0.) return 0;
  }
  return 1;
}

ITG loadpath_census(loadpath *lp,ITG *ipkon,ITG *kon,char *lakon,ITG *ne,
                    const double *xstate,ITG nstate,ITG mi0,
                    const double *stx)
{
  ITG i,*ifacdead=NULL,iconn=1,nreach=0;

  if(!lp->armed) return 1;

  NNEW(ifacdead,ITG,*ne);
  lp->nfacet=0; lp->nfacetdead=0; lp->nfacetshut=0; lp->nlive=0;
  for(i=0;i<*ne;i++){
    if(ipkon[i]<0) continue;
    lp->nlive++;
    if(lakon[8*i]!='U') continue;
    lp->nfacet++;
    /* the facet half of the judgement has an owner; ask it.  Failure alone
       does not settle it - a failed facet in compression still conducts. */
    if((xstate!=NULL)&&(nstate>=4)&&
       damstate_facet_dead(xstate,nstate,mi0,i,3)){
      if(loadpath_facet_open(stx,mi0,i,3)){
        ifacdead[i]=1;
        lp->nfacetdead++;
      }else{
        lp->nfacetshut++;
      }
    }
  }

  if(lp->imode==1){
    FORTRAN(damconnectface,(ipkon,kon,lakon,ne,&lp->nk,lp->nodesa,&lp->na,
                            lp->nodesb,&lp->nb,&iconn,&nreach,ifacdead));
  }else{
    FORTRAN(damconnect,(ipkon,kon,lakon,ne,&lp->nk,lp->nodesa,&lp->na,
                        lp->nodesb,&lp->nb,&iconn,&nreach,ifacdead));
  }
  SFREE(ifacdead);

  lp->connected=iconn;
  lp->nreach=nreach;
  return iconn;
}

/* The latch, separated from the printing so the self test can drive it.

   Two mechanisms, and they answer different worries.

   HYSTERESIS guards against a judgement that flickers.  Severance must be
   confirmed on `nconfirm` CONSECUTIVE converged increments before it
   latches; one increment that reconnects resets the count.  This is the
   honest way to buy robustness - the alternative, a policy switch that
   suppresses the judgement, hides it instead.  Note what is recorded: the
   FIRST increment without a path, not the increment the latch confirmed it,
   because the first is the physical event and the rest is confirmation.

   The LATCH keeps severance from un-happening once confirmed.  A specimen
   that has come apart does not become a specimen again, even though its
   faces may touch and legitimately carry compression afterwards.

   The cost of a false positive is stopping a live run; the cost of
   confirming is `nconfirm` increments against the 178 the defect used to
   waste.  That asymmetry is why the default is small and not one.

   Returns 1 on the increment severance is CONFIRMED, 2 on every increment
   after it, 3 while a loss of path is provisional, 0 while intact. */
ITG loadpath_latch(loadpath *lp,ITG inc,double t)
{
  if(!lp->armed) return 0;
  if(lp->sev_seen){lp->ninc_past++; return 2;}
  if(lp->connected==0){
    if(lp->ndisc==0){lp->sev_inc=inc; lp->sev_time=t;}
    lp->ndisc++;
    if(lp->ndisc>=lp->nconfirm){lp->sev_seen=1; return 1;}
    return 3;
  }
  lp->ndisc=0;
  return 0;
}

/* Latch and report.  Severance is announced ONCE, in full; after that every
   increment carries a one-line stamp, so a number quoted from the log after
   this point states its own status on the same screen. */
void loadpath_note(loadpath *lp,ITG inc,double t)
{
  ITG what=loadpath_latch(lp,inc,t);

  if(what==1){
    printf("\n[LOADPATH SEVERED] inc=%" ITGFORMAT " time=%.12e\n"
           "                   no surviving load path between the grips\n"
           "                   (confirmed on %" ITGFORMAT " consecutive "
           "converged increments; the increment above is the FIRST "
           "without a path)\n"
           "                   endpoints: %s\n"
           "                   nodes still reachable from the reacting "
           "set: %" ITGFORMAT "\n"
           "                   live elements %" ITGFORMAT ", cohesive "
           "facets %" ITGFORMAT " of which %" ITGFORMAT " fully failed\n"
           "                   everything after this increment is two "
           "pieces, not a specimen\n\n",
           lp->sev_inc,lp->sev_time,lp->nconfirm,lp->origin,lp->nreach,
           lp->nlive,lp->nfacet,lp->nfacetdead);
    fflush(stdout);
    return;
  }

  /* The standing census.  02-DIAGNOSTICS.md 6 asked for this: connectivity
     as a quantity the solver prints, cheaply, rather than an offline python
     script somebody remembers to run afterwards.  Every 50 converged
     increments is a handful of lines per run and makes 7 - "bisect the
     deletion history for severance" - unnecessary, because the run already
     said it. */
  if((what==0)&&(inc>0)&&(inc%50==0)){
    printf("[LOADPATH CENSUS] inc=%" ITGFORMAT " time=%.7f connected=%"
           ITGFORMAT " reach=%" ITGFORMAT " live=%" ITGFORMAT
           " facets=%" ITGFORMAT "/%" ITGFORMAT " open-failed, %"
           ITGFORMAT " failed-but-shut\n",
           inc,t,lp->connected,lp->nreach,lp->nlive,
           lp->nfacetdead,lp->nfacet,lp->nfacetshut);
    fflush(stdout);
  }

  if(what==3){
    printf("[LOADPATH PROVISIONAL] inc=%" ITGFORMAT " time=%.7f no load path,"
           " %" ITGFORMAT " of %" ITGFORMAT " consecutive converged"
           " increment(s); not severance until confirmed\n",
           inc,t,lp->ndisc,lp->nconfirm);
    fflush(stdout);
  }

  if(what==2){
    /* Say WHAT joins the pieces, not just that they are two.  05-DEBT.md
       item 1 is stated in exactly those terms - "two separated halves
       joined by 451 facets, 322 at the residual floor" - and a stamp that
       omits it cannot be checked against that claim.  It also carries the
       only reading of the state-dependent judgement that a deck can show:
       failed-but-shut counts facets that are dead AND in compression, which
       still conduct. */
    printf("[LOADPATH PHANTOM] inc=%" ITGFORMAT " time=%.12e is %"
           ITGFORMAT " increment(s) past severance at inc=%" ITGFORMAT
           " (theta=%.7f); this is not a specimen; joined by %" ITGFORMAT
           " open-failed and %" ITGFORMAT " failed-but-shut facet(s), "
           "connected=%" ITGFORMAT "\n",
           inc,t,lp->ninc_past,lp->sev_inc,lp->sev_time,
           lp->nfacetdead,lp->nfacetshut,lp->connected);
    fflush(stdout);
  }
}

/* Print the verdict however the run ends.

   A run that dies at a wall leaves through FORTRAN stop(), which is
   exit(201), and never reaches the bottom of nonlingeo() - so the case that
   most needs the verdict was the one case that did not print it.  Measured:
   fast-wrapped-wall reports "never severed" from its summary while its own
   log carries [LOADPATH SEVERED] inc=96.

   atexit covers every exit path, including the ones not enumerated here. */
static loadpath *lp_exit_target=NULL;
static ITG lp_exit_done=0;

static void loadpath_atexit(void)
{
  if((lp_exit_target!=NULL)&&(lp_exit_done==0)){
    lp_exit_done=1;
    loadpath_summary(lp_exit_target);
    fflush(stdout);
  }
}

void loadpath_report_at_exit(loadpath *lp)
{
  if(lp_exit_target==NULL) atexit(loadpath_atexit);
  lp_exit_target=lp;
}

/* The one-line verdict for the end of a run. */
void loadpath_summary(const loadpath *lp)
{
  if(!lp->armed){
    printf("[LOADPATH] not armed; this run made no statement about "
           "whether the specimen stayed in one piece\n");
    return;
  }
  if(!lp->sev_seen){
    /* Every run states its own final topological condition, not only the
       ones that end badly.  "still connected, 96 facets of which 51 have
       failed" is the reading that says how close a run that stopped for
       some other reason was to coming apart. */
    printf("[LOADPATH] the specimen kept a grip-to-grip load path for the "
           "whole run (%s); at the last converged increment %" ITGFORMAT
           " element(s) were live and %" ITGFORMAT " of %" ITGFORMAT
           " cohesive facet(s) had fully failed\n",
           lp->origin,lp->nlive,lp->nfacetdead,lp->nfacet);
    return;
  }
  if(lp->ninc_past==0){
    printf("[LOADPATH] SEVERED at increment %" ITGFORMAT ", theta=%.7f, and "
           "the run ended there.  %" ITGFORMAT " of %" ITGFORMAT " cohesive "
           "facet(s) had fully failed; nothing past this point would have "
           "been a specimen.\n",
           lp->sev_inc,lp->sev_time,lp->nfacetdead,lp->nfacet);
    return;
  }
  printf("[LOADPATH] SEVERED at increment %" ITGFORMAT ", theta=%.7f, and "
         "the run then continued for %" ITGFORMAT " further increment(s).  "
         "Those increments are two separated pieces joined by %" ITGFORMAT
         " failed cohesive facet(s) at the residual-stiffness floor; any "
         "load they report is that floor, not the specimen.\n",
         lp->sev_inc,lp->sev_time,lp->ninc_past,lp->nfacetdead);
}

ITG loadpath_severed(const loadpath *lp)
{
  return (lp->armed&&lp->sev_seen)?1:0;
}

/* ---------------------------------------------------------------- tests */

static void lp_chk(const char *name,ITG got,ITG want,ITG *nbad)
{
  ITG ok=(got==want);
  printf("   %-40s got=%-4" ITGFORMAT " want=%-4" ITGFORMAT " %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
}

/* A three-link chain.  Element 0 and element 2 are C3D4 and share no node
   with each other; the middle link is what joins the grips, and every test
   below changes only the middle link.

       grip A {1,2} -- e0(1,2,3,4) -- [middle] -- e2(5,6,7,8) -- grip B {7,8}

   The middle link is either a C3D4 on nodes (3,4,5,6) or a UC6 facet on
   (3,4,9,5,6,10).  Nothing here reimplements the walk: the kernel under test
   is the same damconnect the run uses. */
static ITG lp_chain(ITG middle_is_facet,ITG middle_deleted,ITG nfailed,
                    ITG *iconn_out)
{
  ITG ipkon[3],kon[16],ne=3,nk=10;
  ITG iconn,i,nstate=5,mi0=3;
  double xstate[5*3*3];
  char lakon[24];
  loadpath lp;

  for(i=0;i<24;i++) lakon[i]=' ';
  for(i=0;i<16;i++) kon[i]=0;
  for(i=0;i<5*3*3;i++) xstate[i]=0.;

  /* e0 */
  ipkon[0]=0; kon[0]=1;kon[1]=2;kon[2]=3;kon[3]=4;
  memcpy(&lakon[0],"C3D4",4);

  if(middle_is_facet){
    ipkon[1]=4; kon[4]=3;kon[5]=4;kon[6]=9;kon[7]=5;kon[8]=6;kon[9]=10;
    lakon[8]='U'; lakon[8+6]=(char)3; lakon[8+7]=(char)6;
    ipkon[2]=10; kon[10]=5;kon[11]=6;kon[12]=7;kon[13]=8;
  }else{
    ipkon[1]=4; kon[4]=3;kon[5]=4;kon[6]=5;kon[7]=6;
    memcpy(&lakon[8],"C3D4",4);
    ipkon[2]=8; kon[8]=5;kon[9]=6;kon[10]=7;kon[11]=8;
  }
  memcpy(&lakon[16],"C3D4",4);
  if(middle_deleted) ipkon[1]=-1;

  /* mark nfailed of the middle facet's three points as failed */
  for(i=0;i<nfailed;i++) xstate[3+nstate*(i+mi0*1)]=1.;

  loadpath_init(&lp);
  lp.nk=nk; lp.imode=0; lp.armed=1;
  NNEW(lp.nodesa,ITG,2); NNEW(lp.nodesb,ITG,2);
  lp.nodesa[0]=1;lp.nodesa[1]=2; lp.na=2;
  lp.nodesb[0]=7;lp.nodesb[1]=8; lp.nb=2;
  iconn=loadpath_census(&lp,ipkon,kon,lakon,&ne,xstate,nstate,mi0,NULL);
  loadpath_free(&lp);
  *iconn_out=iconn;
  return iconn;
}

ITG loadpath_selftest(void)
{
  ITG nbad=0,iconn;

  printf("[LOADPATH] self test%s","\n");

  /* A-B: the walk itself, with no facet anywhere */
  lp_chain(0,0,0,&iconn); lp_chk("A intact bulk chain is connected",iconn,1,&nbad);
  lp_chain(0,1,0,&iconn); lp_chk("B deleted middle severs it",iconn,0,&nbad);

  /* C-E: the facet judgement.  E is the one that matters - a facet that is
     partly failed is NOT failed, and calling it dead would sever a specimen
     that is still carrying load. */
  lp_chain(1,0,0,&iconn); lp_chk("C live facet conducts",iconn,1,&nbad);
  lp_chain(1,0,3,&iconn); lp_chk("D fully failed facet does not",iconn,0,&nbad);
  lp_chain(1,0,2,&iconn); lp_chk("E 2 of 3 points failed still conducts",
                                 iconn,1,&nbad);
  lp_chain(1,1,0,&iconn); lp_chk("F deleted facet severs it",iconn,0,&nbad);

  /* G-J: deriving the grips from the deck's own boundary cards.  The layout
     is the one every deck in this tree uses:
        FACE_X0  dir 1 = 0      reacting
        FACE_XL  dir 1 = 1.0    driven
        FIXPOINT dir 2,3 = 0    rigid-body suppressors, NOT a grip */
  {
    ITG nodeboun[8]={1,2, 3,4, 5,5, 6,6};
    ITG ndirboun[8]={1,1, 1,1, 2,3, 2,3};
    double xboun[8]={0.,0., 1.,1., 0.,0., 0.,0.};
    ITG *a=NULL,*b=NULL,na=0,nb=0,idir;

    idir=loadpath_derive(nodeboun,ndirboun,xboun,8,10,&a,&na,&b,&nb);
    lp_chk("G driven direction is 1",idir,1,&nbad);
    lp_chk("G reacting set has 2 nodes",na,2,&nbad);
    lp_chk("G driven set has 2 nodes",nb,2,&nbad);
    lp_chk("G the dir-2/3 fixpoints are not a grip",
           (na==2)&&(a[0]==1)&&(a[1]==2),1,&nbad);
    lp_chk("G the driven set is the nonzero one",
           (nb==2)&&(b[0]==3)&&(b[1]==4),1,&nbad);
    if(a!=NULL) SFREE(a);
    if(b!=NULL) SFREE(b);
  }
  {
    /* H: nothing driven anywhere - the deck presents no grip pair and the
       module must refuse to arm rather than invent one. */
    ITG nodeboun[4]={1,2,3,4};
    ITG ndirboun[4]={1,1,2,3};
    double xboun[4]={0.,0.,0.,0.};
    ITG *a=NULL,*b=NULL,na=0,nb=0,idir;
    idir=loadpath_derive(nodeboun,ndirboun,xboun,4,10,&a,&na,&b,&nb);
    lp_chk("H no driven dof -> refuses to arm",idir,0,&nbad);
    if(a!=NULL) SFREE(a);
    if(b!=NULL) SFREE(b);
  }
  {
    /* I: a negative prescribed value is driven too - the closing step of the
       close benchmark drives NRIGHT to -0.12, and reading the sign instead
       of the magnitude would put the driven face on the reacting side and
       make the two lists share a node. */
    ITG nodeboun[3]={1,2,2};
    ITG ndirboun[3]={1,1,1};
    double xboun[3]={0.,-0.12,-0.12};
    ITG *a=NULL,*b=NULL,na=0,nb=0,idir;
    idir=loadpath_derive(nodeboun,ndirboun,xboun,3,10,&a,&na,&b,&nb);
    lp_chk("I a negative prescribed value is driven",idir,1,&nbad);
    lp_chk("I and counted once",nb,1,&nbad);
    lp_chk("I reacting side keeps only the zero node",
           (na==1)&&(a[0]==1),1,&nbad);
    if(a!=NULL) SFREE(a);
    if(b!=NULL) SFREE(b);
  }
  {
    /* J: a node prescribed both zero and nonzero in the same direction must
       not end up in both lists, which would read as trivially connected for
       ever. */
    ITG nodeboun[2]={1,1};
    ITG ndirboun[2]={1,1};
    double xboun[2]={0.,1.};
    ITG *a=NULL,*b=NULL,na=0,nb=0,idir;
    idir=loadpath_derive(nodeboun,ndirboun,xboun,2,10,&a,&na,&b,&nb);
    lp_chk("J over-specified node does not join both sides",
           (idir==0)||(na==0),1,&nbad);
    if(a!=NULL) SFREE(a);
    if(b!=NULL) SFREE(b);
  }

  /* L: a failed facet in COMPRESSION is still a load path.  This is the
     judgement close.inp exists to check, and it used to be a policy switch
     instead of a measurement. */
  {
    ITG mi0=3,e=1,j;
    double stx[6*3*4];
    for(j=0;j<6*3*4;j++) stx[j]=0.;
    lp_chk("L all points at zero traction -> open",
           loadpath_facet_open(stx,mi0,e,3),1,&nbad);
    stx[6*(mi0*e+1)]=+5.;                     /* one point in tension */
    lp_chk("L tension does not close it",
           loadpath_facet_open(stx,mi0,e,3),1,&nbad);
    stx[6*(mi0*e+2)]=-1.e-30;                 /* one point in compression */
    lp_chk("L one compressed point conducts",
           loadpath_facet_open(stx,mi0,e,3),0,&nbad);
    lp_chk("L a neighbouring facet is unaffected",
           loadpath_facet_open(stx,mi0,0,3),1,&nbad);
    lp_chk("L no traction available -> open, as before",
           loadpath_facet_open(NULL,mi0,e,3),1,&nbad);
  }

  /* M: hysteresis.  A loss of path that does not persist is not severance,
     and the increment recorded is the FIRST one, not the confirming one. */
  {
    loadpath lp; ITG w;
    loadpath_init(&lp); lp.armed=1;
    lp_chk("M default confirmation is three",lp.nconfirm,3,&nbad);
    lp.connected=0; w=loadpath_latch(&lp,20,0.2);
    lp_chk("M first loss is provisional",w,3,&nbad);
    lp.connected=1; w=loadpath_latch(&lp,21,0.3);
    lp_chk("M reconnecting clears it",w,0,&nbad);
    lp_chk("M and resets the count",lp.ndisc,0,&nbad);
    lp.connected=0;
    loadpath_latch(&lp,30,0.4); loadpath_latch(&lp,31,0.5);
    lp_chk("M still provisional at two",lp.sev_seen,0,&nbad);
    w=loadpath_latch(&lp,32,0.6);
    lp_chk("M confirmed at three",w,1,&nbad);
    lp_chk("M records the FIRST increment",lp.sev_inc,30,&nbad);
    loadpath_free(&lp);
  }

  /* K: the latch.  Severance must not un-happen when a failed facet is
     driven back into compression and starts carrying load again. */
  {
    loadpath lp; ITG w;
    loadpath_init(&lp); lp.armed=1;
    lp.connected=1; w=loadpath_latch(&lp,10,0.1);
    lp_chk("K intact reports nothing",w,0,&nbad);
    lp.connected=0;                       /* confirm it: hysteresis */
    loadpath_latch(&lp,11,0.2); loadpath_latch(&lp,12,0.3);
    w=loadpath_latch(&lp,13,0.4);
    lp_chk("K severance is announced once",w,1,&nbad);
    lp_chk("K and records the first increment",lp.sev_inc,11,&nbad);
    lp.connected=1;                       /* the facet closes again */
    w=loadpath_latch(&lp,14,0.5);
    lp_chk("K reconnection does not un-sever",w,2,&nbad);
    lp_chk("K past-severance increments are counted",lp.ninc_past,1,&nbad);
    lp_chk("K severed stays severed",loadpath_severed(&lp),1,&nbad);
    loadpath_free(&lp);
  }

  printf("[LOADPATH] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
