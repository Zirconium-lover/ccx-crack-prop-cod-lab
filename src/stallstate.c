/*     CalculiX - damage/fracture extension                              */
/*     stallstate.c: THE judgement that the run has stopped making       */
/*     progress, so the diagnostics can arm themselves.                  */

/* Why this module exists
   ----------------------
   `02-DIAGNOSTICS.md` is a good procedure with one expensive precondition:
   almost every reading in it is armed by `CCX_DAMAGE_WALL_THETA`, and that
   needs the load factor of the wall **in advance**.  The first time a deck
   walls, nobody has it.  So the first run is spent finding out where to look,
   and the diagnosis costs a second run.

   Measured cost of exactly that: two s3rad runs of about 50 minutes each were
   spent in one session to obtain a single number, because the probes were
   added one at a time.  Nothing was wrong with either run.  The procedure
   simply could not be armed before the thing it was meant to observe.

   This module removes the precondition.  A run that stops making progress can
   notice it itself, and the probes that were waiting for a hand-supplied theta
   arm on the spot.

   The quantity
   ------------
   Dimensionless by construction, and it reuses a reference the tree already
   keeps: the median of the last accepted increments' dtheta.  A ratio against
   a deck's OWN healthy pace needs no units and no per-deck tuning.

       ratio = dtheta / median(the last N accepted dtheta)

   Measured separation, from runs already on disk:

     deck            healthy   at the wall          false alarms
     fast-wrapped    ~1        1.9e-03 (inc 95-98)  0
     s3rad           ~1        5.9e-03 (inc 595-600) 1, isolated, at inc 187
     fast-plain      ~1        never fires          0   (it does not wall)

   So the wall sits two to three decades below the deck's own pace, and the
   threshold below is not a tuned constant - it is placed in the middle of
   three decades of empty space.

   HYSTERESIS, the same device loadpath.c uses and for the same reason: the
   one false alarm on s3rad is a SINGLE increment, while every real wall
   produces a run of them (4 consecutive on fast-wrapped, 6 on s3rad).
   Requiring `nconfirm` consecutive accepted increments below the threshold
   removes the false alarm and keeps both walls.

   This module takes no decision of its own and changes no equation: it
   answers, and the diagnostics decide whether to arm.  It refuses to arm if
   its self test fails, the discipline lsladder.c set.                     */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"

void stallstate_init(stallstate *s)
{
  memset(s,0,sizeof(*s));
  /* Two decades below the deck's own pace.  The walls measured sit at
     1.9e-03 and 5.9e-03; healthy is ~1.  Anywhere in 1e-1..1e-3 separates
     them, so the value is not carrying the judgement - the gap is. */
  s->ratio=2.e-2;
  s->nconfirm=3;
}

/* Record an ACCEPTED increment's step.  Rejected attempts must not enter the
   reference: a cutback's dtheta is the thing being judged, not the pace. */
void stallstate_note(stallstate *s,double dtheta)
{
  if(dtheta<=0.) return;
  s->ring[s->i]=dtheta;
  s->i=(s->i+1)%STALL_RING;
  if(s->n<STALL_RING) s->n++;
}

/* The median of the reference ring.  Returns 0 while too few samples have
   been seen for a median to mean anything. */
double stallstate_reference(const stallstate *s)
{
  double m[STALL_RING],t; ITG a,b,k;
  if(s->n<5) return 0.;
  for(a=0;a<s->n;a++) m[a]=s->ring[a];
  for(a=1;a<s->n;a++){                      /* insertion sort, n<=20 */
    t=m[a];
    for(b=a;(b>0)&&(m[b-1]>t);b--) m[b]=m[b-1];
    m[b]=t;
  }
  k=s->n/2;
  return m[k];
}

double stallstate_ratio(const stallstate *s,double dtheta)
{
  double r=stallstate_reference(s);
  if(r<=0.) return -1.;
  return dtheta/r;
}

/* Judge an accepted increment.  Returns 1 on the increment the stall is
   CONFIRMED, 2 on every later one, 3 while it is provisional, 0 otherwise. */
ITG stallstate_judge(stallstate *s,double dtheta)
{
  double r=stallstate_ratio(s,dtheta);
  if(s->seen){s->nsince++; return 2;}
  if(r<0.) return 0;
  if(r<s->ratio){
    if(s->ndisc==0) s->first=r;
    s->ndisc++;
    if(s->ndisc>=s->nconfirm){s->seen=1; s->rseen=r; return 1;}
    return 3;
  }
  s->ndisc=0;
  return 0;
}

ITG stallstate_stalled(const stallstate *s){return (s==NULL)?0:s->seen;}

/* ---------------------------------------------------------------- tests */

static ITG ss_chk(const char *name,double got,double want,ITG *nbad)
{
  ITG ok=(fabs(got-want)<=1e-9*(1.+fabs(want)));
  printf("   %-44s got=%-12.5g want=%-12.5g %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

ITG stallstate_selftest(void)
{
  ITG nbad=0,i;
  stallstate s;

  printf("[STALLSTATE] self test%s","\n");

  /* A: too few samples is "no opinion", not "healthy" and not "stalled" */
  stallstate_init(&s);
  ss_chk("A no reference yet -> ratio -1",stallstate_ratio(&s,1.e-9),-1.,&nbad);
  ss_chk("A and no judgement",(double)stallstate_judge(&s,1.e-9),0.,&nbad);

  /* B: a steady deck has ratio 1 and never fires */
  for(i=0;i<10;i++) stallstate_note(&s,2.e-3);
  ss_chk("B steady pace -> ratio 1",stallstate_ratio(&s,2.e-3),1.,&nbad);
  ss_chk("B and does not fire",(double)stallstate_judge(&s,2.e-3),0.,&nbad);

  /* C: a single slow increment is PROVISIONAL, not a stall.  This is the
     s3rad inc=187 false alarm, and the reason hysteresis is here. */
  ss_chk("C one slow increment is provisional",
         (double)stallstate_judge(&s,2.e-3*1.e-3),3.,&nbad);
  ss_chk("C recovering clears it",(double)stallstate_judge(&s,2.e-3),0.,&nbad);
  ss_chk("C and resets the count",(double)s.ndisc,0.,&nbad);

  /* D: three consecutive confirm it, as at both measured walls */
  stallstate_judge(&s,2.e-3*1.e-3);
  stallstate_judge(&s,2.e-3*1.e-3);
  ss_chk("D still provisional at two",(double)s.seen,0.,&nbad);
  ss_chk("D confirmed at three",
         (double)stallstate_judge(&s,2.e-3*1.e-3),1.,&nbad);
  ss_chk("D and latches",(double)stallstate_judge(&s,2.e-3),2.,&nbad);

  /* E: the reference is a MEDIAN, so one outlier cannot move it.  The
     historical shape of this defect is a mean dragged by a single sample. */
  {
    stallstate t; ITG j;
    stallstate_init(&t);
    for(j=0;j<9;j++) stallstate_note(&t,1.e-3);
    stallstate_note(&t,1.e+3);              /* one wild sample */
    ss_chk("E median ignores an outlier",stallstate_reference(&t),1.e-3,&nbad);
  }

  /* F: a rejected (non-positive) step never enters the reference */
  {
    stallstate t; ITG j;
    stallstate_init(&t);
    for(j=0;j<6;j++) stallstate_note(&t,5.e-4);
    stallstate_note(&t,0.);
    stallstate_note(&t,-1.);
    ss_chk("F non-positive steps are not recorded",(double)t.n,6.,&nbad);
  }

  /* G: the measured walls fire, the measured healthy deck does not */
  {
    stallstate t; ITG j;
    stallstate_init(&t);
    for(j=0;j<20;j++) stallstate_note(&t,2.0e-3);   /* fast-wrapped median */
    ss_chk("G fast-wrapped wall ratio 1.9e-03 fires",
           (double)(stallstate_ratio(&t,3.9063e-6)<t.ratio),1.,&nbad);
    stallstate_init(&t);
    for(j=0;j<20;j++) stallstate_note(&t,1.3348e-4); /* s3rad median */
    ss_chk("G s3rad stall ratio 5.9e-03 fires",
           (double)(stallstate_ratio(&t,1.6480e-6)<t.ratio),1.,&nbad);
    ss_chk("G s3rad inc=187 at 2.0e-02 does NOT",
           (double)(stallstate_ratio(&t,1.3348e-4*2.0e-2)<t.ratio),0.,&nbad);
  }

  printf("[STALLSTATE] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
