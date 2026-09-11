/*     CalculiX - damage/fracture extension                              */
/*     rescuelevel.c: THE judgement of which rescue level can act first. */

/* Why this module exists
   ----------------------
   The rescue ladder is sequenced by one integer, `damage_rescue_used`,
   incremented once per wall.  Its levels are not interchangeable:

     level 1  emergency backtracking   acts only when idamagereeq==1
     level 2  event step               acts only when idamagereeq==1
     level 3  regularised diagonal, IF CCX_DAMAGE_REEQ_RESCUE3 armed it;
              otherwise the dogleg trust region
     level 4  coupled continuation

   A wall reached with `idamagereeq==0` never enters a same-load solve, so
   levels 1 and 2 cannot act on it and simply recompute the identical attempt.
   nonlingeo.c already knew this and already had a skip - but the skip was
   gated on `damage_reg_nlam>0`, i.e. it only fired when the REGULARISATION
   ladder was armed, and it jumped to level 2.

   `test/s3rad/run_s3rad.sh` arms RESCUE2 + TR_DOGLEG and not RESCUE3, so
   `damage_reg_nlam` is 0, the skip never fires, and the level it would have
   jumped to could not act either.  Measured on s3rad, increment 602, three
   consecutive walls on one increment:

     wall 1 (level 1)   0 [DAMAGE BT], 0 [DAMAGE TR]   nothing happened
     wall 2 (level 2)   0 [DAMAGE BT], 0 [DAMAGE TR]   nothing happened
     wall 3 (level 3)  12 [DAMAGE TR]                  the dogleg CONVERGED

   Two of the three attempts were spent recomputing an identical solve, and
   the whole run printed `skipping straight` exactly zero times.  That is
   05-DEBT.md item 3's "three consecutive attempts produced bit-identical
   residual sequences", still live on the configuration the target runner
   uses.

   This module answers one question - given what is armed and how the wall
   was reached, which level is the first that can DO anything - and has a self
   test, because the previous version of this logic was a condition nobody
   could check.  It changes no equation and takes no decision of its own.   */

#include <stdio.h>
#include <stdlib.h>
#include "CalculiX.h"

/* Returns the first level that can act, given:
     reeq     idamagereeq: 1 if a same-load solve is available
     nlam     damage_reg_nlam: regularisation ladder length, 0 if not armed
     dl       damage_dl_mode: 1 if the dogleg is armed
     maxlevel damage_rescue_maxlevel
   Never returns more than maxlevel, and never less than `used`+1's natural
   value of 1: a caller that is already past a level must not be sent back. */
ITG rescuelevel_first_useful(ITG reeq,ITG nlam,ITG dl,ITG maxlevel)
{
  if(reeq!=0) return 1;                 /* 1 and 2 can act; start at 1 */
  if(nlam>0)  return (maxlevel>=2)?2:maxlevel;   /* regularised level 2 */
  if(dl==1)   return (maxlevel>=3)?3:maxlevel;   /* the dogleg          */
  return 1;                             /* nothing above can act either */
}

/* ---------------------------------------------------------------- tests */

static ITG rl_chk(const char *name,ITG got,ITG want,ITG *nbad)
{
  ITG ok=(got==want);
  printf("   %-46s got=%-3" ITGFORMAT " want=%-3" ITGFORMAT " %s%s",
         name,got,want,ok?"ok":"*** FAIL ***","\n");
  if(!ok) (*nbad)++;
  return ok;
}

ITG rescuelevel_selftest(void)
{
  ITG nbad=0;

  printf("[RESCUELEVEL] self test%s","\n");

  /* A: with a same-load solve available, level 1 is useful - no skipping */
  rl_chk("A reeq=1 starts at level 1",
         rescuelevel_first_useful(1,0,1,3),1,&nbad);
  rl_chk("A reeq=1 even with the ladder armed",
         rescuelevel_first_useful(1,4,1,4),1,&nbad);

  /* B: the case the old condition handled - regularisation armed */
  rl_chk("B reeq=0 with regularisation -> level 2",
         rescuelevel_first_useful(0,4,0,3),2,&nbad);

  /* C: THE DEFECT.  reeq=0, no regularisation, dogleg armed: the old
     condition fired not at all and burned levels 1 and 2. */
  rl_chk("C reeq=0, no reg, dogleg -> level 3",
         rescuelevel_first_useful(0,0,1,3),3,&nbad);

  /* D: nothing above can act either, so do not invent a level */
  rl_chk("D reeq=0, nothing armed -> level 1",
         rescuelevel_first_useful(0,0,0,3),1,&nbad);

  /* E: never exceed the ladder that is actually configured */
  rl_chk("E dogleg wanted but maxlevel 2",
         rescuelevel_first_useful(0,0,1,2),2,&nbad);
  rl_chk("E regularisation wanted but maxlevel 1",
         rescuelevel_first_useful(0,4,0,1),1,&nbad);

  printf("[RESCUELEVEL] self test %s (%" ITGFORMAT " failure(s))%s",
         nbad?"FAILED":"PASSED",nbad,"\n");
  return nbad;
}
