# 1. Where the code actually is

Everything here is measured on this branch unless marked otherwise. Numbers
without a source are not in this document.

## Shape

| | |
|---|---|
| `src/nonlingeo.c` | 16140 lines; the `nonlingeo()` body alone is ~14381 |
| `CCX_*` switches the binary reads | **139** |
| of those, explained nowhere but the line that reads them | **63** |
| of those, set by **no test in this tree** | **122** |
| units already extracted with an owner and a self test | `src/lsladder.c`, `src/damstate.c`, `src/damswitch.c`, `src/loadpath.c` |
| units still inline in `nonlingeo.c` | erosion/topology, convergence judgement, the solve, every diagnostic |
| the gate | **11 cases, 47 s**, proven able to go red |

The switch counts are generated, not estimated: `tools/mkswitches.py` scans
the sources, `docs/SWITCHES.md` is its output, and the coverage column comes
from the gate reading back its own `[SWITCHES]` banners.

## The mechanisms that have accumulated

Six overlapping ways of getting an increment to converge, added one wall at a
time: the adaptive damage line-search ladder (`lsladder.c`), transactional
backtracking, Rescue level 1, Rescue level 2 with an event step, a dogleg
trust region as level 3, and dissipation/crack-control path following. Plus
the load-path judgement (AUTOSPC, and its extension to the force residual),
the terminal-deletion machinery, a minimum-increment override, and a
reference-force floor.

Each was justified by the wall it passed. None was removed when the next wall
showed the previous diagnosis had been incomplete. That is the accumulation.

## What is settled, with evidence

**The walls are not one thing.** On the fast wrapped deck the load-path mask
changes nothing and the crack-face kink changes everything; on `s3rad` the
opposite holds. Calling them all "the wall" is part of why fixing one kept
producing the next.

**Second `s3rad` wall — increment 554, `theta=0.2555742`.** The peak force
residual sat on node 1246: a fragment held by ONE live bulk element with all
six of its cohesive facets failed and open. Equilibrating it needs a
displacement of order one against a grip that has moved 0.2556, so its
residual fell 0.1% per Newton iteration. That component was inherited by
every following increment and grew — 15.7, 50.4, 74.9, 94.5, 105.7% of
tolerance over increments 551–555 — **while `dtheta` collapsed 48x**. Cutting
the step made equilibrium worse, because the part that will not reduce is not
proportional to the step. `CCX_DAMAGE_AUTOSPC_FORCE=1` extends the load-path
judgement to the force residual: increment 554 → 930, `theta` 0.2556 →
0.5575, `dtime` recovered 1000x, pre-wall history bit-identical over all 1147
attempts.

**The crack-face closure kink.** `cohesive_uc6.f` leaves the normal traction
continuous at `deltal(1)=0` but jumps its slope from `g*kn` to `kn` — a
factor of `1/gmin`. Verified directly on one facet, not argued from reading
the source: far-compression `dT/dd = 999999` against `Kn=1e6`, far-tension
`1` against `g*Kn=1`, **ratio 1e+06**. `CCX_UC6_CONTACT_SMOOTH` blends the
two slopes over a penetration band. On the fast wrapped deck the wall goes
(rc 201→0, increment 99→515, active-set transitions 19→0 at every ladder
rung, `|r0|2` 3.48e-02→1.89e-10) and the fracture does not move: the same 64
elements, same order, largest shift in a deletion time 7.4e-04 in `theta`.

**The regulariser is a byte-for-byte no-op where no damaged facet closes.**
Plain deck, both arms: 542 attempts, 1283 Newton iterations, `m.sta` and
`m.damage` equal byte for byte. The reason is exact — the blended expression
reduces to `kn*d` identically at `g=1`, so the two laws differ only where a
facet is both damaged and in compression.

**On `s3rad` the regulariser moves nothing physical.** Severance `theta`
0.3411981 (sharp) against 0.3418522 (blended) — **0.19%**. The 28-increment
difference in stopping point is entirely inside the post-severance regime and
means nothing. This is a *confirmation that it is harmless at full scale*,
not a demonstration that it helps there.

**One owner for the load-path judgement.** `src/damstate.c`, verified
bit-identical to the pre-refactor binary on a deck where the predicate
actually decides (14 increments), and over 849 attempts / 3135 committed
deletions of `s3rad`.

**Nothing in `s3rad` is adrift.** The per-node diagonal test is blind by
construction to a piece that is internally stiff but attached to nothing.
Measured with `test/s3rad/fragments.py`: at every facet-stiffness threshold,
**0 elements adrift**. The free-fragment class does not occur here.

**Runs are not reproducible across thread counts.** `MKL_CBWR=COMPATIBLE`
fixes reproducibility across instruction sets, not across thread counts. Two
runs identical for 482 attempts diverged at 483 on a 6-against-5 iteration
count; twenty increments later they were on different walls.

**Feature off is byte-identical on the TARGET deck, not just the fast ones.**
A single-variable A/B — pre-change binary (`225439c`, no `loadpath.c` in it at
all) against HEAD, same committed deck, same `CCX_DAMAGE_AUTOSPC_FORCE=1`,
4 threads both arms — is **byte-identical across all 1253 status lines** of
`m.sta`, at 42807 elements, up to the point both arms stop. This is the
strongest form of the feature-off obligation and it now holds where it
matters.

**The second wall does not yield at 4 threads, and this is not a regression.**
Both arms stall at the *same* state: increment 604, `theta = 0.255032`,
attempt 3U, `dtime = 1.23596e-06`, no new status line in 120 s. The recorded
`CCX_DAMAGE_AUTOSPC_FORCE=1` result — increment 554 → 930, `theta` 0.2556 →
0.5575 — was measured at a different thread count, and `05-DEBT.md` §5 says
exactly what to expect from that. **So `s3rad` severance was not reached in
this session**, and every full-scale statement here stops at `theta=0.255`.

**Correction to the recorded second-wall result: `AUTOSPC_FORCE`'s 554 -> 930
was measured on TWO threads and is not robust.** `run_s3rad.sh` defaults to
`OMP_NUM_THREADS=6` but the recorded gain was taken at two. Re-run at **four**,
same deck, same switch, same PARDISO, the run does not pass the second wall at
all: it stalls at **increment 604, `theta = 0.255032`**, `dtime` 1.24e-06, no
new status line in 120 s. Both the pre-change and post-change binaries stall
there identically, so this is a property of the mechanism and not of any
patch. The entry in this file that reads "increment 554 -> 930, `theta` 0.2556
-> 0.5575" should be read as **"at two threads"**, and `05-DEBT.md` §5 is the
reason: this is exactly the thread-count non-reproducibility it warns about,
showing up in a result the branch had been treating as settled.

**What the §4 signal actually points at.** Node 1246 is a fragment held by ONE
live bulk element with all six of its cohesive facets failed **and open**. So:

- **removing dead facets cannot help it.** Its diagonal is held by the
  tetrahedron, not by the facets — the root fix of `05-DEBT.md` §1
  (`CCX_DAMAGE_FACET_DELETE`) would change nothing for this node;
- what the signal points at instead is two of the judgements `06-TARGET.md`
  lists as having **no owner**: *what counts as eroded* (should a nearly dead
  element holding a fragment be removed?) and *what counts as converged*
  (should a node whose `need_du` exceeds a physical length be judged at all —
  `NEXT_TASK.md` already proposes exactly this).

`AUTOSPC_FORCE` is a third answer to the same question, given at the level of
the norm rather than the judgement, which is why it excludes 24.4x tolerance
and still does not converge. Fixing the norm cannot fix a state the model
should not be in.

**The exclusion report is showing §4's "this fix has become a fiction"
signal.** At that stall, in BOTH arms identically:

| | |
|---|---|
| node | 1246 — the node the second wall was always about |
| excluded residual | `4.4835e-02` |
| residual actually judged (`ram[0]`) | `4.1680e-02` |
| tolerance | `1.8355e-03` (`0.1221 x qam`) |
| **excluded / tolerance** | **24.4x** |

`02-DIAGNOSTICS.md` §4 records the healthy peak as **1.7x**, returning to
machine zero once the fragment resolved. Here it is 24.4x, it is *larger than
the residual being judged*, and it is not returning to zero. Earlier in the
same run it does return to machine zero (`1.3e-13`, `6.4e-14`), so the
mechanism is not broken everywhere — it has stopped discriminating at this
state. By the criterion the project wrote for itself, `CCX_DAMAGE_AUTOSPC_FORCE`
must be revisited before any further weight is put on it.

**The grips can be read off the deck.** `loadpath.c` derives them from the
`*BOUNDARY` cards — the direction carrying the largest prescribed magnitude,
driven nodes against nodes held at zero in that same direction. On the target
deck that gives **181 reacting against 180 driven**, which is exactly
`FACE_X0_NSET` (181) and `FACE_XL_NSET` (180), with `FIXPOINTA` and
`FIXPOINTB` correctly excluded because they hold directions 2 and 3. Measured
on `m12_s3rad_gc24_w.inp` itself, not inferred from the fast decks. So the
judgement no longer depends on anybody having configured
`CCX_FRACTURE_TERMINATION`; that switch is now an override, not the arming
condition.

**A specimen now says when it has stopped being one.** `src/loadpath.c` owns
the judgement, asks it once per converged increment, derives the grips from
the deck's own `*BOUNDARY` cards when nobody names them, and stops the run at
severance. The three gaps it closed — a hand-named set pair, a second opt-in
switch, and a call site inside the bulk-deletion transaction — are in
`05-DEBT.md` §1 with the measurements.

**The phantom regime was worse than "wasted increments".** On the fast plain
deck the specimen separates at increment 65, `theta=0.1175`, and the run used
to continue to 507. Over that stretch the grip reaction is a **straight line
through the origin**, `Fx/theta` constant to **0.04%** across `theta`
0.19…1.0: two detached halves joined by 36 failed facets at `gmin*Kn`. The
reported load falls to 1.2% of peak at severance and climbs back to **6.5% of
peak** by `theta=1`. The model did not merely keep going, it reported a
recovering specimen.

**~~The fast-wrapped wall is past severance.~~ WITHDRAWN — it was an artefact
of a load-path judgement that ignored compression.** The claim was that the
deck severs at increment 96 and walls three increments later, so its wall and
the regularised arm's `theta=1` were phantom. Making the judgement
state-dependent removes the severance entirely: **`fast-wrapped` never loses
its grip-to-grip load path**, and its wall at increment 99, `theta=0.158766`,
happens on a connected specimen. Measured per increment:

```
inc=88   52 open-failed,  2 failed-but-shut   connected=1
inc=93   53 open-failed,  1 failed-but-shut   connected=1
inc=94   54 open-failed,  0 failed-but-shut   connected=1
inc=96   53 open-failed,  1 failed-but-shut   connected=1   (live 2193 -> 2192)
```

`open-failed + failed-but-shut` is 54 at every one of those increments, so the
classification is complete and consistent rather than noise. At increment 96 a
bulk element is deleted; the old judgement lost the path there, and the one
failed-but-compressed facet carries it. That single facet is the whole
difference between "the regulariser was validated against a phantom" and "it
was validated against a specimen".

The lesson is the one this branch keeps re-learning: **a wall called phantom
on a judgement nobody had measured is not a finding.** Two of the three walls
this work called phantom — this one and the second `s3rad` wall — turned out
to be real once the judgement was correct. Only the fast-plain 442 increments
survive, and that deck has `0 failed-but-shut` at every census.

**The census fires exactly when it should, checked against the deck's own
state output.** On `mixed.inp` — 120 bulk elements in two halves joined by two
UC6 facets and nothing else — the `.dat` internal-state block shows all six
integration points carrying the failure flag `0.0` at increment 2499
(`theta=0.4998`) and `1.0` at increment 2500 (`theta=0.5000`). The census
reports severance at increment 2500. No off-by-one, no premature firing, and
the check is against `*EL PRINT` output rather than against the same code
path that makes the judgement.

**The second `s3rad` wall is NOT past severance — the third one is.** The
census now separates the two classes at full scale. Stock `run_s3rad.sh`
configuration, 4 threads, PARDISO, HEAD binary: the run reaches the recorded
second wall at `theta = 0.255023` against the recorded `0.2555742` (0.2%
apart, which is what a different thread count buys), and at the last census
before it — increment 550, `theta = 0.24175` — the specimen is **still one
piece**: `connected=1`, 1711 of 5318 facets failed, 39148 elements live. The
grip reaction there is 66.4 against a peak of 3112.8 at `theta=0.1575`, i.e.
**2.1% of peak**.

So the second wall is a real convergence failure of a real specimen, and
deserved the work that went into it. The third wall is not, and did not. Only
a census that runs every increment can tell them apart, and until now nothing
did. (That run was stopped by hand once `*ERROR: increment size smaller than
minimum` had appeared and the rescue ladder was spiralling at
`dtime=1.5e-06`; it was a control, not a completion.)

**Severance is not always a reason to stop.** `test/pathfollow/close.inp` is
two blocks joined by two facets and nothing else; it drives them past `df` so
`g -> gmin` and then closes them, because the compressive branch of a crack
face is what it measures. A dead facet in compression is a real load path.
This rejected the obvious design before it was written.

## What is believed but not established

- That the remaining three discontinuities (damage initiation, loading/
  unloading, element deletion) are tractable the same way the fourth was.
  Nobody has measured a wall attributable to them.
- That the `s3rad` wall at 931 has a single identifiable cause. The only
  evidence is a residual peak on node 3053 with four live bulk elements, no
  facets, and a 0.1%-per-iteration contraction — the same *signature* as the
  second wall with a different *mechanism*. It has never been read with the
  ladder census.
