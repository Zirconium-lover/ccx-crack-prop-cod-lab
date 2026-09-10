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

**The fast-wrapped wall is past severance.** That deck — the one the
crack-face regulariser was built and validated against — severs at increment
**96**, `theta = 0.1587578`; its wall is at increment 99, `theta = 0.158766`,
**three increments later**. The wall is unchanged and kept as a regression
(`fast-wrapped-wall`, byte-identical to the pre-change run), and the kink is
still real and still measured on the law itself. But it is a convergence
failure of a specimen that had already come apart, and the regularised arm's
`theta=1` is **420 increments** of phantom. Measured where the specimen still
exists, the two arms sever 7.4e-04 apart in `theta` — the same figure already
on record as the largest shift in any deletion time.

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
