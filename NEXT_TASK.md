# Next task

## Read this first

Three things in this tree are cheap and were not, and they change how you
should work:

| | |
|---|---|
| `test/regress/run.py` | **11 cases, 47 s on four cores.** Run it before and after every change. It is proven able to go red. |
| `test/fast/run_fast.sh` | a 2160-element specimen that reaches a node with **no bulk support left**, in 24 s |
| `[SWITCHES]` at the top of every `run.log` | states the run's own configuration and names any `CCX_*` set that the binary does not read |

`s3rad` still costs ~2.5 h on 2 threads. It is the last gate, not the first.

## Settled, with evidence. Do not re-derive.

**The walls are not one thing.** Cross the two mechanisms on the fast wrapped
deck and the mask changes nothing while the kink changes everything; on
`s3rad` the opposite holds. Calling them all "the wall" is part of why fixing
one kept producing the next.

**Second `s3rad` wall (increment 554, `theta=0.2555742`).** The peak force
residual sat on node 1246, a fragment held by ONE live bulk element with all
six of its cohesive facets failed and open. Equilibrating it needs a
displacement of order one against a grip that has moved 0.2556, so its
residual fell 0.1% per iteration. The component was inherited by every
following increment - 15.7, 50.4, 74.9, 94.5, 105.7% of tolerance over
increments 551-555 - **while `dtheta` collapsed 48x**. Cutting the step made
equilibrium worse, because the part that will not reduce is not proportional
to the step. `CCX_DAMAGE_AUTOSPC_FORCE=1` extends the load-path judgement to
the force residual: increment 554 -> 930, `theta` 0.2556 -> 0.5575, `dtime`
recovered 1000x, pre-wall history bit-identical over all 1147 attempts.
Refuted for this wall, each because they repair the step and the step was
never the problem: Newton-Krylov (stops EARLIER at 0.191828), event
truncation, collapsed-node projection of the step, and a smaller increment.

**The crack-face closure kink.** `cohesive_uc6.f` leaves the normal traction
continuous at `deltal(1)=0` but jumps its slope from `g*kn` to `kn`, a factor
of `1/gmin`. On the fast wrapped deck's wall, seven UC6 points sit on that
kink and change category at EVERY line-search rung down to `eps=6.1e-5`,
while the residual grows strictly linearly in `eps` and never falls below its
base value. `CCX_UC6_CONTACT_SMOOTH=1.e-2` blends the slopes over a
penetration band: transitions 19 -> 0, `|r0|2` 3.48e-02 -> 1.89e-10, the deck
goes from `theta=0.158766` to 1.0, and the SAME 64 elements are deleted in the
same order. Flat over two decades of the band.

**The load-path judgement has one owner**, `src/damstate.c`, verified
bit-identical on a deck where the predicate actually decides and over 849
attempts / 3135 deletions of `s3rad`.

**Nothing in `s3rad` is adrift.** The diagonal test is blind by construction
to a piece that is internally stiff but attached to nothing. Measured with
`test/s3rad/fragments.py` rather than assumed: at every facet-stiffness
threshold, 0 elements adrift. The same table shows the model in **two
pieces, one at each grip**, when only facets above `g=0.5` count - the
severance result again, by connectivity instead of reaction force.

**The metal severs at `theta=0.3411981`**, and the run continues past it
because `cohesive_uc6.f` pins `g` at `gmin` while terminal deletion scans
`C3D4` only, so a dead facet reads as a load path for ever.
That is FIXED: `src/loadpath.c` owns the judgement, asks it every converged
increment, and stops the run at severance.  `CCX_FRACTURE_DEADFACET` is
retired -- it made the judgement opt-in, and the one runner that armed the
connectivity test explicitly unset it.  The fast plain deck demonstrated the
same blindness in 46 s and now stops at increment 65 instead of 507.

**Runs are not reproducible across thread counts.** `MKL_CBWR=COMPATIBLE`
fixes reproducibility across instruction sets, not across thread counts. Two
runs identical for 482 attempts diverged at 483 on a 6-against-5 iteration
count. Fix `OMP_NUM_THREADS` and `MKL_NUM_THREADS` on both arms of any A/B.

## Refuted, so you do not spend a run on them again

| | why |
|---|---|
| a "node with at most one live bulk element" gate | node 3053 at the third wall has **four** live bulk elements and no facets; the signature was the same and the mechanism was not |
| a stiffness-ratio threshold tuned until a node is included | 76 nodes below 1e-3, 166 below 1e-2, 381 below 1e-1 - that is the chain of thresholds this project exists to avoid |
| a connected-component / free-fragment unit | measured: nothing is adrift at any threshold |
| a small wrapped inclusion, to manufacture a fragment | 0 of its 24 elements erode: the wrap debonds at `Tn0=300` before ZRH yields at 600, and an unloaded inclusion cannot damage |
| a wrap stronger than the phase it wraps | the phase erodes, but the stripped nodes are then held by INTACT facets, flooring the ratio at 2.5e-2 |
| more mesh, to make the fast deck wall | four variants, none walls; the mechanism is where the facets are, not how many elements there are |

## The method that has worked

1. State a falsifiable hypothesis and **name the one measurement that can
   reject it**, before writing code.
2. Reproduce the wall on the fast deck if you can. 24 s beats 2.5 h, and
   diagnosis is what the speed buys.
3. Read `CCX_DAMAGE_WALL_THETA`'s ladder. If the residual is linear in `eps`
   and the active-set transition count does not decay, you are on a kink, not
   a stiffness loss. If the residual peaks on nodes with a healthy diagonal
   ratio and tiny `need_du = |R|/k`, likewise.
4. Regularise the kink; do not damp the step. Then prove the answer did not
   move: same deletion set, same order, bounded shift in deletion time.
5. Feature off must be bit-identical. Run `test/regress/run.py`.
6. Only then spend `s3rad`.

## In progress: an owner for "what counts as converged"

`src/convstate.c` exists and **measures only**. It reports, at the peak
residual of every iteration:

    need_du = |R| / k        and    need_du / grip

The grip travel comes from the driven node set `loadpath.c` owns, so the two
modules cannot disagree about what a grip is. It decides nothing, on purpose:
`03-TEST-LADDER.md` records that no deck reproduces the class this judgement
is for, so changing the criterion now would change a judgement no gate case
can be put against.

Measured so far:

| deck, at its wall | node | `need_du/grip` |
|---|---|---|
| `fast-wrapped`, increment 99 | 440 | **3.8e-03** |
| the class this is for | — | **of order 1** |

**Measured, and it revealed a defect in the probe rather than an answer.**
One run to the stall (4 threads, `AUTOSPC_FORCE=1`, increment 604,
`theta=0.255032` — reproduced exactly) reported:

```
[CONVSTATE] inc=604 peak residual node 8305 |R|=4.168e-02 k=6.996e+03
            need_du=5.96e-06 grip=0.2550 need_du/grip=2.34e-05
excluded:   largest excluded residual 4.4835e-02 at node 1246
```

The probe read the peak of `ram[0]` — and `ram[0]` is what remains AFTER
`AUTOSPC_FORCE` excludes node 1246. **So the one node the exclusion exists for
was the one node never measured.** Fixed: the exclusion report now also reports
the excluded node as a length. The measurement itself still has to be re-taken.

Worth keeping from that run anyway: the peak of the JUDGED residual sits on a
**healthy** node — `k = 7.0e+03`, `need_du/grip = 2.3e-05` — at 22.7x
tolerance. By `02-DIAGNOSTICS.md` §2 that is the signature of a **kink**, not
of a stiffness loss, which is a different diagnosis from the one the second
wall is on record for. One measurement, not yet corroborated.

**So the measurement that matters is still outstanding:** `need_du/grip` at
the stall, on node 1246. The prediction is that it is of order
one, and that is what would make the criterion dimensionally sound and free of
any tuned threshold. It has NOT been taken — the `s3rad` runs in this session
used a binary built before `convstate.c` existed, so their logs carry no
`[CONVSTATE]` line. Taking it costs one run to `theta=0.2556` at TWO threads
(at four the run stalls at increment 604 and never reaches the wall).

Only after that number exists should the criterion itself change, and the
change should go through `convstate.c` rather than through another norm-level
mask like `CCX_DAMAGE_AUTOSPC_FORCE` - which answers the same question with a
stiffness ratio, is a fiction by `02-DIAGNOSTICS.md` 4's own criterion at the
stall, and is not robust across thread counts.

## Open

- **The other three discontinuities** the audit lists: damage initiation
  (`deff > d0`), loading/unloading (`deff` vs `dmax`), and element deletion.
  The first two are kinks in the tangent and the method above applies; the
  third changes the mesh and is partly irreducible.
- **Viscous stabilization** in the `*STATIC, STABILIZE` sense. Note that the
  class it is normally for - a free floating fragment - was measured NOT to
  occur here, so build it against a measured failure, not against the name.
- **94 of 139 switches** have no prose anywhere but the line that reads them
  (`docs/SWITCHES.md`, generated). Retire or document.
- ~~**`src/ccx_2.22`** is a 6.5 MB executable tracked in git.~~ Removed from
  HEAD; it stays in history.

## Validation and success

Same deck, same PARDISO binary family, same environment; the only functional
difference is the change under test. Do not change physical parameters, the
deck, tolerances, solver, viscosity, tangent mode, `DEADALL` or AUTOSPC inside
a causal A/B. Validate first with crack-control engagement off.

Compare more than the stopping increment: `theta` and reaction history,
`deffmax`, process-zone and failed UC6 counts, cumulative deletions and the
exact batches, residual and iteration and cutback counts, and connectivity.

**The Newton-Krylov arm is the cautionary case: it led in `theta` at every
increment from 91 to 141 and still lost the run.** Being ahead early is not
evidence. Only the stopping state is.

## Compute budget

A full `s3rad` run to the wall is about 7400-8800 s on 2 threads, so an A/B
pair is roughly four hours. Print the budget you have chosen before launching.
A sensible default is four full runs per investigation. An extra run is
allowed when it tests a new falsifiable hypothesis - say why first. Unit and
benchmark runs do not count. Do not cap accepted increments arbitrarily and
then call reaching the cap completion.

## Delivery

Report separately: reproduced facts, new measurements, the actual cause,
functional changes, regressions, and whether `s3rad` completed or physically
separated. Passing one wall is not the same as completing the specimen. If
GitHub write access is unavailable, deliver a verified bundle with its base
and head commits.
