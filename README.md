# CCX crack-propagation COD workbench

**This branch is the base for the architecture rebuild.** Start with
[`PROMPT.md`](PROMPT.md), then [`handover/`](handover/README.md).
[`NEXT_TASK.md`](NEXT_TASK.md) is the previous investigation's brief and is
kept for its measurements, not as a statement of what to do next.

This repository is a public experimental workbench for continuing a
CalculiX 2.23 fracture calculation through a snap-back wall.  It is not a
production release of CalculiX.

The current code contains a bordered path-following driver, a proven
Mode-I crack-opening-control (COD) benchmark, a mixed-mode control
coordinate built on the effective separation the UC6 law itself advances
along, and a mixed-mode benchmark with a closed form.  The engineering
target is the real mixed-mode `s3rad` model included in this repository.

## Current evidence

On the clean cohesive benchmark, stock displacement control stops at the
limit point.  COD follows 798 consecutive accepted post-peak increments,
with the load factor and loaded-end displacement both turning backwards and
the reaction matching the closed-form bilinear branch.

That result proves the continuation mechanism on a small Mode-I problem.  It
does **not** prove the mixed-mode `s3rad` calculation.

On the mixed-mode benchmark - a bar cut by an INCLINED cohesive plane, so
a fixed non-trivial mode mix is loaded through the whole branch - stock
control stops after 62 increments and the mixed-mode control follows 2400
consecutive post-peak increments with the reaction on the closed form to
5e-07 and the constraint residual at 2e-18.  The same passes at 36% and
at 90% shear.

On `s3rad` the front was measured, not assumed: 63% of the process-zone
integration points carry more than half of `deff^2` in shear and 39%
carry more than 90%, so normal opening is not the coordinate there.

**The `s3rad` wall itself turned out not to be a continuation problem at
all.**  It is the adaptive damage line search: its backtracking ladder is
`{1.0, 0.5, 0.1}`, in 46.5% of its activations the best step lies below
that floor, and when nothing contracts it takes the last rung anyway - at
the wall that rung raises the residual for five iterations running while a
step of 0.03 lowers it every time.  With the ladder deepened and the best
rung taken, the same binary on the same deck with `DEADALL=1.e-2`
unchanged goes from `theta=0.212177` to `theta=0.255574`, deleting 3734
elements in 415 batches against 2416 in 234.  Three topological
explanations - orphan degrees of freedom, a detached component, and
amplification of the near-null space - were each refuted by measurement
first.  See `test/s3rad/README.md`.

**The metal severs at `theta=0.3412`, and only past the second wall.**
Replaying the committed deletion history against bulk-only connectivity puts
the last increment with a solid-material path between the grips at 752 and
the first without one at 753.  Stock died at increment 554, `theta=0.2555742`
- 0.086 of grip displacement short - so it could never have reached severance
whatever it was given.  What still joins the halves afterwards is 430
cohesive facets, and they cannot let go because `cohesive_uc6.f` pins
`g = max(gmin, 1-D)` at `gmin=1.e-5` while terminal deletion scans `C3D4`
only.  `CCX_FRACTURE_DEADFACET` was the switch that let the run say so; it
is retired.  `src/loadpath.c` now asks the question every converged increment
without being asked to, and stops the run at severance -- see
`handover/05-DEBT.md` item 1.

**The second wall is passed.**  `CCX_DAMAGE_AUTOSPC_FORCE=1`, one flag, with
the deck, `DEADALL=1.e-2`, viscosity, tangent, AUTOSPC threshold, convergence
criteria and PARDISO all unchanged, takes the run from increment 554 at
`theta=0.2555742` to increment 598 at `theta=0.2606` and still going, with
**15 new deletion batches** where stock had committed none for the last four
increments of its life.  The run reproduces the stock `.sta` attempt for
attempt over all 1147 pre-wall attempts and diverges for the first time at
the wall itself.

What the wall was: the peak force residual sat on **node 1246, a fragment
held by ONE live bulk element with all six of its cohesive facets failed and
open**.  Equilibrating it needs a displacement of order one against a grip
that has moved 0.2556, so its residual fell 0.1% per Newton iteration and was
irreducible in practice.  That component was inherited by every following
increment and accumulated - 15.7, 50.4, 74.9, 94.5, then 105.7% of the
convergence tolerance over increments 551 to 555 - **while `dtheta` collapsed
48x**.  Cutting the step 48x made the converged equilibrium 6x worse, because
the part that will not reduce is not proportional to the step.  It was also a
deadlock: the increment that could not converge was the one in which that
fragment's last element would damage and delete.  AUTOSPC already identified
such nodes and already removed them from the DISPLACEMENT norm; the fix
extends that one judgement to the force residual and nothing else.  Full
evidence, including the two honest qualifications, in `test/s3rad/README.md`.

**How the second wall was measured, before it was passed.**  It is
reproduced exactly - increment 554, `theta=0.255574`, 6780 live bulk
elements, 415 batches, 3734 elements deleted, 566 line-search activations -
and there the tangent is CONSISTENT: the linear-model defect falls like
`O(eps)` over five halvings.  What is wrong is the step.  It is
`|p_N|inf = 1.08`, a displacement of order one on a specimen whose grip has
moved 0.2556, because the front now carries 76 nodes that have lost 99.9%
of their diagonal stiffness.

**99.96% of that step sits on the three degrees of freedom of ONE of them.**
Node 1244 is held by one live bulk element and six cohesive facets that are
still in the mesh with zero stiffness; its assembled diagonal has collapsed,
its own residual is only 12% of the peak by then, and the operator inverts
that small residual into a displacement of 1.08.  (The increment BEGINS with
the imbalance on node 1244 - `|R|inf=6.78` there - and resolves most of it,
8.32 down to 0.68, in four iterations.  What it cannot do is finish: once
that node's residual is small the correction does not shrink with it, it
grows, 0.043 to 0.61 to 1.22, onto the same node.)  The remaining 29405 degrees of
freedom carry 0.04% of the step between them.  The direction is a genuine
descent direction - the transpose identity holds to twelve digits - but the
descent is spent moving a node that is barely attached, while the force
imbalance, which is at healthy nodes 1245 and 5282, hardly moves.  **The best
residual reduction available along the exact Newton direction, at ANY step
length, is 1.03%.**

**That was then measured directly, and the answer is half a fix.**  A third
pass on the linearisation ladder walks the same rungs along `p_N` with the
AUTOSPC-masked nodes zeroed.  Five controls, where the mask carries 0.2% to
88.8% of the step, show the two ladders agreeing to three or four digits - so
the probe is inert when it should be.  Where the collapsed nodes DO carry the
step, removing them converts a direction that multiplies the residual by
723.4 at full length into one whose best rung IS full length, reducing it to
0.5784.  But **at the wall itself, where the mask removes 99.9999% to 100% of
the step, the best available reduction only improves from 0.76% to 5.1%** -
a real 6.7x that still cannot converge an increment.  The nearly-detached
nodes are therefore where the step goes and part of why it is unusable, and
removing them is not by itself enough to pass the wall.

The 6069 integration points that change branch along the step - 4511 of them
UC6 loading/unloading - are a consequence of dragging that node's
neighbourhood an order-one distance, not an independent cause: none of them
crosses below `eps=0.0078`, and 5480 are in the last half of the step.  That
also refutes the obvious remedy before it was built: truncating at the first
crossing buys 0.64% against the 1.03% best and the ~0.4% the line search
already gets.  The topology is sound there, re-measured and not carried
over: one component, no floating piece, no orphan dof, no isolated equation,
and the residual orthogonal to the near-null space to twenty-one digits.

A SEPARATE first-order tangent defect was found and quantified around
increment 92, where the run first cuts back: there the defect ratio is
constant instead of `O(eps)`, the Newton iteration's linear rate equals it
to six decimals, and its size tracks the viscous factor `dtime/(eta+dtime)`,
so the missing term is damage-rate related.  It is NOT the damage-consistent
rank-1 term: that term's `dD/d(eps)` is present at max 36.4 and mean 4.18
over 6850 elements, and its action on the Newton step is 0.27% of the
operator's - two to three orders below the defect.  The missing derivative
is therefore still unnamed.  It costs iterations and cutbacks all the way
up, and it is NOT what stops the run at `theta=0.2556`: by then the
viscosity has damped it to 1.7%.

**A Newton-Krylov corrector was built for that defect, and it is REFUTED.**
Taking the Jacobian action from the residual and using the assembled tangent
only as a preconditioner does repair the increment-92 defect - it holds
`dtime` at maximum where the plain arm cuts it 16x, and it is ahead in
`theta` at every increment from 91 to 141 - but the same binary with the
flag on stops the run at `theta=0.191828`, against `0.255574` with it off,
short of even the old wall.  It helps where the defect is large and costs
where the defect is small, and the hard part of this run is where the defect
is small.  The functional patch is reverted; the A/B table is in
`test/s3rad/README.md`.  **No fix in this tree passes the second wall.**

## Repository map

- `src/` — CalculiX 2.23 sources and the experimental path-following code.
- `src/pathfollow.c` — the bordered algebra and its self test.
- `src/crackcontrol.c` — the mixed-mode control coordinate and its self test.
- `src/lsladder.c` — the backtracking ladder of the damage line search,
  extracted with a regression test that pins the defect it used to have.
- `src/topodiag.c` — connectivity, orphan-dof, rank and residual-projection
  measurements; diagnostic only, nothing on a solution path.
- `src/Makefile.ubuntu2404` — reproducible open build using SPOOLES.
- `src/Makefile.ubuntu2404.mkl`, `src/build_mkl.sh` — reproducible
  PARDISO/MKL build, out of tree so its objects can never be mixed with
  the SPOOLES ones.
- `test/pathfollow/` — analytical bar, Mode-I cohesive and mixed-mode
  benchmarks, with `check_mixed.py` and `plot_control.py`.
- `test/s3rad/` — the self-contained target deck, its baseline
  provenance, `uc6census.py` (front census from the deck's own output)
  and `mkpilotdeck.py` (reduced-output derivation with a printed diff).
- `CLAUDE.md` — working context and experimental gates for a new coding
  session.

## Quick start

On Ubuntu 24.04:

```sh
sudo apt-get update
sudo apt-get install -y gcc gfortran make libblas-dev liblapack-dev \
    libarpack2-dev libspooles-dev python3
make -C src -f Makefile.ubuntu2404 -j4
```

Reproduce the small COD benchmark before changing the implementation:

```sh
cd test/pathfollow
../../src/ccx_2.23 -i cohesive
CCX_PATHFOLLOW=1e30 CCX_PATHFOLLOW_COD=2e-6 \
    ../../src/ccx_2.23 -i cohesive
```

Then the mixed-mode gate:

```sh
../../src/ccx_2.23 -i mixed
CCX_PATHFOLLOW=1e30 CCX_CRACK_CONTROL=4e-6 CCX_CRACK_CONTROL_ENGAGE=1 \
    CCX_PATHFOLLOW_DTHETA=2e-4 ../../src/ccx_2.23 -i mixed
python3 check_mixed.py .
```

See `test/pathfollow/README.md` for the complete evidence and diagnostics.

## Important solver distinction

The historical `s3rad` wall was measured with PARDISO and the target deck
contains `*Static, Solver=Pardiso`.  The open Ubuntu makefile builds SPOOLES,
which is useful for unit tests and exploratory comparisons but is not a
silent substitute for the recorded PARDISO baseline.  A coding session must
either establish a reproducible PARDISO/MKL build or clearly label a copied
SPOOLES deck as exploratory.  That build now exists:

```sh
sudo apt-get install -y intel-mkl        # universe/multiverse, 2020.4.304-4
./src/build_mkl.sh                       # -> build-mkl/ccx_2.23_pardiso
```

Every `s3rad` number in this repository was produced with it.  No SPOOLES
run is quoted against the PARDISO baseline.

Generated binaries, solver outputs, logs and plots must remain untracked.
