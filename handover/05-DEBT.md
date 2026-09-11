# 5. The debt, in priority order

## 1. A dead facet is a load path for ever — DONE, and what it cost

**Fixed in `src/loadpath.c`.** Kept here because the anatomy is the useful
part: the judgement was not missing, it was unreachable, and it was
unreachable in three independent ways at once.

`cohesive_uc6.f` pins `g = max(gmin, 1-dvisc)` and terminal deletion scans
`C3D4` only, so a fully failed facet never disappears and never stops
carrying `gmin*Kn`. The code already knew how to notice. It could not:

| gap | measured |
|---|---|
| the test ran only when `CCX_FRACTURE_TERMINATION` named two node sets by hand | off for any deck nobody configured |
| a dead facet still counted as conducting unless `CCX_FRACTURE_DEADFACET` was **also** set | `run_s3rad.sh` explicitly `unset` it |
| the call site was inside `if(damage_tent_count>0)` — the **bulk deletion transaction** | an interface-dominated fracture could never reach it however carefully it was configured |

The third is the one that hides the other two: the guard could be configured
correctly and still never run. That is why the source comment recording the
DHC1 interface case ("severance at t=0.6100 while the run was driven on to
0.7124") sits right next to a fix that could not have fired there.

### What it looked like on a 46-second deck

The fast plain deck is the same defect without the 2.3 h. 1 thread, PARDISO:

| | |
|---|---|
| grip-to-grip connectivity through live bulk lost at | increment **65**, `theta = 0.1175` |
| the run continued to | increment **507**, `theta = 1.0` |
| deletions after severance | **0** — all 90 are complete by increment 65 |
| grip reaction, `Fx/theta` over `theta` 0.19…1.0 | **constant to 0.04%** |

That last row is the reading to keep. Post-severance the model reports a
*perfect linear spring*: the reaction falls to 1.2% of peak at severance and
then climbs back to **6.5% of peak** by `theta=1`, because two detached halves
are joined by 36 failed facets at `gmin*Kn` and the grip keeps pulling. The
defect does not merely waste increments — it manufactures a recovering
load-displacement curve.

### The fix

One owner, on the `damstate.c`/`lsladder.c` pattern: a self test, a refusal to
arm if the test fails, and consumers that ask rather than recompute. It
composes rather than re-derives — the facet half is still
`damstate_facet_dead`, the walk is still `damconnect` — and changes no
equation.

- asked **once per converged increment**, unconditionally;
- endpoints from the deck's own `*BOUNDARY` cards when nobody names them: the
  driven direction against the reacting one. Reproduces the hand-named pair
  exactly on every deck here and correctly ignores the direction-2/3
  rigid-body fixpoints;
- severance **latches**, is announced once, stamps every later increment, and
  the run prints its own verdict at the end;
- on severance the run **stops**;
- the old mechanism and `CCX_FRACTURE_DEADFACET` are **deleted**, not kept
  beside it.

Proof obligations, all met and all measured rather than assumed:

- **feature off is byte-identical**: with `CCX_FRACTURE_PAST_SEVERANCE=1` both
  wrapped cases reproduce the pre-change `m.sta` and `m.damage` exactly;
- every truncated case is **byte-identical up to severance**, and every
  deletion history is byte-identical in full. What was removed is the phantom
  tail and nothing else;
- the gate goes **red** when the defect is reinjected — all 9 cases, and the
  self test trips too, so the module refuses to arm rather than misjudging.
  Two further red proofs exist since: ignoring compression (9/9) and ignoring
  load DIRECTION (3/3 on the wrapped cases).

### The one that must not be generalised

`CCX_FRACTURE_PAST_SEVERANCE` exists because **severance is not always a
reason to stop**, and the deck that proves it is already in the tree.
`test/pathfollow/close.inp` is two blocks joined by two cohesive facets and
nothing else; it drives them past `df` on purpose so `g -> gmin` (measured:
`xstate(4)=1.0` on both facets from `t=1.38`), and then **closes** them,
because the compressive branch of a crack face is what the benchmark
measures. A dead facet in compression is a real load path.

The naive design — sever means stop, always — was rejected by that deck
before any of it was written. It is the cheapest hypothesis rejection in this
document: one 0.3 s run.

### What is still open here

The **root** fix named in the original entry is still not done. A fully
failed facet is still in the equations carrying `gmin*Kn`; what changed is
that the run now *knows* and stops. `CCX_DAMAGE_FACET_DELETE` already exists
for the equation-level removal and is still default OFF and still unvalidated
— it changes the equations and carries the risk that sank `damfloatface`.
Doing it properly would let the halves become free bodies and the matrix
become singular, so the model would find out on its own rather than being
told.

## 2. 139 switches, 122 that no test sets

Generated, not estimated (`docs/SWITCHES.md`). 63 have no prose anywhere but
the line that reads them. The 122 include every rescue, corridor,
backtracking and diagnostic knob — the combinatorial space where nobody can
say what any combination does.

A switch that no test sets and no prose explains **has no defenders**.
Retiring one is now cheap to justify: make its behaviour the default, or
delete it. Prefer either to leaving it.

## 3. Six overlapping globalization mechanisms

The adaptive damage line-search ladder, transactional backtracking, Rescue
level 1, Rescue level 2 with an event step, a dogleg trust region as level 3,
and dissipation/crack-control path following. Each was added for one wall.
None was removed when a later diagnosis showed the earlier one had been
incomplete.

Evidence they do not all discriminate: at one wall, **three consecutive
attempts produced bit-identical residual sequences** — two rescue levels ran
and did nothing.

**That was not history — it was live on the configuration the target runner
uses, and it is now fixed.** The ladder already had a skip for it, but the
skip was gated on `damage_reg_nlam>0`, i.e. it only fired when the
REGULARISATION ladder was armed, and it jumped to level 2.
`test/s3rad/run_s3rad.sh` arms RESCUE2 + TR_DOGLEG and **not** RESCUE3, so
`damage_reg_nlam` is 0, the skip never fired, and the level it would have
jumped to could not act either. Measured on `s3rad`, increment 602:

| attempt | level | what it did |
|---|---|---|
| wall 1 | 1 | 0 `[DAMAGE BT]`, 0 `[DAMAGE TR]` — **nothing** |
| wall 2 | 2 | 0 `[DAMAGE BT]`, 0 `[DAMAGE TR]` — **nothing** |
| wall 3 | 3 | 12 `[DAMAGE TR]` — the dogleg **converged** |

The whole run printed `skipping straight` exactly **zero** times.

`src/rescuelevel.c` now owns "which level can act first", with a self test,
because the previous version was a condition nobody could check. Measured on
`fast-wrapped`, which has the same defect at increment 97:

```
before   97  2U   6 iter  theta=0.158758     <- did nothing
         97  3U   6 iter  theta=0.158758     <- did nothing
         97  4   11 iter  theta=0.158762     <- converged
after    97  2   11 iter  theta=0.158762     <- the same, immediately
```

Same accepted increments (98), byte-identical `m.damage`, two attempts fewer.
The gate pins it as `attempts`, so reinstating the dead levels goes red.

## 4. A 14000-line function

`nonlingeo()` carries the solve, the convergence judgement, the erosion and
topology bookkeeping, the load-path judgement, and every diagnostic, with the
diagnostics interleaved among the physics. Three units have been extracted
(`lsladder.c`, `damstate.c`, `loadpath.c`, `convstate.c`; `damswitch.c` is the
switch registry, not a judgement) and the pattern is proven: one owner, a self
test, and the mechanism refuses to arm if its test fails.

The obvious remaining units: erosion/topology, the convergence judgement, the
solve, the diagnostics.

## 5. Reproducibility is only within a fixed thread count

**This is no longer a theoretical worry about A/Bs — it has invalidated a
headline result.** The recorded `CCX_DAMAGE_AUTOSPC_FORCE=1` gain, increment
554 -> 930, was measured on **two** threads. At four, the same configuration
stalls at increment 604, `theta=0.255032`, and never passes the second wall.
Nothing in the tree stated the thread count next to the result, so it read as
a property of the switch. Quote the thread count with every trajectory number
or the number does not mean anything.


`MKL_CBWR=COMPATIBLE` buys reproducibility across instruction sets, not
across thread counts — a distinction that was being relied on silently. Two
runs identical for 482 attempts diverged at 483 on a 6-against-5 iteration
count and were on different walls twenty increments later.

Partly irreducible: quasi-brittle fracture with erosion is a sequence of
topology jumps, and a 1e-16 difference decides which side of a threshold an
integration point falls on. But nothing in the tree *states* this at run
time, and an A/B across thread counts is silently meaningless.

## 6. Three discontinuities left

The fourth — crack-face closure — was measured and regularised; the method is
in `02-DIAGNOSTICS.md` §1 and the result is that the far field on both sides
is unchanged and only the transition differs. Remaining:

| switch | jump |
|---|---|
| damage initiation, `deff > d0` | onto the softening branch |
| loading/unloading, `deff` vs `dmax` | the rank-1 softening term appears and disappears |
| element deletion, `D > 1-DEADALL` | **the mesh topology changes** |

The third is partly irreducible. The first two are not, and nobody has yet
measured a wall attributable to either — **do not regularise them on faith**.

## 7. Smaller items

- ~~`src/ccx_2.22`, a 6.5 MB ELF executable, is **tracked in git**.~~ **Done.**
  It was a build output of the inherited upstream `src/Makefile`, which
  targets `ccx_2.22.c` — a file that is not in this tree. Nothing on any build
  path used it. Removed from HEAD and added to `.gitignore`, which had only
  ever matched `ccx_2.23*`. It remains in history; removing it from there
  needs a rewrite nobody should do to a shared branch.
- Diagnostic probes are compiled into the solve path and gated by
  environment variables rather than isolated behind one interface.
- `CCX_DAMAGE_AUTOSPC` is silently clamped at `1.e-1`. Legitimate as a safety
  rail, invisible as a behaviour — a deck whose worst node sits at `1.04e-01`
  masks nobody and gives no indication why.
