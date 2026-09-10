# Task: rebuild the architecture of the CalculiX fracture branch

You are taking over a working but structurally unsound fracture extension to
CalculiX 2.23. It produces correct physics and it has passed several hard
walls, but it did so by accumulating point fixes, and the accumulation is now
the main obstacle to further progress. Your job is to fix that — the
structure, not the next wall.

**Read `handover/` before changing anything.** Six documents, in order:

| | |
|---|---|
| `handover/01-STATE.md` | what exists, what is measured, what is believed |
| `handover/02-DIAGNOSTICS.md` | **what to measure and how to read it** — the most useful document here |
| `handover/03-TEST-LADDER.md` | the test-calculation strategy, and the rung that is missing |
| `handover/04-REFUTED.md` | measured and rejected, so you do not spend a run re-deriving it |
| `handover/05-DEBT.md` | the crutches, in priority order |
| `handover/06-TARGET.md` | what the architecture should become |

## What matters, in order

**1. Remove the defects, crutches and dead ends.** This is the priority. The
code has 139 environment switches, 122 of which no test ever sets; a
14000-line function; six overlapping globalization mechanisms; and at least
one defect that silently invalidates every result past a certain point
(`handover/05-DEBT.md`, item 1 — read that first, it is the one that
manufactured a wall nobody needed to pass).

**2. Build a structure that catches whole classes of error, not instances.**
Every wall so far has been the same physical state surfacing through a
consumer nobody had wired up. Fixing them one at a time is guaranteed to
produce a new one each time — and did, three times. The pattern that works
is: give a judgement ONE owner, give that owner a self test, and make every
consumer ask it. Two units are already like this (`src/damstate.c`,
`src/lsladder.c`); most of the module is not.

**3. Make validation cheap enough to use, without making it dishonest.**
A single check must not cost 2–3 hours. But the opposite failure is real
too: a model that runs in one second reproduces neither real defects nor
real walls, and proves only that a change is harmless. The balance point has
been found empirically once — a 2160-element specimen that manufactures a
node with no load path and walls, in 24 seconds — and
`handover/03-TEST-LADDER.md` records exactly what made it work, what did not,
and which rung is still missing. Building that missing rung is a concrete,
high-value task.

**4. Make error-finding a procedure, not an inspiration.** When a run
fails, there should be a defined set of quantities to look at and a defined
reading of each. `handover/02-DIAGNOSTICS.md` is the current version of that,
assembled from what actually worked. It should become better and, where
possible, automatic.

## How to work

You have wide latitude. This is git — everything is revertible, so prefer
the bold, well-founded change over the timid one. Rewrite a subsystem if
that is what it needs. Delete a mechanism nobody can justify. Change
defaults when the evidence supports it.

What is NOT negotiable is the evidence discipline, because it is the only
reason anything here is trustworthy:

- **Before a material change, state a falsifiable hypothesis and name the one
  measurement that can reject it.** Then run that measurement and believe it.
- **A rejected hypothesis gets its functional patch reverted**, not another
  workaround stacked on top. Record the rejection — see `04-REFUTED.md` for
  the form; every entry there is a run somebody else will not have to spend.
- **Feature-off must be bit-identical**, and you must check it, not assume it.
- **Do not change physical parameters, the deck, tolerances, the solver,
  viscosity, tangent mode, `DEADALL` or AUTOSPC inside a causal A/B.** Change
  one thing.
- **Runs are not reproducible across thread counts** — measured, not
  theoretical. Fix `OMP_NUM_THREADS` and `MKL_NUM_THREADS` on both arms.
- **Report faithfully.** If a fix moves a wall by 0.1% of load factor, say
  0.1%. If a specimen already broke 178 increments before the run stopped,
  say that, and say what the run after that point is worth.

Run `test/regress/run.py` before and after every change. It is eleven cases,
under a minute on four cores, and it is proven able to go red.

## The state of the target problem, stated exactly

`test/s3rad` is the real specimen: 42807 elements, ~2.3 h on 2 threads.

- It **physically breaks at increment 753, `theta = 0.3411981`** — the metal
  loses connectivity between the grips. Verified by bisecting the deletion
  history for grip-to-grip connectivity.
- The run then **continues to increment 931, `theta = 0.5575`**, and stops
  with `rc=201`. Everything in between is two separated halves joined by 451
  cohesive facets, 322 of them at the residual stiffness floor.
- So the "third wall" is a convergence failure of a **phantom
  configuration**, not of a specimen. It is not worth passing. It is worth
  making impossible.

Do not take "drive `s3rad` further" as the goal. The goal is a code base in
which a result like the one above is impossible to produce without noticing.

## Delivery

Report separately, and do not merge them: reproduced facts, new
measurements, the actual cause, functional changes, regressions, and what is
still open. Passing one wall is not the same as completing the specimen, and
a green test suite is not the same as a sound architecture.
