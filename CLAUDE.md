# Working brief for coding agents

**Start with [`PROMPT.md`](PROMPT.md), then [`handover/`](handover/README.md).**
This branch exists to rebuild the architecture, not to pass the next wall.

## What you are free to do

Wide latitude, deliberately. This is git and everything is revertible, so
prefer the bold, well-founded change to the timid one:

- restructure or rewrite a subsystem;
- delete a mechanism nobody can justify — the test for keeping one is *name
  the failure it addresses and the gate case that would go red without it*;
- change a default when the evidence supports it, and say so;
- add whatever diagnostics you need, including expensive ones, as long as
  they are off by default and the noisy or unsafe ones do not reach the final
  patch;
- build new test specimens, and retire ones that have stopped earning their
  runtime.

## What is not negotiable

The evidence discipline, because it is the only reason anything here is
trustworthy.

- **Before a material change, state a falsifiable hypothesis and name the one
  measurement that can reject it.** Then run it and believe it.
- **If a hypothesis is rejected, revert its functional patch.** Do not stack
  another workaround. Record the rejection in `handover/04-REFUTED.md` — every
  entry there is a run somebody else does not have to spend.
- **Feature off must be bit-identical**, and you must check it.
- **Change one thing in a causal A/B.** Not the deck, tolerances, solver,
  viscosity, tangent mode, `DEADALL`, AUTOSPC and the thing under test at
  once.
- **Pin `OMP_NUM_THREADS` and `MKL_NUM_THREADS` on both arms.** Runs are not
  reproducible across thread counts — measured, not theoretical.
- **Compare arms at a physically meaningful event, not at whichever increment
  the solver gave up on.** Breaking this rule once here inverted a conclusion.
- **Use PARDISO for any quoted `s3rad` comparison.** SPOOLES is fine for small
  unit tests.
- **Keep binaries and generated solver output out of git.**

## Before and after every change

    CCX_EXE=/path/to/ccx_2.23_pardiso test/regress/run.py -j 4

Nine cases, under a minute on four cores, exit status is the number of
failures. It is proven able to go red three separate ways, each caught first
by a self test: reinjecting the dead-facet defect (9/9 red), ignoring
compression (9/9), and ignoring load DIRECTION (3/3 on the wrapped cases). Regenerate the
switch registry with `tools/mkswitches.py` if you add or remove a `getenv`;
the gate's preflight fails if it is stale.

Preserve, unless you have measured a reason not to and said so: feature-off
behaviour, the six self tests (`DAMAGE TR`, `DAMSTATE`, `LSLADDER`,
`LOADPATH`, `CONVSTATE`, `STALLSTATE`), the analytical Mode-I and mixed-mode benchmarks, the
constitutive law checks, and the old-wall line-search A/B.

## Delivery

Report separately, and do not merge them: reproduced facts, new measurements,
the actual cause, functional changes, regressions, and what is still open.
Passing one wall is not the same as completing the specimen, and a green test
suite is not the same as a sound architecture. If GitHub write access is
unavailable, deliver a verified git bundle with its base and head commits.
