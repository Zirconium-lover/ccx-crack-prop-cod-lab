# 6. What the architecture should become

Not a plan to follow step by step — a set of properties to judge a change
against. Where a property already holds somewhere, that place is named, so
the pattern can be copied rather than reinvented.

## One owner per judgement

Every wall so far has been the same physical state surfacing through a
consumer nobody had wired up. The load-path judgement was computed inline and
told to consumers piecemeal: the displacement norm learned it one way, the
force norm another, the termination connectivity a third, and the solve was
never told at all. The second wall was exactly the gap between the first two.

**Pattern, already working in `src/damstate.c` and `src/lsladder.c`:** the
judgement lives in one module, the module has a self test, the mechanism
**refuses to arm if the test fails**, and every consumer asks the module
rather than recomputing. Extending this is the highest-leverage structural
work available.

Judgements that still have no owner: what counts as eroded, what counts as
converged, what counts as connected, what counts as a load path *including
facets*.

## A model that knows when it is finished

The worst defect in the tree (`05-DEBT.md` §1) is not a wrong answer, it is a
plausible one: the run continues 178 increments past physical separation
because a dead facet still reads as a load path. Connectivity should be a
first-class, cheap, **online** quantity — a census the solver prints, not a
python script somebody remembers to run afterwards. Then a run past severance
is impossible to produce without noticing.

## Discontinuities regularised, not damped

The pattern that worked: **measure the kink, regularise it, prove the answer
did not move.** Contrast it with the pattern that accumulated: add another
mechanism to survive the kink.

The proof obligations, all of which were met for the one case that is done,
and all of which should be met for the next:

- feature off is **bit-identical**;
- the far field on both sides of the transition is **unchanged**, measured on
  the law itself, not on a structure;
- the fracture does not move — the same elements deleted, in the same order,
  with a bounded shift in when;
- the result is **flat over decades** of the regularisation parameter. A
  result that is sensitive to it is a tuned constant, not a regularisation.

## Validation that is cheap, honest, and layered

The ladder in `03-TEST-LADDER.md`: a law-level check in under a second, a 3-D
specimen with real erosion in under a minute, a specimen that manufactures
the pathological state in under half a minute, and the real deck as the last
gate rather than the first. Eleven cases, 47 s on four cores, **proven able
to go red** — and the proof is on record: reinjecting the dead-facet defect
turns all eleven red and trips the `LOADPATH` self test.

Two properties matter as much as the speed:

- **cross-case relations, not just per-case numbers**: deletion sets compared
  element by element, status files compared byte for byte. "The same numbers"
  is much weaker than "the same bytes".
- **a stated purpose per case**. A case whose purpose nobody can state is a
  case nobody will fix when it goes red.

## Every run states what it is

`[SWITCHES]` at the top of every log: what is in force, and what is set that
this binary does not read. The failure it prevents is an A/B that is not
wrong but **uninformative**, and that is the expensive kind.

Generalise it: a run should be able to state its own configuration, its own
mesh and material provenance, and its own reproducibility conditions —
including the thread count, since results are not reproducible across them.

## Fewer mechanisms, each justified

Six globalization mechanisms exist and at one wall three consecutive attempts
produced bit-identical residual sequences — two of them did nothing.

The test for keeping one: **name the failure it addresses, and the case in
the gate that would go red if it were removed.** A mechanism that cannot pass
that test is not insurance, it is a place for the next bug to hide.

## Diagnostics as an interface, not as scattered probes

They currently live interleaved with the solve, gated by environment
variables, printing to stdout in ad-hoc formats. What is needed is one
interface with a stable, machine-readable output, so that `02-DIAGNOSTICS.md`
can become a program rather than a document — a run that fails should be able
to say which of the readings in that document apply to it.

## What "done" means for this work

Not that `s3rad` reaches `theta=1`. It may not, and it may be right that it
does not.

Done is: **a code base in which the result described in `05-DEBT.md` §1 —
running 178 increments past separation and calling the stop a wall — cannot
be produced without noticing.** Every other property here is in service of
that one.
