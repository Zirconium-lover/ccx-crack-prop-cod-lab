# 3. The test ladder, and the rung that is missing

The constraint that shaped everything else: **validation used to cost 2.5
hours**, so nothing could be checked cheaply and every fix was made blind.
That is the mechanism behind "we fix one wall and another appears".

The opposite failure is just as real and less often noticed: a model that
runs in one second reproduces neither real defects nor real walls. It can
prove a change is *harmless*. It cannot prove one *works*.

## The ladder as it stands

| rung | cost | what it can prove | what it cannot |
|---|---|---|---|
| `test/pathfollow/close.inp` | **0.7 s** | a constitutive law, point by point, against its closed form | anything about a structure |
| `test/pathfollow/mixed.inp` | 23 s | the mixed-mode post-peak branch, snap-back, continuation | anything 3-D; its facet never closes |
| `test/fast` **plain** | 46 s | bulk damage, terminal deletion, cutbacks, 3-D topology, 90 deletions, and severance at increment 65 | it never masks a node and never walls |
| `test/fast` **wrapped** | **24 s** | **a node with no bulk support, held by facets alone**, severance at increment 96, and a wall three increments later | the `s3rad` convergence class |
| `test/s3rad` | **2.3 h** | everything | cheaply |

`test/regress/run.py` runs the first four as **eleven cases in 47 s** on four
cores, and is proven able to go red: reinjecting the dead-facet defect turns
all eleven red *and* trips the `LOADPATH` self test, so the module refuses to
arm rather than misjudging quietly.

## What made the 24-second rung work — and what did not

This is the transferable part, because the first four attempts failed.

**Size is not the lever.** Four bigger meshes were tried; none walls, none
masks a node. More elements and more deletion do not produce the mechanism.

**Where the facets are relative to where the crack goes IS the lever.** The
plain deck's crack runs along a free surface two cells away from the cohesive
plane; its worst node has 9 of 12 tetrahedra deleted and **no facet at all**,
so it either keeps bulk support or goes to exactly zero — it never lands in
between, which is where the interesting state is.

**The ordering is the lever.** For a node to end up with no bulk support and
failed facets, **the facets must fail first and the bulk erode afterwards.**
Every arrangement that erodes the phase while its wrap is intact leaves the
stripped node well supported by that wrap — the diagonal ratio floors at
2.5e-2 and never reaches 1e-3.

**Half-width wrapping is what bought the right order.** Wrapping the eroding
phase across the full section makes it the only load path, so it debonds
wholesale and unloads. Wrapping half the width lets the wrap fail
progressively under a stress gradient while the surrounding matrix keeps the
phase loaded.

Failures worth not repeating, each cheap to re-check and expensive to
rediscover: a small wrapped inclusion erodes **0 of its 24 elements** (the
wrap debonds at `Tn0=300` before the phase yields at 600, and an unloaded
inclusion cannot damage); a wrap made *stronger* than the phase erodes it but
then holds the stripped nodes on **intact** facets; moving the brittle band
behind the crack tip changes nothing at all.

## What stopping at severance did to the ladder's cost

Same machine, same runner, `-j4`, one solver thread per case, before against
after — the only difference is stopping at severance:

| case | before | after |
|---|---|---|
| `fast-plain` | 29.5 s | **9.4 s** |
| `fast-wrapped` | 15.8 s | **13.7 s** |
| `fast-wrapped-smooth` | 33.7 s | **12.3 s** |

(wall clock under 4-way concurrency, so it jitters by a few tenths; the
`fast-plain` ratio is the one that is not noise.)

Nothing was lost to buy it: every deletion in both decks is complete at or
before the severance increment, and every status line up to it is
byte-identical to the pre-change run. The time that disappeared was being
spent solving two separated halves joined by failed facets.

The whole suite went from nine cases in 45.5 s to **eleven cases in 47 s**,
so two new cases were added for free.

This is the rare case where correctness and speed are the same change, and it
scales: a fifth of `s3rad`'s runtime was the same regime.

## The rung that is missing

**Nothing reproduces the `s3rad` convergence class.** The wrapped deck
manufactures the *state* — a node with no load path — but the wall it hits is
a different mechanism: `CCX_DAMAGE_AUTOSPC_FORCE=1` does not move it, and its
residual misses tolerance by 520x where the `s3rad` wall missed by 5.7%.

The class to reproduce is: **a node with a healthy diagonal whose residual is
irreducible** — `need_du = |R|/k` of order the grip displacement, so no
physically meaningful step can equilibrate it, while the diagonal test sees
nothing wrong.

**That statement is now a number, not a claim.** `convstate.c` reports
`need_du` and `need_du/grip` at the peak residual of every iteration, so the
distance between the two wall classes can be read off a 24-second run:

| deck, at its wall | node | `k` | `need_du/grip` | which class |
|---|---|---|---|---|
| `fast-wrapped`, increment 99 | 440 | 1.63e+01 | 3.8e-03 | neither, cleanly |
| `s3rad` stall, increment 604 | **1246** | **1.06e+00** | **1.7e-01** | no load path |
| `s3rad` stall, increment 604 | **8305** | **7.00e+03** | **2.3e-05** | **THIS one** |

**Read the `k` column before the ratio, because the two `s3rad` rows are two
different patologies and an earlier version of this file conflated them.**

- Node **1246** has a diagonal that has COLLAPSED — 1.06 against a healthy
  6996 — so its large `need_du` follows from a small `k`. That is "this node
  has no load path", and `AUTOSPC` and `loadpath.c` already own it. It is
  **not** the missing rung, because the diagonal test sees it perfectly well.
- Node **8305** is the one this section is about: diagonal **healthy**, and
  its residual nevertheless sits at **22.7x tolerance** and will not come
  down. `need_du/grip` is 2.3e-05 — *small*, because `k` is large. By
  `02-DIAGNOSTICS.md` §2 that reading is a **kink**.

So `need_du/grip` is a good detector of the first class and says almost
nothing about the second. A candidate rung for THIS class has to reproduce:
a healthy diagonal, a residual stuck at tens of times tolerance, and a
`need_du/grip` that stays small while it happens. The screening number is
therefore not "`need_du/grip` of order 1e-1" — it is the pair
(`k` healthy, `|R|/tol` ≫ 1), with the eps ladder of §1 as the confirming
reading.

`fast-wrapped` reaches neither: its residual misses tolerance by 520x, which
is a different failure again.

Building this rung is the highest-value test work available. It would make
the one remaining wall class diagnosable in seconds instead of hours, and it
is the last thing standing between the gate and being a real substitute for
`s3rad`.

## Rules for a new rung

- **Generate the deck, do not commit it**, so it cannot drift from its
  generator. Record the sha256 in the run's provenance.
- **Lift material and section cards verbatim from the target deck.** Restating
  them by hand was tried and was wrong.
- **State what the rung is FOR** in one sentence, in `cases.json`. A case
  whose purpose nobody can state is a case nobody will fix when it goes red.
- **Prove it can fail.** A suite that cannot go red is decoration. Force a
  behaviour change and check the right cases turn red with the right numbers.
- **Pin the trajectory, not just the outcome**: increments, deletion counts,
  and where possible the deletion **set** compared element by element and the
  status file compared byte for byte.
- Fix the thread count. Runs are not reproducible across thread counts.

## Where the time actually goes

`s3rad` at 42807 elements is 2.3 h on 2 threads. It should stay the last
gate, not the first. But note that a large part of its cost is spent **after
the specimen has broken** — increments 753 to 931 are a phantom regime. Fix
that (`05-DEBT.md`, item 1) and the target deck itself gets substantially
cheaper, which is the rare case where correctness and speed are the same
change.
