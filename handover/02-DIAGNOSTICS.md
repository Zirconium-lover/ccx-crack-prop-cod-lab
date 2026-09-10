# 2. What to measure, and how to read it

This is the document to reach for when a run fails. Everything in it was
assembled from what actually resolved a wall; the readings are not
theoretical.

The general shape of the procedure is:

1. **Is the run still physical?** (§6, §7) — if the specimen broke 200
   increments ago, stop and fix that instead.
2. **Is it a kink or a stiffness loss?** (§1) — this single reading splits
   the two big classes and takes one armed run.
3. **Where does it live?** (§2, §3) — which node, and does it have a load
   path.
4. **Is the fix real or a fiction?** (§4) — the number that tells you a fix
   has quietly stopped working.

---

## 1. The line-search ladder census — `CCX_DAMAGE_WALL_THETA=<theta>`

**The most valuable single diagnostic in the tree.** Arms at a load factor so
the probes cost nothing for the first N hundred increments. It walks a ladder
of step fractions `eps` along the Newton direction and prints, for each:

```
pass 2 eps=0.003906250  |res|2/|r0|2=437.4   active-set transitions 18
pass 2 eps=0.000061035  |res|2/|r0|2=  6.67  active-set transitions 19
```

**Read the two columns together:**

| ratio as `eps` → 0 | transitions as `eps` → 0 | what it is |
|---|---|---|
| → 0 | → 0 | the linearisation is right; a step-size or globalization problem |
| **linear in `eps`, never below 1** | **does not decay** | **a kink** — the iterate sits on a switching surface and any step crosses it |
| erratic | → 0 | suspect the operator: assembly, a tangent that is not the differential of the residual |

A transition count that stays flat at, say, 19 down to `eps=6.1e-5` is not
"nearly converged", it is "there is no step, however small, that does not
change the model". Newton cannot converge on a kink it cannot step off.

**The sub-counts name WHICH kink**: `UC6 loading/unloading`, `UC6
initiation`, `UC6 viscous`, `UC6 failure`, `UC6 tension/compression`, `bulk
plastic`, `bulk initiation`, `bulk damage growth`. In the case that was
solved, `tension/compression` was 7 and `loading/unloading` 12 — and fixing
the first took both to zero, because the iterate stopped bouncing.

Also printed at each armed iteration: where the residual, the Newton step and
the linear-model defect *live*, split by whether a node touches a live UC6
facet, softening bulk or plastified bulk. Useful for ruling out a region.

## 2. Residual peak attribution — `Rpeak`, same switch

Per node at the residual peak: `R`, `addiag`, `addiag0`, their `ratio`,
`spc_masked`, and `need_du = |R|/k`.

| reading | meaning |
|---|---|
| `ratio` 0.75–1.0, `need_du` 1e-4…1e-6 | a **healthy** node needing a microscopic displacement it never gets. Not a stiffness problem — look for a kink |
| `ratio < 1e-3` | the node has no load path; this is the mask class |
| `need_du` of order the **grip displacement** | the node cannot be equilibrated by any physically meaningful displacement. This is the dimensionally sound criterion and it needs no tuned threshold |

`need_du` against a length scale is the criterion to prefer over any
stiffness ratio. Tuning a ratio until a particular node is included is the
chain of thresholds this project exists to avoid — at one wall there were 76
nodes below 1e-3, 166 below 1e-2 and 381 below 1e-1, and picking among them
proves nothing.

## 3. The stiffness census — printed whenever AUTOSPC is armed

```
[DAMAGE STIFFNESS] inc=818 below_1e-3=85 below_1e-2=170 below_1e-1=374
                   nonpositive=0 worst=node_495_at_3.6664e-07
```

Counts of nodes whose assembled diagonal has fallen below a fraction of
**their own intact value**, so it is dimensionless and needs no mesh-wide
median. `worst` is the extreme.

Worth knowing: `3.7e-07 ≈ gmin × (Kn·A / K_bulk)`. A node at that ratio has
lost **all** its bulk and is held by **failed** facets. A node at ~0.04 has
lost its bulk but its facets are still intact. A node at 0.1–0.3 still has
live elements. The number tells you the topology.

**`CCX_DAMAGE_AUTOSPC` is clamped at `1.e-1` in the source.** A deck whose
worst node never goes below that cannot exercise the predicate at any legal
setting — which is exactly what made the first fast deck useless for
validating it.

## 4. The exclusion report — `CCX_DAMAGE_AUTOSPC_FORCE=1`

Prints the excluded residual next to the criterion it was excluded from.

**This is the number that tells you the fix has become a fiction.** At the
wall it fixed, the excluded residual peaked at `0.0086 × qam` — 1.7x the
tolerance — on one node for eight increments, then returned to machine zero
once the fragment resolved. **If it stops returning to zero, the mechanism is
hiding a real imbalance and must be revisited.**

## 5. The convergence record — `<job>.cvg`

`RESID` and `CORR` per iteration, per attempt.

| reading | meaning |
|---|---|
| two attempts with **identical** residual sequences | the second attempt did nothing — a rescue level that is not discriminating |
| residual flat while the correction shrinks | stalling: either a kink (§1) or an irreducible component |
| residual halves when `dtime` halves | a **step-proportional** residual Newton cannot remove — cutting back will not help and may hurt |

## 6. Connectivity — now printed by the solver, `[LOADPATH]`

**This is no longer a script somebody has to remember.** Every run with a
damage material or a cohesive element prints, unprompted:

```
[LOADPATH] grip-to-grip connectivity is checked every converged increment;
           endpoints derived from *BOUNDARY, direction 1: 49 reacting node(s)
           against 49 driven.
[LOADPATH CENSUS] inc=100 time=0.1817500 connected=1 reach=9954 live=42785
                  facets=1055/5400 failed
[LOADPATH SEVERED] inc=65 time=1.175e-01
                   live elements 2106, cohesive facets 36 of which 36 fully failed
[LOADPATH] SEVERED at increment 65, theta=0.1175000, and the run ended there.
```

The census line every 50 increments answers "is the specimen still one piece"
and "how far through the interface is it" at a glance; `[LOADPATH SEVERED]`
fires once, and the closing `[LOADPATH]` line is the run's verdict on itself,
so a number quoted from a log is quoted next to its own status.

**The verdict prints however the run ends**, including at a wall. It is an
`atexit` handler, because a run that walls leaves through `stop()` — that is
`exit(201)` — and never reaches the bottom of `nonlingeo()`. The case that
most needs the verdict was briefly the one case that did not print it:
`fast-wrapped-wall` reported "never severed" while its own log carried
`[LOADPATH SEVERED] inc=96`.

The gate pins the verdict (`severed_inc`, `severed_theta`) on every case. On
a case that deliberately runs past severance that is the *only* check on
where severance was, since the stopping increment no longer says.

An element conducts iff it is not deleted and it is not a cohesive facet whose
every integration point has failed. On severance the run stops;
`CCX_FRACTURE_PAST_SEVERANCE=1` continues and stamps every later increment
`[LOADPATH PHANTOM]`.

### The offline cross-check — `test/s3rad/fragments.py <rundir>`

Still worth running for the questions the census does not answer: whether
anything is **adrift** (a component reaching no grip), and how the picture
changes as the facet-stiffness threshold is swept.

Connected components of the live bulk, with a cohesive facet counted as a
load path only above a given stiffness fraction. Answers two questions the
per-node diagonal cannot:

- **Is anything adrift?** A component reaching no grip has six rigid-body
  modes with no stiffness, and every node in it looks healthy to a diagonal
  test. This is the diagonal test's structural blind spot.
- **Is the specimen in two pieces?** Sweep the threshold and see where it
  falls apart. On `s3rad` it is two pieces, one at each grip, as soon as only
  facets above `g=0.5` count.

## 7. Severance — the run says so itself

Formerly: bisect `m.damage` in time order for loss of grip-to-grip
connectivity. That is now unnecessary — §6 prints it, and by default the run
does not go past it. Keep the bisection only to re-derive a number from an
**old** run log that predates `loadpath.c`.

**Still true, and the reason all of this exists:** A wall 178 increments
after separation is a wall in a configuration that is not a specimen, and no
amount of convergence work on it means anything.

Corroborate with the reaction from `m.dat` (`total force`): 2676.81 at
`theta=0.1955`, 27.32 at 0.2609 (**1.0%**), 1.51 at 0.4064 (**0.056%**).
Below a per cent of peak, the answer is bookkeeping.

## 8. The deletion history — `m.damage`

`element, step, increment, step_time, total_time, material, damage,
critical_ip, batch`. Cross-reference the element ids against the elsets to
see which phase erodes and where the crack actually goes. On `s3rad`: 3007 ZR
(10.8% of the matrix) and 699 ZRH (7.4% of the hydride).

Two runs' deletion **sets** compared element by element is a far stronger
statement than their counts. The gate does this.

## 9. The configuration report — `[SWITCHES]`, top of every log

Every switch in force, and every `CCX_*` name set that this binary does **not**
read. A misspelt switch produces an A/B that is not wrong but
**uninformative**, and that is the expensive kind — it costs a whole run to
notice. Check this block before believing any comparison.

## 10. The self tests — `DAMAGE TR`, `DAMSTATE`, `LSLADDER`

Run on every job that arms the relevant mechanism, and the mechanism
**refuses to arm** if its test fails. Prefer this to a test that runs
elsewhere: a judgement that is not the one that was tested should not be
allowed to decide anything.

## 11. Constitutive laws — `test/pathfollow/check_close.py`, `check_mixed.py`

Verify a law point by point against its closed form, from printed
separations, damage and tractions, with no bulk model or structural branch in
the way. A failure there is a failure of the constitutive routine and nothing
else. Tolerance `1e-5` is the `.dat` file's seven significant digits, not the
identity being approximate.

---

## What is missing from this list

- Nothing measures the **operator** directly. A tangent that is not the exact
  differential of the residual would show up in §1 as an erratic ratio, but
  there is no clean check. A directional-derivative comparison would be one.
- Nothing reports **why an attempt was abandoned** in a machine-readable
  form. `m.cvg` has to be read by eye.
- ~~The connectivity check (§6) is offline python.~~ **Done** — `loadpath.c`
  prints it every increment and stops the run at severance. What is still
  offline is the *adrift* question (a component reaching no grip), which the
  grip-to-grip walk does not answer.
- Nothing reports the load-path census in a **machine-readable** file the way
  `m.damage` and `m.sta` are. It is log text — the gate scrapes it, which
  works but is a parser against prose. A small `.loadpath` history, one line
  per census, would let two arms be compared on connectivity the way they are
  already compared on deletion sets.
