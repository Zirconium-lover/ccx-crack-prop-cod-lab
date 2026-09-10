# Fast fracture regression specimen

The point of this deck is stated in the architecture audit
(`docs/ARCHITECTURE_AUDIT.md`, section 6): **validation currently costs 2.5
hours**, so no refactor of the 14000-line damage module can be checked
cheaply, and that single fact blocks every structural improvement.

| deck | elements | wall clock | exercises |
|---|---|---|---|
| `test/s3rad` | 42807 | **~2.5 h** | everything, far too slowly |
| `test/pathfollow/*.inp` | 1-D chains | seconds | cohesive law, continuation |
| **this, plain** | **2160** | **74 s** | bulk damage, deletion, cutbacks, 3-D topology |
| **this, wrapped** | **2160** | **24 s** | all of that, **plus a node with no load path left** |

The 1-D chains cannot stand in for the wall classes: a chain has no node that
can end up hanging by one tetrahedron, which is the mechanism behind the
second `s3rad` wall.  This is the smallest thing that can.

## Build and run

```sh
CCX_EXE=/path/to/ccx_2.23_pardiso test/fast/run_fast.sh <run-dir>
```

That generates the deck and runs it under exactly the environment
`test/s3rad/run_s3rad.sh` sets, then prints the load-path census.
`FAST_VARIANT=plain` selects the original specimen; the default is `wrapped`.
The deck is generated rather than committed so it cannot drift from the
generator, and its sha256 is recorded in `provenance.txt`.

## What it is

A bar of tetrahedra, split at mid-span by a **partial** cohesive plane - a
pre-crack over the first half of the section - with a brittle band ahead of
the crack tip to localise it, pulled in displacement control.

The partial plane matters and was found by measurement.  A full-section
interface makes the wrong specimen: it fails first, the bulk unloads and
never damages, and the run reaches `theta=1` with **zero deletions**.  With a
pre-crack the remaining ligament has to tear through the bulk, which is what
produces deletion and fragments.

Material and section cards are **lifted verbatim from the target deck** by
the generator rather than restated, so the physics cannot drift from what
`s3rad` uses.  (Restating them by hand was the first attempt and was wrong -
the real cards are `*Damage Initiation, Criterion=Ductile,
Evolution=Displacement, Npoints=6` and a `*Depvar 4` cohesive material.)

## Measured behaviour of the plain variant, at `4394403`

| | |
|---|---|
| return code | 0, reaches `theta=1.0` |
| increments / attempts | 507 / 544 |
| elements terminally deleted | **90** |
| `too slow convergence` events | 2 |
| wall clock, 1 thread | **58 s** on an idle box, 74 s beside a running `s3rad` |

So it fractures with real deletion and exercises the cutback machinery, in
under a minute.  That makes it usable as a characterization baseline: any
refactor must reproduce these numbers exactly before it is allowed near
`s3rad`.

## The plain variant demonstrated the severance blindness, and no longer has it

This deck used to run to increment 507, `theta=1`, having come apart at
increment 65.  It now stops itself, with no switch to remember:

```
[LOADPATH SEVERED] inc=65 time=1.175000000000e-01
                   no surviving load path between the grips
                   endpoints: derived from *BOUNDARY, direction 1: 49 reacting
                   node(s) against 49 driven
                   live elements 2106, cohesive facets 36 of which 36 fully failed
```

All 90 deletions are complete by increment 65, and every status line up to it
is byte-identical to the old run.  `CCX_FRACTURE_DEADFACET`, which used to be
what let the run say so, is retired; `CCX_FRACTURE_PAST_SEVERANCE=1`
reproduces the old walk to `theta=1` if you want it.

Without it, the identical deck walks on to `theta=1.0` at increment 507 -
**442 increments after the specimen has actually separated**.  That is the
same blindness measured on `s3rad`, where the metal severs at `theta=0.3412`
and the run continues to `0.5575`, and it is the consequence of
`cohesive_uc6.f` pinning `g` at `gmin` while terminal deletion scans `C3D4`
only: a dead facet reads as a load path for ever.

Here it costs a minute to see instead of two and a half hours.

## The wrapped variant: a node with no load path, in 24 seconds

The plain bar above was **not enough**, and the reason is exact rather than
vague.  `nonlingeo.c` clamps `CCX_DAMAGE_AUTOSPC` at `1.e-1`:

```c
if(damage_spc_g>1.e-1) damage_spc_g=1.e-1;
```

and over the plain bar's whole run the worst node reaches

```
worst=node_415_at_1.0406e-01
```

So the plain deck misses the **largest threshold the code will accept** by 4%.
It cannot mask a single node at any legal setting, and therefore cannot
exercise the judgement that AUTOSPC, `DEADFACET` and the wall fix all ask.
Measured, not assumed: `below_1e-1=0` for every increment of a completed
`theta=1` run, and zero exclusion events at `AUTOSPC=3.e-1` as well.

Node 415 is at `x=2.8, y=2.0`, on the free surface, with **9 of its 12
tetrahedra deleted and no cohesive facet at all**.  That is the whole story of
the plain deck: its crack never goes near the cohesive plane, so no node it
strips is ever held by facets, and a node with no facets either keeps bulk
support or goes to exactly zero - it never lands in between, which is where
the mask lives.

    FAST_VARIANT=wrapped test/fast/run_fast.sh <run-dir>     # the default

puts the eroding phase (`ZRH`) where the crack actually goes and wraps its
boundary in cohesive facets, which is the arrangement the target deck has -
`MATRIX` and `PLATETANGENTIAL` with `INTERFACE` between them - and the one
this generator did not reproduce.  It spans **half** the width in z, so the
specimen is not forced to tear through it and the wrap debonds progressively
instead of all at once.

| | plain | wrapped |
|---|---|---|
| elements | 2160 | 2160 |
| wall clock, 1 thread | 74 s | **24 s** |
| return code | 0, `theta=1` | 201 at `theta=0.158766` |
| elements terminally deleted | 90 | 64 |
| worst diagonal ratio | 1.0406e-01 | **4.4641e-04** |
| nodes masked at the standard `AUTOSPC=1.e-3` | **0** | **1** (node 440) |
| increments with a masked node | 0 | 14 |

Node 440 sits at `x=2.8, y=2.0, z=1.0` on the wrap surface with a bulk star of
**4 elements, all four deleted**, held by **3 cohesive facets and nothing
else**.  That is the state the target deck reaches at node 495
(`worst=node_495_at_3.6664e-07`) and it is the state the whole load-path
judgement exists to recognise.

### What this validates

The `damstate` refactor moved that judgement out of `nonlingeo.c`.  Run the
wrapped deck under the binary from `a2ee478` (the last commit before the
refactor) and under `8edfec6`:

| file | result |
|---|---|
| `m.sta`, `m.damage`, `m.cvg`, `m.de1stats`, `m.dat` | identical |
| `m.frd` | identical except the `UTIME` wall-clock stamp |
| `[DAMAGE STIFFNESS]` census, all 99 increments | identical |
| `run.log` otherwise | differs only by the 14-line self-test banner |

with the predicate **actually deciding** - 14 increments where it masks node
440 - rather than never being asked, which is all the plain deck could prove.

### What the wall turned out to be

Being able to run the wall in 24 seconds is what made it diagnosable.  The
line-search ladder at the wall reads:

```
pass 2 eps=0.003906250  |res|2/|r0|2=437.4   active-set transitions 18
pass 2 eps=0.001953125  |res|2/|r0|2=218.7   active-set transitions 21
pass 2 eps=0.000976562  |res|2/|r0|2=109.2   active-set transitions 19
pass 2 eps=0.000488281  |res|2/|r0|2= 54.5   active-set transitions 19
pass 2 eps=0.000244141  |res|2/|r0|2= 27.1   active-set transitions 19
pass 2 eps=0.000122070  |res|2/|r0|2= 13.5   active-set transitions 19
pass 2 eps=0.000061035  |res|2/|r0|2=  6.67  active-set transitions 19
```

The residual is **strictly linear in the step and never below 1**: halving the
step halves the damage it does, and the Newton direction never helps.  The
transition count does **not** decay - 12 UC6 loading/unloading plus 7 UC6
tension/compression flip at every rung, so those points sit exactly on a
switching surface and any step at all crosses it.  That is a kink, not a
stiffness loss, and the residual peaks are on nodes with `ratio` 0.75 to 1.0
and `need_du = |R|/k` of 1e-4 to 1e-6: healthy nodes needing a microscopic
displacement they never get.

`cohesive_uc6.f` is where the kink is.  The normal traction is continuous at
`deltal(1)=0` but its slope jumps from `g*kn` to `kn` - a factor of `1/gmin`,
five orders.  `CCX_UC6_CONTACT_SMOOTH=1.e-2` blends the two slopes over a
penetration band of 1% of `d0`:

| | sharp | smoothed |
|---|---|---|
| return code | 201 | **0** |
| last increment / `theta` | 99 / 0.158766 | 515 / **1.0** |
| active-set transitions at the smallest rung | **19** | **0** |
| `|r0|2` at that iterate | 3.480060e-02 | **1.885763e-10** |
| elements deleted | 64 | 64, **the same 64, in the same order** |
| largest shift in a deletion time | - | 7.4e-04 in `theta` |
| worst diagonal ratio | 4.4641e-04 | 4.4639e-04 |

It costs nothing to leave armed.  On the plain deck the two arms are
**byte-identical** - 542 attempts, 1283 Newton iterations, `m.sta` and
`m.damage` equal byte for byte, 46 s against 47 s - because the blend differs
from the sharp law only where a **damaged** facet is in compression: at `g=1`
the blended expression reduces to `kn*d` exactly.

and the result is flat over two decades of the band (`1.e-2`, `1.e-1`, `1.e0`
give 515, 514, 514 increments), which is what a regularisation should look
like rather than a tuned constant.

### Which mechanism actually does the work

Cross the two switches on this deck and the answer is unambiguous:

| `CCX_DAMAGE_AUTOSPC` | `CCX_UC6_CONTACT_SMOOTH` | rc | last increment / `theta` | deleted |
|---|---|---|---|---|
| `1.e-3` | off | 201 | 99 / 0.158766 | 64 |
| `1.e-3` | `1.e-2` | **0** | 515 / **1.0** | 64 |
| **off** | off | 201 | 99 / 0.158766 | 64 |
| **off** | `1.e-2` | **0** | 515 / **1.0** | 64 |

Masking the collapsed node changes **nothing** here - the two AUTOSPC rows are
identical to the digit.  The kink decides everything.  So this deck's wall was
never the load-path problem even though it manufactures a load-path node, and
that separation is worth keeping in mind: the walls are not one thing.  On
`s3rad` the opposite holds, where `CCX_DAMAGE_AUTOSPC_FORCE=1` alone moved the
wall from increment 554 to 930.

(The damstate A/B above is unaffected by this: it compares the census the
predicate computes, which is identical between the two binaries.  What this
table adds is that on this deck the *consumer* of that mask does not change
the outcome.)

### The wall it hits without that switch is NOT the s3rad wall

It stops with `increment size smaller than minimum` at increment 99, and
`CCX_DAMAGE_AUTOSPC_FORCE=1` does **not** move it: same increment, same
`theta`, to the digit.  The residual there is `ram[0]=1.826` against a
tolerance of `3.5e-3` - a factor of 520, not the 5.7% miss that the `s3rad`
wall was.  So this deck reproduces the *collapsed-support state* cheaply; it
does not yet reproduce the *convergence* failure that state caused in the
target deck.  Do not use it to claim a wall fix works.

### What was tried and did not work

Recording these because each one costs a minute to repeat and an hour to
rediscover.

| attempt | result |
|---|---|
| `--band-back 1`, `--band-back 2` (eroding phase straddling the crack tip) | worst `1.0405e-01`, 90 deletions - **identical** to no change |
| a small 2x1x2 inclusion wrapped in facets | **0 of its 24 elements eroded**: the wrap debonds at `Tn0=300` long before `ZRH` reaches its yield of 600, and an unloaded inclusion cannot damage |
| `--anchor` with a wrap stronger than the phase (`Tn0=1500, Gc=6`) | the phase does erode, and 11 nodes lose their bulk - but they are then held by **intact** facets, which floors the ratio at `2.5e-2`.  A node with no bulk left simply follows its facets, so nothing ever opens them |
| the same, stiffer (`Kn=2.e6`) | worse: `1.4e-01`, the stiffer wrap holds more |
| the same, weaker (`Tn0=800`) | worse: `1.9e-01`, debonds early, only 36 deletions |
| a wrapped block spanning the **full** width | `theta=1`, worst `2.97e-01`: the block is the only load path, so it debonds wholesale and unloads |

The order matters and is the point: **the facets must fail first and the bulk
erode afterwards.**  Every arrangement that erodes the phase while its wrap is
still intact produces a node that is well supported by that wrap.  Half-width
is what buys the right order - the wrap fails progressively under a stress
gradient while the surrounding matrix keeps the phase loaded.

## What the plain variant still cannot do, and why that was worth recording

The plain bar completes rather than walling.  A geometry sweep was run to try
to make it wall, and it failed:

| variant | rc | wall clock | deleted | walls? |
|---|---|---|---|---|
| `--nx 10 --ny 6 --nz 6 --seed 0.5` | 0 | 58 s | 90 | no |
| `--seed 0.3` | 0 | 54 s | 126 | no |
| `--nx 14 --ny 8 --nz 8 --seed 0.4` | 0 | 208 s | 216 | no |
| `--nx 16 --seed 0.35` | 0 | 115 s | 126 | no |

More mesh and more deletion do not produce it, and the wrapped variant above
now says why: the mechanism is not a quantity of damage, it is **where the
facets are relative to where the crack goes**.  The plain deck's crack runs
along the free surface two cells past the cohesive plane and never touches a
facet, so no amount of mesh will ever leave a node hanging on one.

## Limits

- The wrapped variant reproduces the collapsed-support **state**, not the
  `s3rad` convergence failure; see above.
- Both variants are one-thread baselines.  Runs are not reproducible across
  thread counts (measured; `test/s3rad/README.md`), so an A/B must fix
  `OMP_NUM_THREADS` and `MKL_NUM_THREADS` on both arms.
- The wrapped variant's box position `6,3,0,2,3,3` was chosen from the
  measured crack path of the plain deck, so it is tied to
  `--nx 10 --ny 6 --nz 6`.  Change the mesh and the box has to be re-placed.
