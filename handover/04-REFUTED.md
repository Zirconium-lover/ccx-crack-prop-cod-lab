# 4. Measured and rejected

Every entry is a run somebody does not have to spend again. Add to it in the
same form: what was tried, and the measurement that killed it.

## About the walls

| tried | measurement that rejected it |
|---|---|
| a Newton–Krylov corrector for the second wall | stops **earlier**, `theta=0.191828`, worse than no fix. It led in `theta` at every increment from 91 to 141 and still lost the run — **being ahead early is not evidence** |
| event truncation (stop the step at the first active-set crossing) | rejected from the eps ladder before being built: the first crossing buys 0.64% and the best rung 1.03%, against the ladder's own 0.4% |
| projecting the step onto the non-collapsed nodes | 6.7x on `|R|2`, nothing on `|R|inf` — and `|R|inf` is what is judged |
| a smaller increment | the wall's residual is **not proportional to the step**: `dtheta` collapsed 48x while the inherited residual grew from 15.7% to 105.7% of tolerance |
| global snap-back as the cause | `qa` declines smoothly to 58.3% of peak; no snap-back |
| contact chatter as the cause of the *second* wall | `ncomp=0` through increments 553–555 — no UC6 point in compression at all |
| a "node with at most one live bulk element" gate | node 3053 at the third wall has **four** live bulk elements and no facets. The signature was the same and the mechanism was not |
| a stiffness-ratio threshold tuned until a node is included | 76 nodes below 1e-3, 166 below 1e-2, 381 below 1e-1. Tuning among them is the chain of thresholds this project exists to avoid |
| the crack-face regulariser, **on `s3rad`** | moves the wall 28 increments, 0.10% of load factor — and the whole difference is inside the post-severance regime. It works on the fast deck and does nothing here |

## About the load-path judgement

| tried | measurement that rejected it |
|---|---|
| **H3: making the judgement state-dependent would make `CCX_FRACTURE_PAST_SEVERANCE` unnecessary on `close.inp`** | it does not. `close.inp` reaches `dback=1.0` with the normal traction still **positive** (`+1.008e-02` at `theta=0.84`, `+1.032e-02`, `+1.056e-02` on the next two), i.e. fully failed and OPEN, for at least three increments. The path really is gone there; the faces only meet in step 2, after the run would have stopped. The switch was **not** a workaround for a wrong judgement, and saying it was, was wrong |
| **a shut dead facet as a load path regardless of load direction** | the grip **reaction**. On `s3rad` the census read `connected=1` while the grip carried **0.044% of peak**, and the run walked to increment 931, `theta=0.5575` — the recorded phantom numbers, reproduced by a different wrong judgement. A shut facet carries compression, not tension; connectivity is directional |
| a facet's failure flag alone as "not a load path" | **withdrawn.** `cohesive_uc6.f` gives `deltal(1)<0` the FULL `kn`, not `g*kn`: a dead facet in compression carries load, and that is how a real crack transmits compression. On `fast-wrapped` exactly one such facet is the difference between a severed and a connected specimen at increment 96 |

## About the model's structure

| tried | measurement that rejected it |
|---|---|
| a free-fragment / connected-component unit (the `*STATIC, STABILIZE` class) | at **every** facet-stiffness threshold, **0 elements adrift** in `s3rad`. The class does not occur here. Build it against a measured failure, not against the name |
| viscous stabilization, by extension | same measurement. It would have been a seventh globalization mechanism beside six existing ones, for a failure this model does not have |

## About the test decks

| tried | measurement that rejected it |
|---|---|
| a full-section cohesive plane | fails first, the bulk unloads and never damages, `theta=1` with **zero deletions** |
| bigger meshes, to make the fast deck wall | four variants, none walls, none masks a node |
| moving the brittle band behind the crack tip | worst diagonal ratio `1.0405e-01`, 90 deletions — **identical** to no change |
| a small (2x1x2) wrapped inclusion | **0 of its 24 elements erode**: the wrap debonds at `Tn0=300` before the phase yields at 600 |
| a wrap **stronger** than the phase it wraps | the phase erodes, 11 nodes lose their bulk — and are then held by **intact** facets, flooring the ratio at 2.5e-2 |
| the same, stiffer wrap (`Kn=2e6`) | worse: 1.4e-01 |
| the same, weaker wrap (`Tn0=800`) | worse: 1.9e-01, only 36 deletions |
| a wrapped block spanning the **full** width | `theta=1`, worst 2.97e-01: it is the only load path, so it debonds wholesale and unloads |

## About the specimen-level load-path judgement

| tried | measurement that rejected it |
|---|---|
| "severance means stop, always" | `test/pathfollow/close.inp` is two blocks joined by **two facets and nothing else**, drives them past `df` on purpose (`xstate(4)=1.0` on both from `t=1.38`) and then **closes** them. Stopping at severance destroys the compressive step the benchmark exists to measure. A dead facet in compression is a real load path — 0.3 s to reject |
| arming the judgement inside the damage-switch block, next to everything else | that block is `if(*ndmat_>0)`. `close.inp` has no bulk damage material, so it was **never judged at all** — the same "judgement behind a gate that has nothing to do with the judgement" the module exists to remove. Arm on a property of the deck (a damage material **or** a cohesive element), not on a switch |
| re-arming per step with `loadpath_init` | `nonlingeo()` is entered once per STEP, so the severance latch was cleared at the step boundary and a specimen that came apart in step 1 was judged intact at the top of step 2 |
| trusting `--check` on the generated switch registry after regenerating it | the registry was current and the **binary** still reported the new switch as `UNKNOWN`: `damswitch.o` had no dependency on the generated header, so it was never recompiled. The `[SWITCHES]` banner caught itself. The makefiles now carry the dependency |

## Corrections to earlier conclusions, kept because the reasoning recurs

- **"The residual peaks at node 1244."** No — the *correction* peaks there.
  Both statements are true at different states (increment start vs iteration
  8) and conflating them wasted a run.
- **"The census shows N transitions relative to the base state."** No — it is
  referenced to the FULL step. Misreading the direction inverted a conclusion.
- **"DEADFACET changed the trajectory."** No — its code was read-only.
  Discriminating properly showed **thread count** was the cause. (The switch
  itself is now retired: `loadpath.c` applies the judgement unconditionally.)
- **"The fast deck validated the refactor."** No — it never masked a single
  node, so the decision branch was never exercised. The threshold is clamped
  at `1.e-1` in the source and that deck's worst node reaches `1.0406e-01`,
  missing it by 4%.
- **"The fast-wrapped deck walls on the crack-face kink."** True, and
  incomplete in the way that matters. That deck **severs at increment 96,
  `theta = 0.1587578`**, and the wall is at increment 99, `theta = 0.158766` —
  **three increments and 0.005% of load factor past severance**. The wall is
  real, the kink is real and still measured on the law itself
  (`close-sharp`, ratio 1e+06), and the regulariser still removes it — but it
  is a convergence failure of a specimen that had already come apart. It is
  kept as a regression (`fast-wrapped-wall`, byte-identical to the old run)
  rather than as a frontier.
- **"The regulariser takes the wrapped deck to `theta=1`."** It does, and
  **420 of those increments are past severance**. Compared where the specimen
  still exists, the two arms sever at `theta` 0.1587578 (sharp) against
  0.1595 (blended) — a shift of **7.4e-04**, which is the same number already
  recorded as the largest shift in any deletion time, arrived at a second way.
  That is the honest statement of what the regulariser does here.
- **"The regulariser was rejected by `s3rad`."** Also no, and this is the
  subtler one: the stopping point was compared, but the stopping point is in
  the phantom regime. Compared where the specimen still exists — severance
  `theta` — the two arms agree to **0.19%**, so the run *confirmed the change
  is harmless at full scale*. **Compare arms at a physically meaningful
  event, not at whichever increment the solver happened to give up on.**
