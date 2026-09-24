# Modern-population ground truth for DRIFT — study specification

**Status:** design + feasibility measurements complete (2026-09-22). No study cells run yet.
**Purpose:** establish, before any further creationist-model claim is made, whether DRIFT's
engine reproduces what population genetics expects of a *modern* population — large *N*,
modern mutation rate, long time frame — and to say precisely where it does not and why.

This is the companion to roadmap §6h. §6h asked "does the engine reproduce neutral theory at
all?" and answered yes-with-one-characterized-deviation on a small, compressed-life-history,
pool-only scenario. This study asks the harder question: **does the engine still comport when
every parameter is pushed to modern human values, and across the four axes Rob named —
mutations on/off, growing/static, normal lifespans, frontloaded/empty founder genomes?**

---

## 0. DRIFT has two genomes, and they divide the labour by design

**The separation is intentional** (Rob, 2026-09-22): the **mutation pool** exists for the case
where **selection** is invoked — it carries per-mutation fitness effects, dominance, and a DFE,
which a presence/absence bitfield cannot. The **bitfield** is the binary
presence/absence substrate, which is all that neutral variation needs and which every
haplotype-aware analysis (LD, ROH, IBD, Fst, EHH) is built on.

What matters for this study is narrower and is a matter of *current implementation*, not
architecture: **the bitfield receives no de-novo mutational input today**. Digital mutations
can be added to it easily — presence/absence is binary and applies directly to the digital
genomes — and §0a argues that doing so is the highest-value change this study could motivate.

DRIFT's two bookkeeping systems:

| | **Bitfield** (`pop.Chromosomes`) | **Mutation pool** (`pop.IndMutations` + `pop.MutationPool`) |
|---|---|---|
| Seeded by | `init_heterozygosity` / `seed_style` at `seed_year` | nothing — starts empty |
| De-novo input | **none, ever** | Poisson(`mu`) per child |
| Gated by | `track_DNA` | `track_mutations` |
| Long-run behaviour | monotone decay to loss/fixation (**today** — see §0a) | mutation-drift **equilibrium** θ = 2·Ne·`mu` |
| Read by | `CalculateSFS` (`track_SFS`), Fst, joint-SFS, het-distance, LD, ROH, IBD, `AvHet` and every `_results.csv` genetics column | `ComputeNeutralStats` (`-validate`, `track_validation`), EHH/iHS/XP-EHH, DFE/load |

**Verified, not assumed.** `pop.Chromosomes` is written only by the seeders
(`pkg/methods/seeding_*.go`) and by meiosis (`meiosisFrom`, `pkg/simulation/birth.go:462`).
No mutation path touches it. Measured directly: a 300-year run at `mu`=5 with
`track_mutations`=1 and `track_DNA`=1 reports **bitfield S = 0** while the pool holds
**S = 1630** (θ_W = 329).

> **Units — what θ is, and DRIFT's convention.** θ is the population-scaled mutation rate, the
> compound parameter neutral coalescent theory says governs how much variation a population
> holds at equilibrium. The textbook form is θ = 4·Ne·µ with µ *per gamete per generation*.
> **DRIFT's `mu` is the Poisson mean number of new mutations per DIPLOID individual per
> generation**, i.e. per two gametes, so the autosomal expectation is
> θ = 4·Ne·(`mu`/2) = **2·Ne·`mu`**, and implied Ne = θ_W/(2·`mu`) — which is what the
> validation CSV reports as `ImpliedNe_thetaW_over_2mu`. See `pkg/analysis/neutral.go:33-37`.
> θ here is also **genome-wide, not per-site**, because `mu` is a whole-genome per-individual
> count; θ_W = 669 means 669 expected segregating sites genome-wide, not a per-bp diversity.

### Three consequences that must be designed around

1. **As implemented today, the 2×2 of "mutations" × "frontloaded haplotypes" has no
   interaction term.** The two knobs act on two substrates that never meet, so a cell with
   both on is the *sum* of two independent processes. Worth confirming as a regression guard,
   but it is not a biological interaction and must not be reported as one. **This is the
   consequence §0a would remove.**

2. **`TajimasD` means two different things in two different files.** The `TajimasD` /
   `numSegSitesSFS` / `TajimaPi` / `WattersonsTheta` / `Fis` columns in `<model>_results.csv`
   are **bitfield** statistics. The D in `neutral_validation.csv` and in the
   `track_validation` report is a **pool** statistic. On a mutations-only run the former are
   all zero while the latter is fully populated; on a frontloaded-no-mutations run the
   reverse. Any table that puts them side by side without saying which substrate each came
   from is wrong.

3. **The mainstream arm must be pool-only — today.** A standard neutral equilibrium
   population is what θ describes, and only the pool reaches it as things stand. So the
   mainstream ground-truth arm runs `init_heterozygosity = 0`, `track_mutations = 1`.
   Frontloading is a *separate, creationist-model axis* — a decaying founder component — not
   an ingredient of the mainstream comparison. Mixing them produces a population that is
   neither. **Again, §0a would remove this constraint.**

---

## 0a. The highest-value change this study motivates: de-novo mutation on the bitfield

Adding a digital mutation operator to the bitfield — flip a bit at a random position in the
child's genome, Poisson(`mu_bitfield`) times per birth — is small work, and it unlocks
disproportionate value:

- **It puts the entire haplotype-aware analysis suite on an equilibrium population.** Fst,
  joint-SFS, LD decay, ROH, IBD, het-distance, `track_SFS` and the VCF export all read the
  bitfield. Today every one of them, on a long run, is reading a substrate that has drifted
  toward fixation with nothing replenishing it. With mutational input they would describe a
  real mutation-drift equilibrium and become directly comparable to published human numbers.
- **It makes the frontloading axis a genuine contrast rather than a separate universe.**
  "Frontloaded vs not" becomes a question about *initial conditions on one genome* —
  which is the question actually being asked — instead of a comparison between two
  disjoint bookkeeping systems.
- **It gives the ARGweaver line a defensible input.** See below.
- **It costs nothing when off**, and the neutral case needs no fitness effects, so none of the
  pool's per-mutation machinery has to be duplicated.

Design notes if it is built: positions must not silently collide (the bitfield is
presence/absence, so a second mutation at an occupied position is invisible — the same
identity problem `vcf_merged.go` documents), which argues for sizing `genome_bits` well above
the expected segregating-site count, and for an infinite-sites guard. The pool remains the
substrate whenever selection is on.

### And the consequence for the TMR4A / ARGweaver line

The merged VCF export (W6) is the only object carrying both substrates, and what an ARG
inference tool reads is therefore **a decaying non-equilibrium founder component plus an
equilibrium de-novo component**. For a young-single-origin model that is arguably the right
object. But it is not a standard neutral population, and any statement about what ARGweaver
recovers has to name which component the signal came from. See §4 cell X1.

---

## 1. Feasibility, measured

All timings on this machine, `go test` driving the real engine.

- **Cost is linear in N × years**, at **≈ 3–5 µs per individual-year** (measured: 2.59 µs/ind-yr
  at N=120, 2.77 at N=400, 3.49 at N=1000, 5.47 at N=2000).
- **Emergent generation time under a normal human life history is ≈ 31 years** (measured
  29.0–35.2 across a fecundity grid at lifespan 85 / maturity 20 / spacing 1–2 /
  birth_prob 1–2). This is reassuringly human and is *not* the `generation_time` parameter,
  which nothing in the engine reads — use `analysis.PedigreeGenerationTime`.
- **Equilibrium needs ~4·Ne generations.** With Ne ≈ 0.31·N (the §6h refcount baseline) and
  genT ≈ 31 y, that is ≈ **38·N simulated years**, so cost ≈ **1.5 × 10⁻⁴ · N² seconds** at µ=1.

| N | years to equilibrium | wall (µ=1) | wall (µ=70) |
|---|---|---|---|
| 500 | 19,000 | ~40 s | ~5 min |
| 1,000 | 38,000 | ~10 min | ~1.3 h |
| 2,000 | 76,000 | ~40 min | ~5 h |
| 10,000 | 380,000 | ~4.5 h | **~36 h + memory blowup** |

- **A modern mutation rate costs ~8× wall time and ~15× memory** (measured at N=500/400 y:
  µ=1 → 2.0 s / 2 MB; µ=10 → 3.2 s / 7 MB; µ=70 → 15.9 s / 58 MB). At N=10⁴ and µ=70 the
  equilibrium pool would hold ~8.7 × 10⁵ segregating sites with every individual carrying most
  of them — tens of GB. **N = 10⁴ at modern µ to equilibrium is not reachable here.**

### The lever that makes the study affordable

**Tajima's D, θ_π/θ_W, SFS shape and Ne/N are all µ-invariant** — µ scales S but not the
shape. Measured at N=500: D = −2.523 / −2.542 / −2.513 and Ne/N = 0.124 / 0.119 / 0.114 for
µ = 1 / 10 / 70. So:

> Run every **shape** test cheaply at µ = 1, and spend the modern-µ budget on a single
> **linearity** cell that confirms θ_W ∝ µ with shape unchanged.

That converts "modern mutation rate" from a cost multiplier on the whole study into one cell,
and puts N = 1,000–2,000 comfortably in reach.

---

## 2. Design principle: test differences, not absolutes

§6h established that DRIFT's neutral Tajima's D is **not** the Wright-Fisher 0 but a
characterized **−0.66**, because its monogamous, high-fecundity, birth-then-cull,
overlapping-generations reproduction has far higher offspring-number variance than WF. That
offset is real and is not going away.

So an absolute test ("is D ≈ 0?") fails for a reason that is already understood and tells us
nothing new. Every prediction in this study is therefore framed as a **differential or a
ratio** against DRIFT's own static baseline at the same N:

- growth is tested by **how much D moves** relative to static, not by where D sits;
- frontloading is tested by the **decay rate** of H, not by H₀;
- modern µ is tested by **linearity**, not by the value of θ;
- life history is tested by whether **Ne/N tracks measured offspring-number variance**.

Each of those is a real, falsifiable population-genetics prediction that survives the baseline
offset.

---

## 3. Gaps that must be closed before the study can run

Found while measuring; none is large, but the study is blocked on the first two.

- **G1 — the §6h harness cannot run the frontloading axis at all.** `RunNeutral`'s loop
  (`pkg/validation/validation.go`, `runOneReplicate`) never calls `events.Seed`, which in the
  real engine fires from `drift.go:527`. Measured: `init_heterozygosity` 0 vs 0.02 through the
  harness gives *byte-identical* pool statistics — the seeder never ran. **Fix:** either wire
  the seed step into `runOneReplicate`, or (preferred) run the frontloading cells through the
  real CLI path so the study validates the path users actually run.

- **G2 — the harness's synthetic actuarial curve kills a normal life history.**
  `buildNeutralModel` builds its own mortality curve, `risk = (0.02 + 0.010·ageGroup/5)·scale`
  — ~2 %/y at birth rising to 19 %/y at 85. Under it, lifespan 85 / maturity 20 /
  spacing 2 / birth_prob 2 goes **extinct at year 195**. The real `static/actuarial_table.csv`
  (the supplied table, header row: "WHO LIFE TABLE FOR 1999: AFR D" — a better approximation to
  world history than a modern first-world life table; **used verbatim, not to be modified**) is far
  gentler through the reproductive years (0.0055 at age 5), and the
  real path sustains the same life history fine. Measured sustainable envelope through the
  harness: `MortalityScale ≤ 0.5`. **Decision: run this study on the real actuarial table via
  the CLI, not on the harness's synthetic curve.** The synthetic curve is a fast-turnover
  convenience for §6h and should not be what "normal max lifespans" is tested against.

- **G3 — the heterozygosity-decay Ne estimator is dead code.**
  `calculateNeFromHeterozygosityDecline` (`pkg/simulation/save.go:575`) reads
  `model.HetHistory` / `model.TimeHistory`, and **nothing anywhere ever appends to either**
  (they are allocated in `initializemodel.go:117` and `manager.go:53` and never written). It
  returns −1 on every run. Two further problems if it is revived as written: its regression is
  in **years**, not generations, so it would report Ne inflated by the generation time (~31×),
  and the `ne > 50000` guard would then silently reject the inflated answer as well.
  **Not blocking:** `AvHet` is already a per-`save_interval` column in `<model>_results.csv`,
  so H(t) can be fitted post-hoc today. Fix the estimator, but do not wait for it.

- **G4 — `-validate` has no knobs.** The CLI flag is a bare bool that runs
  `DefaultNeutralConfig()` (N=120, µ=1, compressed life history). Every cell of this study
  therefore has to be driven either from Go or from a base model. Adding `-validate-n`,
  `-validate-mu`, `-validate-years`, `-validate-reps` would make the study reproducible from
  the command line; worth doing once the cells are settled.

- **G5 — the flagship modern model has mutations off.** `static/basemodels/OoA/parameters.csv`
  has `track_mutations,0` with `init_heterozygosity,0.1` — i.e. DRIFT's closest thing to a
  "modern population" is a **pure drift-decay model with no mutational input**, which by §0
  cannot reach mutation-drift equilibrium and cannot be compared to θ at all. This
  study needs a new base model (`Modern`) rather than a tweak to `OoA`.

---

## 4. The study

Two families, because §0 says the substrates do not mix. Baseline for every cell: normal human
life history (lifespan 85, maturity 20, spacing 2, birth_prob 2, real actuarial table),
`num_demes`=1, random mating, `f_neutral`=1, `linkage_model`=arm, `recombination_model`=map,
replicate seeds, sample n=40–80.

### Family P — mainstream arm, pool substrate (`init_heterozygosity`=0, `track_mutations`=1)

| cell | configuration | prediction | readout |
|---|---|---|---|
| **P1** | static, N = 1,000, µ = 1, to equilibrium | θ_W and θ_π both plateau; implied Ne = θ_W/2µ stable across seeds; F_IS ≈ 0; D at the characterized baseline, **not** −2 (a D near −2 means not yet equilibrated) | `track_validation` report; θ(t) trace |
| **P2** | P1 at µ = 70 (modern) | **θ_W scales ×70; D, θ_π/θ_W, Ne/N unchanged** | the modern-µ linearity test |
| **P3** | P1 but growing (logistic to K = 10·N₀) | D moves **further negative** than P1 by the amount expansion theory predicts; θ_π/θ_W falls (θ_W tracks the recent larger Ne, θ_π lags) | ΔD vs P1 |
| **P4** | P1 at compressed vs normal life history | Ne/N differs, and the difference tracks **measured** offspring-number variance, not census N | Ne/N + offspring-number variance |
| **P5** | P1 at N = 500 / 1,000 / 2,000 | Ne/N approximately constant; θ ∝ N | Ne/N vs N |

### Family B — creationist arm, bitfield substrate (`init_heterozygosity` > 0, `track_DNA`=1)

| cell | configuration | prediction | readout |
|---|---|---|---|
| **B1** | static, `seed_style`=population, `init_het`=0.1, no mutations | **H(t) = H₀·exp(−t_gen/2Ne)** — log H linear in generations, slope −1/2Ne | fit `AvHet` vs year ÷ genT |
| **B1h** | as B1 | **haplotype-diversity decay gives the same Ne.** With no mutational input, distinct founder haplotypes are lost at the drift rate, so haplotype diversity decays at 1/(2Ne) per generation just as heterozygosity does — an estimator that needs no mutation rate at all | haplotype counts per arm from the bitfield |
| **B2** | B1 but growing | decay **slows**; fitted Ne = harmonic mean of Ne(t) | same fit |
| **B3** | B1, long | neutral allele fixation probability = initial frequency; mean fixation time ≈ 4Ne generations | `allelesLost` / `allelesFixed` |
| **B4** | `seed_style`=max_het | H = 1 at every site but only **2 distinct founder haplotypes**; LD complete at t=0 and decaying at the recombination rate | LD decay, ROH |

> **Caveat to state in any write-up of Family B.** Neither seeder produces a realistic founder
> genome. `seeding_population` draws each site's alt frequency as ~Exp(0.005) floored at 0.005
> and then assigns alleles **independently per individual per strand**, so founder variation
> starts at very low frequency and with **essentially zero linkage disequilibrium**.
> `seeding_max_het` gives every founder strand-0 all-ones and strand-1 all-zeros, so
> heterozygosity is 1 everywhere but there are only two ancestral haplotypes and LD is perfect.
> Real founder variation is neither. The *decay rate* predictions above survive this; any claim
> about the *starting* allele-frequency spectrum or haplotype structure does not.

### Cross-family

| cell | configuration | prediction |
|---|---|---|
| **X1** | both substrates on, static, to equilibrium | **Three independent Ne estimates must agree: pool θ_W/(2·`mu`), bitfield heterozygosity-decay slope (B1), and bitfield haplotype-decay slope (B1h).** Three unrelated routes, one engine. The θ route needs `mu`; the two decay routes do not — so agreement also cross-checks the mutation rate itself. This is the single strongest test in the study. |
| **X2** | X1 vs P1 vs B1 at the same seed | pool statistics identical to P1 and bitfield statistics identical to B1 — no interaction (regression guard for §0) |

---

## 5. Suggested order

1. **G2 decision + `Modern` base model** — normal life history on the real actuarial table,
   `init_het`=0, `track_mutations`=1, `track_validation`=1. (Blocks everything.)
2. **P1 at N = 1,000** — establishes the equilibrated modern baseline and the real
   Ne/N under a normal life history (§6h's 0.313 was measured on the *compressed* one and
   should not be assumed to carry over).
3. **P2** — the modern-µ linearity cell. One run, settles the "modern mutation rate" axis.
4. **P3 / P5** — growth and N-scaling, both differential against P1.
5. **G1 + G3 fixes**, then **B1** — the frontloading axis.
6. **X1** — the two-substrate Ne agreement test. The headline result.
7. **P4**, then fold the whole thing into `-validate` knobs (**G4**) so it is re-runnable.

---

## 5a. Results

### P1 — EXECUTED 2026-09-22. **PASS.**

`drift -base-model=Modern`: N = 1,000 static, µ = 1, 40,000 y (≈ 1,290 generations ≈ 4.2·Ne),
normal human life history on the supplied actuarial table (used verbatim, unmodified), pool-only, sample n = 80.
Wall time **4 m 26 s**. Report: `results/Modern_validation_run1_year40000.csv`.

| statistic | observed | expectation | verdict |
|---|---|---|---|
| θ_W | 668.94 | — | — |
| θ_π | 584.61 | — | — |
| **θ_π/θ_W** | **0.874** | 1.0 under WF; 0.81 at the §6h compressed baseline | PASS |
| **Tajima's D** | **−0.414** | 0 under WF; −0.66 at the §6h compressed baseline | PASS |
| **SFS vs neutral 1/i** | **χ²/dof = 1.46** (dof 158) | ~1 is a good fit | PASS |
| **F_IS** | **−0.014** | ~0 (HWE) | PASS |
| Fu & Li's D | +0.240 | 0 | PASS |
| **implied Ne = θ_W/2µ** | **334** | — | — |
| **Ne/N** | **0.335** | 0.313 at the §6h compressed baseline | — |

**Three things this establishes.**

1. **Ne/N transfers.** 0.313 (N=120, compressed life history, synthetic mortality curve) →
   **0.335** (N=1,000, normal life history, real actuarial table). The emergent effective size
   is a robust property of DRIFT's reproduction, not an artefact of §6h's compressed settings.
   It is safe to plan runs with Ne ≈ 1/3 · N.

2. **A normal human life history moves DRIFT *toward* Wright-Fisher, not away.** D goes
   −0.66 → **−0.414** and θ_π/θ_W goes 0.81 → **0.874**. The §6h rare-variant excess is a
   consequence of high reproductive variance, and the compressed life history (maturity 5,
   lifespan 24, `spacing`=1 and `birth_prob`=1) has *more* of it than a realistic one. Breeding
   is never simply annual: an individual must clear a minimum gap of `spacing` years and then
   conceives with probability 1/`birth_prob` in each year after it, so `spacing`=1 +
   `birth_prob`=1 is the one setting where it does become every-year — and that is exactly
   what §6h used. `Modern` uses 2 and 2. So the
   §6h baseline is the pessimistic end of the range, and the configuration this study actually
   cares about is the better-behaved one. The remaining offset is real and still must be
   accounted for in any D-based inference, but it is roughly half what §6h reported.

3. **The SFS is a good fit to neutral 1/i** (χ²/dof = 1.46 vs the §6h compressed baseline's
   4.77). This is the strongest single piece of evidence in the run: shape, not just moments.

**Caveat — the plateau was not directly verified in this run.** `nMuts` in `_results.csv` is a
cumulative creation counter, not standing variation, so it cannot show a plateau (it rises
linearly by construction). The equilibrium evidence here is indirect but strong: a
non-equilibrated expanding pool gives D ≈ −2.5 and a badly skewed SFS (measured on this exact
model: **D = −2.858 at 600 y**), and this run gives D = −0.414 with χ²/dof = 1.46. Cells `Moderny20000` / `Moderny30000` test the plateau directly by comparing θ_W
at 20 k / 30 k / 40 k years.

---

## 6. Risks

- **R1 — mistaking a non-equilibrated run for a result.** D ≈ −2 at 2,000 y is *not* DRIFT's
  neutral baseline, it is a population only ~65 generations old. Every Family-P cell must show
  a θ plateau before it is read. This is the easiest way to publish a wrong number.
- **R2 — reading the wrong D.** See §0 consequence 2.
- **R3 — memory at large N × modern µ.** See §1. Fail fast: check heap at 10 % of the planned
  run length and extrapolate.
- **R4 — assuming §6h's constants transfer.** Ne/N = 0.313, D = −0.6611 and SFS χ²/dof = 4.77
  were measured at N=120 on the compressed life history. Re-measure under the modern
  configuration before using any of them as an expectation.
