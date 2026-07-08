# DRIFT — Roadmap & TODO

Feature backlog for evolving DRIFT from a working base model into robust population-modeling
software. Items are grouped by theme and roughly ordered by research payoff. Stubs are fine —
the goal is to capture direction so work can proceed incrementally.

Status key: `[ ]` not started · `[~]` in progress · `[x]` done · `[?]` needs design decision

---

## 0. Quick wins / foundations (do first)

> **Phase 0 landed 2026-06-25.** Runs are now reproducible under a fixed seed, and genetics is
> active on the user-model (GUI) path. Verified: two runs with the same `rng_seed` produce
> byte-identical results CSVs; genetics-active baseline (n=100, t=500, no animations) ≈ 0.48 s.

- [x] **Deterministic, per-run RNG.** Replaced all 8 `rand.Seed(time.Now().UnixNano())` reseeds and
  every direct `math/rand` call with a single shared, seeded generator (`pkg/utils/rng.go`:
  `SeedRNG`, `RandIntn/Float64/Shuffle/Poisson`). Seeded once per run in `drift.go` from a new
  `rng_seed` param (auto-generated + logged when unset; each run gets `baseSeed + run`). The gonum
  `distuv.Poisson` draw was replaced with a seed-driven `RandPoisson` (gonum dep removed via
  `go mod tidy`). Unit tests in `pkg/utils/rng_test.go`.
  - [x] **Map-iteration determinism.** A fixed seed alone was *not* enough — Go randomizes map
    iteration order, and the sim consumed RNG while ranging over `pop.IndData` / `ChromosomeArms` /
    `CumulativeProb`. Sorted the ID/chromosome/age lists before RNG consumption in `Mating`,
    `Birth` (incl. `createMask` meiosis), `Death` (`generateKeyList`), `SetupPopDefault`, and the
    single/population seeders. This was the actual blocker to byte-identical runs.
- [x] **Fix `genome_bits` on the user-model path.** Root cause was deeper than a missing call:
  `LoadChromosomesFromPath` (user-model path) parsed a **5-column-per-chromosome** format while the
  actual `chromosome_data.csv` files are **4-column-per-arm**, so it silently skipped every row →
  empty `ChromosomeArms` → `genome_bits` = 0 → genetics a no-op. Fixed by making both loaders share
  one parser (`parseChromosomeRecords` in `pkg/utils/csv.go`) plus a `SetGenomeBits` helper called
  by every load path. Also initialized `FreeParameters["seed"] = -1` (+ indID/mutID/last_pop_size)
  in `manager.go` `LoadModel`, without which the founder genome was never even seeded on this path.
- [ ] **Fix or retire the `-base-model` CLI path.** `InitializeModelWithBaseModel` loads base params
  then calls `utils.LoadParameters` (`pkg/utils/csv.go`), which *resets* `model.Parameters` to
  `parameter_defaults.csv`, discarding the base model's overrides. Symptoms: wrong model name
  (writes `Default_results.csv`), wrong flags (`track_map` re-enabled → empty-map load → panic).
  Only the user-model path (`-username/-model`) is currently usable. Found 2026-07-02.
- [x] **Activate mutation dominance.** `Mutation.Dominance` was hardcoded to `0` and ignored, and
  the old `CountFitnessAndMutations` was **broken** — it ranged over `pop.IndMutations[child]` whose
  keys are strand indices (0/1), not mutation ids, so it summed the Effect of mutation ids 0/1 only
  and reported `numMutations` as the strand count (≤2). Rewrote it as **per-locus genotype fitness
  with dominance** (`pkg/utils/mutation.go`): each distinct mutation id is one biallelic locus
  (infinite-sites); carried on both inherited strands ⇒ homozygous ⇒ full `Effect`; on one strand ⇒
  heterozygous ⇒ `h·Effect`. Load sums additively (fitness = 1 + load). `Mutation.Dominance` is now
  populated at creation from a new global **`fitness_dominance`** param (h·100; default **0.5** =
  additive, 0 = recessive, 1 = dominant) and read back per-locus, so h travels with each mutation
  (a future DFE-linked per-mutation h is a one-line change). Homozygosity = the *same* mutation id
  on both strands (identity by descent) — the mechanism that makes recessive load express under
  bottlenecks/founder events. **Compatibility:** no new RNG is drawn and neutral runs (Effect=0)
  give load=0 ⇒ byte-identical (§6h validation unaffected); only non-neutral mutation runs change
  (and they were broken before). Unit tests in `pkg/utils/mutation_test.go`; verified e2e
  (recessive hides het load, dominant depresses fitness). Landed 2026-07-07.

## 1. Genetics realism

- [x] **Per-locus genotype fitness with dominance** (see §0) instead of strand-summed effects.
  Done 2026-07-07 alongside the §0 dominance item: `CountFitnessAndMutations` now scores each
  mutation locus by zygosity (hom = full `Effect`, het = `h·Effect`) rather than the old
  strand-keyed sum. Epistasis / multiplicative-vs-additive across loci remains the next §1 item.
- [x] **Epistasis / synergistic load.** A new **`fitness_model`** string param
  (`additive` (default) | `multiplicative` | `synergistic`) toggles how the per-locus contributions
  from `CountFitnessAndMutations` combine into relative fitness (`pkg/utils/fitness.go`:
  `CombineFitness`/`FitnessModel`). With c_i = locus i's contribution (full `Effect` if homozygous
  derived, h·`Effect` if het; from the §0 dominance scoring) and L = Σc_i the additive load:
  **additive** w = 1 + L; **multiplicative** w = Π(1 + c_i) (independent loci — can't cross zero, so
  always ≥ additive for deleterious load); **synergistic** w = 1 + L + β·L·|L| (deleterious load
  *accelerates* — each further mutation costs more than the last, the genetic-entropy claim), where
  β is a new **`epistasis_coefficient`** param (default 1; β = 0 reduces synergistic exactly to
  additive; for uniform effect s this is the standard quadratic w = 1 + n·s − β·n²·s²).
  `CountFitnessAndMutations` now takes the model and returns the *combined fitness* (was the additive
  load); `birth.go` uses it directly. **Compatibility:** the additive path is unchanged arithmetic
  (fitness = 1 + Σc_i, same sorted-id summation order, no new RNG); with no non-neutral mutations
  every c_i = 0 ⇒ all three models return exactly 1.0 ⇒ neutral runs byte-identical (§6h `-validate`
  still PASS). Both new params added to `parameter_defaults.csv` (Mutation group). Unit tests
  (`mutation_test.go`: per-model combination values, additive-unchanged, and an explicit
  zero-load-⇒-1.0 no-op guard for all three). **Verified end-to-end** (three models identical but
  for `fitness_model`, `selection_mode=none` so genotypes stay byte-identical — nMuts=3190 in all):
  mean fitness ordered synergistic 0.97372 < additive 0.97453 < multiplicative 0.97488, the expected
  severity ordering. Landed 2026-07-08.
- [x] **Mutation-class spectrum.** `GenerateNewMutations` (`pkg/utils/mutation.go`) now draws a
  configurable spectrum of mutation classes (point / indel / CNV / large deletion, …) instead of a
  single uniform `mu` + one Weibull DFE. Each class (`pkg/utils/mutation_class.go`: `MutationClass`)
  carries its own **rate** (Poisson mean), effect distribution (**reuses the existing Weibull
  machinery** — per-class `shape`/`scale`/`Weibull_adj` + `f_neutral`/`f_beneficial`), dominance,
  and a **target `size`** (bits). The kernel iterates classes in spec order, drawing
  `Poisson(class.rate)` per class and stamping each mutation with its **class index** (new
  `core.Mutation.Class`) and **size** (`core.Mutation.Size`) so downstream stats can filter/partition
  by class. Configured via a new **`mutation_classes`** string param — a `;`-separated list of
  `<name> key=value …` entries (fields whitespace/`;`-separated, never commas ⇒ no CSV quoting;
  omitted fields inherit the model globals, e.g. `point rate=8; indel rate=1 scale=0.15 size=10`).
  **Compatibility:** with `mutation_classes` unset the spectrum collapses to a single "point" class
  whose rate/DFE are the model globals, so the kernel draws the **identical RNG stream** (one
  `Poisson(mu)`, then the same per-mutation draws) and builds the identical pool as before — a strict
  no-op (`TestDefaultSpectrumByteIdentical` asserts field-for-field pool equality + aligned RNG tail
  vs an inlined copy of the legacy kernel; `drift -validate` still **PASS**, D=−0.89, π/W=0.740,
  Ne/N=0.188, unchanged). Param added to `parameter_defaults.csv` (Mutation group, empty default).
  Unit tests (`mutation_class_test.go`: byte-identity, spec parsing + field inheritance/overrides +
  malformed-token skip, and a real-kernel run over 3000 individuals confirming counts track rates
  and mean |effect| tracks scale). **Verified end-to-end** through the real engine (temp
  `smoke/MutClassTest`, 3 classes point/indel/large_del, checkpoint→reload→tally): at year 400 the
  segregating pool partitions by class with counts 91:18:10 (rates 0.7:0.2:0.1), sizes 1/10/1000, and
  mean |effect| ordered 9e-9 < 1.5e-7 < 3.9e-7 (scales 0.01/0.15/0.5) — distinguishable rates and
  effects. Landed 2026-07-08. NOTE: `core.Mutation` gained two fields; gob tolerates the addition so
  old checkpoints still decode (Class=0/point, Size=0), no `SchemaVersion` bump.
- [ ] **Variable mutation rate.** Regional rate variation / hotspots instead of one uniform `mu`
  with uniformly-random position (`pkg/utils/mutation.go:48-54`).
- [ ] **Linkage-aware mutation positions.** Tie mutation position to the `ChromosomeArms`
  structure so LD, recombination, and selection interact properly.

## 2. Selection & demography

- [x] **Resolve fitness/survivorship TODO.** Fitness now acts on birth, survival, both, or neither,
  selected by a new **`selection_mode`** string param (`fecundity` (default) | `viability` | `both`
  | `none`), plumbed via `pkg/simulation/selection.go`. `birth.go`'s fecundity gate is now
  `fecunditySelection(model)`. **Also fixed an inverted viability bug:** `deathStandard` multiplied
  death risk by raw `fitness` (≤1 for deleterious load), which made *less-fit* individuals die
  *less*; it now multiplies by `viabilityHazardFactor(fitness) = clamp₀(2 − fitness)` (= 1 − load),
  so lower fitness ⇒ higher hazard, gated on `viabilitySelection(model)`. The `selection` dropdown
  param and `pkg/methods/death_fitness.go` (empty stub) were both dead and are left as-is.
  **Compatibility:** neutral runs have fitness = 1 ⇒ hazard factor 1 and birth prob 1 under every
  mode ⇒ byte-identical (§6h unaffected); only non-neutral runs change (and the old viability
  direction was a bug). Unit tests in `pkg/simulation/selection_test.go`; verified e2e (modes
  diverge; viability raises mean fitness 0.99411→0.99438 by purging load). Landed 2026-07-07.
- [ ] **Density dependence / carrying capacity.** Add a logistic carrying-capacity term (ideally
  local to map density) so runs aren't pure exponential-growth-or-extinction. Enables realistic
  post-bottleneck recovery curves.
- [ ] **Scriptable environmental events.** Generalize the `events` package (currently just `Seed`)
  into a registry of scheduled events: famine, epidemic, migration pulse, etc.

## 3. Spatial & ecological depth

- [ ] **Habitat suitability maps.** Let map cells carry carrying-capacity / mortality multipliers
  (builds on existing `Lat`/`Lon` + `Wander`).
- [ ] **Barriers and corridors** affecting movement.
- [ ] **Deme / island models** with migration matrices — enables real population structure and Fst.

## 4. Analysis, validation & reproducibility

- [ ] **Deterministic seeding** (see 0 — also a prerequisite for everything in this section).
- [ ] **Fst / population-structure statistics** to complement existing SFS / LD / coalescence.
- [ ] **Runs-of-homozygosity** analysis.
- [ ] **Fixation / mutation-load time series** output.
- [ ] **Parameter sweeps / batch experiments.** Expand a parameter grid into many queued runs
  (the webserver already has a `queue`), with aggregated cross-run summary output. Turns DRIFT
  from "one sim" into an experiment platform.

## 5. Software hardening

- [ ] **Unified, validated config schema** with clear errors. Currently three loading paths
  (base-model / user-model / legacy) in `drift.go:66-92`.
- [ ] **Checkpoint / resume** for long runs (serialize `Pop` to disk).
- [ ] **Test suite** for genetics kernels (meiosis, mutation inheritance, SFS) so refactors are safe.
- [ ] **Structured run manifests.** Write the full resolved parameter set + seed + git SHA
  alongside every result for traceability.

## 6. Population-genetics statistics & out-of-Africa testing

The strategic goal: make DRIFT output directly comparable to (a) real human data and (b) the
published OoA demographic models, expressed in the field's own statistics, reproducibly. A result
is only defensible against critics if it's stated in a stat that can't be dismissed as
non-standard. Items below are ordered along the recommended sequencing thread.

Serves three project goals — **(1)** creationist models of human history, **(2)** testing the
out-of-Africa story, **(3)** countering critics.

### 6a. Real-data interoperability (credibility multiplier) — goals 2, 3

- [~] **VCF import/export.** Lets the *same* pipeline (PLINK, vcftools, ADMIXTURE, EIGENSOFT) run
  on simulated and on real 1000 Genomes / HGDP data — critics can't dismiss results produced with
  their own tools. **Highest-leverage single feature.**
  - [x] **Export** (`pkg/analysis/vcf.go`: `ExportVCF`), wired into `drift.go` behind the
    `export_VCF` flag (default off; requires `track_DNA`). Emits VCF v4.2: each genome bit becomes
    one biallelic SNP with a fixed REF=A (ancestral, bit 0) / ALT=T (derived, bit 1) convention;
    genotypes are **phased** (`|`) since the two strand copies are tracked through meiosis; POS is
    per-chromosome 1-based (base = min arm start; contig length = arm span) derived from
    `ChromosomeArms`; INFO carries AC/AN. `vcf_sample_size` bounds columns (0 = all); segregating
    sites only unless `vcf_include_fixed=1`. Reuses `SampleIDs`. Unit tests in `vcf_test.go`;
    verified end-to-end (QuickTest, seed_style=1/`init_heterozygosity`>0 so standing variation
    persists) → structurally valid VCF (field counts, POS-in-contig, GT form all check out).
    NOTE: population seeding with `init_heterozygosity=0` (or a lone single-founder seed that drifts
    to loss) yields empty `pop.Chromosomes` and thus an empty VCF — genome tracking needs standing
    variation to export.
  - [ ] **Import** — parse real 1000G/HGDP VCF into DRIFT structures (larger, separate task).
- [ ] **Read real recombination maps and real human chromosome structure** (builds on the existing
  `ChromosomeArms` scaffold) so LD patterns are comparable to real data.
- [ ] **Tree-sequence (tskit) output** (NEAR-TERM — promoted from longer-term; see §7). Becoming the
  field lingua franca and a massive scaling/compression win (millions of individuals × whole
  genomes become tractable). Opens the `tskit`/`msprime`/`Relate`/`tsinfer` ecosystem. **Three of
  the five §7 cutting-edge bets depend on this** — it is the enabling substrate, not a nicety.

### 6b. Differentiation & admixture statistics (the modern OoA toolkit) — goal 2

- [x] **Fst** between populations — Hudson's and Weir & Cockerham estimators. (See also §4.) The
  OoA story is fundamentally about continental differentiation. **Landed 2026-07-06.**
  - **Deme infrastructure (island model).** DRIFT was single-population, so §6b needed a definition
    of "populations." Chosen (over on-the-fly Lat/Lon binning) to add a first-class `Deme` field to
    the `IndData` enum (`pkg/core/types.go` + `pkg/individual/individual.go`, kept in lockstep).
    Founders are partitioned across `num_demes` demes round-robin, and — when a real map exists —
    also placed in distinct spatial blocks (`partitionLandByDeme` in `setup_pop_default.go`) so
    demes are geographically separated too. Deme is inherited maternally at birth (`createChild`).
    `Mating` now, for `num_demes > 1`, first migrates individuals between demes with prob
    `deme_migration_rate` (the island-model gene-flow channel), then dispatches the selected mating
    module **once per deme** — reusing every existing mating module unchanged. `num_demes = 1`
    (default) is a strict no-op: byte-identical to pre-deme behavior (verified). Checkpoint
    `SchemaVersion` bumped 1→2 (the new field lengthens the serialized `IndData` slice; old saves
    rejected loudly).
  - **Estimators** (`pkg/analysis/fst.go`: `ComputeFst`/`SaveFst`, behind the `track_Fst` flag).
    Both computed as a ratio-of-sums across sites (Bhatia et al. 2013): **Hudson's Fst** (pairwise)
    and **Weir & Cockerham 1984 θ** (pairwise + a global value over all demes, using observed
    heterozygosity from the phased strand copies). Reuses `SampleIDs`/`buildContigs`/`bitAt`; a
    shared `demeMembers`/`countSite` substrate that f-stats + joint-SFS will build on. Demes with
    <2 sampled diploids are dropped; <2 eligible demes → empty result. Writes
    `<model>_fst_run%d_year%d.csv` (one row per deme pair + a global-θ row). Params: `num_demes`,
    `deme_migration_rate`, `track_Fst`, `fst_sample_size` in `parameter_defaults.csv`.
  - **Tests** (`fst_test.go`): hand-computed fixtures — complete differentiation → Fst = 1;
    all-heterozygous identical demes → W&C = 0 exactly, Hudson ≤ 0; intermediate → Hudson = 1/3,
    W&C = 1/2; plus the guard cases. **Verified end-to-end** (QuickTest, seed_style=1,
    init_heterozygosity=0.1, 3 demes): isolated demes (m=0) → global θ ≈ 0.37; high migration
    (m=0.5) → θ ≈ 0.0004 with all 3 demes persisting; `num_demes=1` → no Fst; fixed `rng_seed` →
    byte-identical results (determinism intact).
  - NOTE: with no migration and a shared carrying capacity, a small deme can drift to extinction —
    expected island dynamics, but tune `num_demes`/`start_pop_size`/`max_pop_size` if all demes must
    persist to the end.
- [x] **f-statistics: f2, f3, f4, and D-statistics (ABBA-BABA).** The Patterson/Reich tools used to
  argue admixture, tree topology, and "Africa as outgroup." Non-negotiable for testing/countering
  OoA — it's how the claims are actually made. **Landed 2026-07-06.**
  - **Estimators** (`pkg/analysis/fstats.go`: `ComputeFStats`/`SaveFStats`, behind `track_fstats`).
    All sample-size-unbiased and computed as a mean/ratio of sums across polymorphic sites
    (Patterson 2012 / Peter 2016), reusing the fst.go substrate (`demeMembers`/`countSite`) and
    `buildContigs`/`bitAt`. Per site: **f2**(A,B) = (pA−pB)²−pA(1−pA)/(nA−1)−pB(1−pB)/(nB−1);
    **f3**(T;B,C) = (pT−pB)(pT−pC)−pT(1−pT)/(nT−1) (bias term on the target only; **< 0 ⇒ admixture
    in T**); **f4**(A,B;C,D) = (pA−pB)(pC−pD) (no bias term for four distinct demes); **D-statistic**
    (ABBA-BABA) in the symmetric form D = Σ(ABBA−BABA)/Σ(ABBA+BABA), sharing f4's per-site kernel
    with opposite sign. Frequencies are sample derived-allele frequencies; n = 2×sampled diploids.
  - **Deme-slot selection (hybrid, chosen with Rob).** f2 is symmetric → always all unordered pairs.
    f3/f4/D are driven by the `fstats_demes` string param: empty ⇒ **auto** (every deme as an f3
    target × every other pair; every 4-deme quartet in its 3 independent f4/D pairings); `"T;B,C"`
    ⇒ a single f3; `"A,B;C,D"` ⇒ a single f4 + its D. (In a CSV param file the comma-bearing value
    must be quoted, e.g. `"0,1;2,3"`.) One uniform CSV `<model>_fstats_run%d_year%d.csv`
    (Stat,DemeA..D,Value,NumSites).
  - **Tests** (`fstats_test.go`): exact per-site kernel checks + end-to-end fixtures — complete
    differentiation ⇒ f2 = 1; midway target ⇒ f3 = −1/3 (admixture); clean ABBA ⇒ f4 = −1, D = +1;
    auto-enumeration counts (4 demes ⇒ 6 f2 / 12 f3 / 3 f4-D); explicit-spec and malformed-spec
    (falls back to f2-only) guards. **Verified end-to-end** (QuickTest base, seed_style=1,
    init_heterozygosity=0.1): 3 demes/m=0.05 → auto 3 f2/3 f3/0 f4; 4 demes/m=0.02 + `"0,1;2,3"`
    → explicit 6 f2/0 f3/1 f4-D; fixed `rng_seed` → byte-identical.
- [x] **Joint / multi-population SFS (2D/3D).** Current SFS is single-population; the OoA model is
  fit to the *joint* SFS of African/European/Asian samples (dadi, fastsimcoal, momi).
  **Landed 2026-07-06.** (`pkg/analysis/joint_sfs.go`: `ComputeJointSFS`/`SaveJointSFS`, behind
  `track_joint_sfs`.) Generalizes sfs.go's single-population derived-allele counting to a per-deme
  matrix: a cell (i,j[,k]) counts sites at which each axis deme has that many derived alleles in
  its sample; only sites segregating in the combined axis-deme sample are counted. Axis demes from
  `joint_sfs_demes` (2 or 3; empty ⇒ first eligible demes). Long-format CSV (one row per non-empty
  cell) that loads directly into numpy/dadi. Tests (`joint_sfs_test.go`): 2D/3D cell placement,
  explicit-deme selection, single-deme guard, CSV shape. Verified end-to-end (3D over 3 demes).
  - **Deme inheritance is now PATERNAL by default** (was maternal), parameterized via
    `deme_inheritance` (`paternal`|`maternal`) in createChild — patrilocal island model, consistent
    with the Y-line. No effect when num_demes=1 (all demes 0), so the single-population no-op holds.

### 6c. Run their model in your engine — goals 1, 2

- [x] **Implement canonical OoA demographic models as scenarios** — Gutenkunst 2009 / Gravel 2011
  three-population model (bottleneck, split times, growth, migration). Lets DRIFT *run* the
  standard model and show exactly what it predicts. **Landed 2026-07-06.**
  - **Scenario format (chosen with Rob): a JSON demography file + a per-year scheduler.**
    `pkg/demography` parses a `demography.json` authored in the *published* units of the field
    (diploid Ne, times in generations, per-generation migration rates) as an ordered list of
    **epochs**; each epoch declares its active demes (with optional exponential `growth`), an
    optional `split`, and the pairwise `migration` matrix. `Scenario` implements a new
    `core.DemographyScheduler` interface; `manager.LoadModel` attaches it when a model ships a
    `demography.json` (user-config override, else the base model's). The run loop calls
    `Apply(model,pop)` once per year (before Birth). Same mechanism will run a Babel model.
  - **Forward-time / units mapping.** DRIFT ticks in calendar years; the published models are in
    WF generations. Two header knobs bridge them: `generation_time` (yrs/gen, default 25) and a
    **rescale `Q`** (default in the shipped file = 10) — the standard sim rescaling (divide every
    Ne and generation count by Q, multiply rates by Q). So an epoch of D generations spans
    `round(D/Q * generation_time)` years, Ne→census cap `round(Ne/Q)`, and per-gen migration `m`→
    per-year prob `m*Q/generation_time`. Q is the single speed⇄fidelity knob (Q=1 = full scale).
  - **Splits = relabel a random subset** (chosen with Rob, over seeding fresh founders): a split
    carves the child deme out of the parent by relabelling `round(child_founders/Q)` random living
    members (capped at half the parent), preserving shared ancestry so Fst/f-stats/joint-SFS show
    the real tree. All RNG consumed in sorted-id order → runs stay byte-reproducible.
  - **Engine integration.** `core.Model` gained `DemographyScheduler` + a `DemeCaps` runtime map.
    `death.go` culls **per-deme** to `DemeCaps` when a scenario is active (supersedes the global
    `max_pop_size`/`max_growth_rate`); `mating.go` mates within each *present* deme and skips its
    scalar island-migration (the scheduler owns migration). Both are gated so non-demography runs
    are a **strict no-op** (byte-identical — verified). No checkpoint schema bump: the scheduler is
    rebuilt from the JSON on load and caps are recomputed each year (nothing new serialized).
  - **Shipped model:** `static/basemodels/OoA/demography.json` = Gutenkunst 2009 OutOfAfrica_3G09
    (N_A 7300 → N_AF 12300; OOA bottleneck N_B 2100; EUR 1000·e^{0.004t}, ASN 510·e^{0.0055t};
    m_AF-B 25e-5, m_AF-EU 3e-5, m_AF-AS 1.9e-5, m_EU-AS 9.6e-5; T_AF/T_B/T_EU_AS = 8800/5600/848
    gen). At Q=10 the whole history is 22,500 DRIFT years. `OoA/parameters.csv` set up to drive it
    (num_demes=1 at setup → scheduler activates demes; init_heterozygosity=0.1 + seed_style=1 for
    standing variation; track_Fst/fstats/joint_sfs on; generation_time=25). `generation_time` added
    to `parameter_defaults.csv`.
  - **Tests** (`pkg/demography/demography_test.go`): timeline/epoch-boundary math, cap rescaling +
    exponential growth (+clamp past present), per-gen→per-year migration conversion (+clamp),
    deterministic capped split, symmetric migration, `Apply` split-firing + cap-setting, and a load
    of the shipped OoA file (4 epochs, both splits, 22,500 yrs). **Verified end-to-end** with a
    compact scenario (`users/smoke/models/OoAtest`, generation_time=1/Q=1, ~130 yrs, fixed
    rng_seed=777): both splits fire, 3 demes persist, and the stats recover the OoA topology
    (AFR,(EUR,ASN)) — Fst(AFR-EUR)≈Fst(AFR-AS) > Fst(EUR-AS); all f3 > 0 (no admixture); 3D joint
    SFS populated. Two runs byte-identical; a non-demography model (QuickTest) attaches no scenario.
    CSV gotcha still applies (comma-bearing string params like `fstats_demes` must be quoted).
- [x] **Serial founder effect vs. single-origin dispersal.** The flagship OoA signature is the
  near-linear decline of heterozygosity with distance from Addis Ababa (Ramachandran 2005,
  Prugnolle 2005). Tests whether a single-origin (Babel) dispersal produces the *same* gradient —
  a direct head-to-head, run in the same engine under the same fixed seed. **Landed 2026-07-06.**
  - **Distance axis = serial-split RANK (chosen with Rob), not geography.** The demography engine
    is aspatial (abstract demes; `doSplit` relabels `individual.Deme` only, never Lat/Lon), and the
    published gradient's geographic distance is itself just a proxy for the *number of cumulative
    founder events*. So the causally-correct, spatial-structure-free axis is each deme's serial-
    founder rank: origin demes = rank 0, and a split's child = parent rank + 1. Exposed via
    `demography.Scenario.DemeRanks()` (a forward pass over epoch splits); the analysis reads it
    through an OPTIONAL `demeRanker` interface (type-asserted), so `core.DemographyScheduler` is
    unchanged and mocks/ordinary schedulers need not implement it. No geographic/`Wander` machinery
    added — that would be the follow-up if a literal km-from-origin axis is ever wanted.
  - **Founder-size decay = a Babel-scenario PARAMETER (chosen with Rob), not a scheduler feature.**
    A decaying serial chain is just successive split epochs with smaller `child_founders` in the
    JSON — exactly how OoA specifies 2100/510. Zero new scheduler code; the model stays an auditable
    published-units file. `doSplit`'s existing half-parent cap naturally reinforces the founder
    effect on the deeper splits.
  - **Analysis** (`pkg/analysis/hetdistance.go`: `ComputeHetDistance`/`SaveHetDistance`, behind
    `track_het_distance`). Per deme, over the *shared* pooled-polymorphic site set (the Fst NumSites
    convention, so demes share a denominator): **expected heterozygosity He** = Nei unbiased gene
    diversity mean_sites (2n/(2n-1))·2p(1-p), and **observed het Ho** = mean_sites het/n (reusing
    fst.go's `demeMembers`/`countSite` + vcf.go's `buildContigs`/`bitAt`). Each deme tagged with its
    rank; an OLS line He~rank is fitted — its **slope is the head-to-head number** (diversity lost
    per founder step), with R² for linearity. Writes two CSVs: `<model>_hetdistance_*` (per-deme
    table) and `<model>_hetgradient_*` (one-row slope/intercept/R²/HasRanks summary). Params
    `track_het_distance`/`het_distance_sample_size` in `parameter_defaults.csv`.
  - **Shipped Babel model:** `static/basemodels/Babel/{demography.json,parameters.csv}` (registered
    in `registry.csv`) — one origin (Ne 10000) dispersing sequentially 0→1→2→3→4 with decaying
    founders (2500→1500→900→500) and nearest-neighbour stepping-stone migration; Q=10, gt=25 to
    match OoA, 4500 DRIFT years. `track_het_distance=1` also added to `OoA/parameters.csv` so both
    emit the gradient.
  - **Tests** (`demography_test.go` `TestDemeRanks` — OoA-shaped 0/1/2 and serial 0/1/2/3;
    `hetdistance_test.go` — declining-gradient fixture with slope exactly −1/3, no-scheduler
    label-fallback, single-deme guard, CSV shape). **Verified end-to-end** with compact
    `users/smoke/models/{OoAtest,Babeltest}` (Q=1, gt=1, fixed rng_seed=777): both scenarios run in
    the same engine, all splits fire, ranks resolve (`HasRanks=true`), and both emit comparable
    per-deme He and an He~rank slope. Reruns byte-identical (determinism intact); non-demography
    models leave the flag off (strict no-op). NOTE: the compact smoke runs are too short/small to
    develop a strong gradient (sampling noise dominates the smallest deme) — the shipped base models
    (Babel 4500 yr / OoA 22,500 yr at Q=10) are the ones that produce the research-grade decline.

### 6d. Effective population size & the bottleneck question — goals 1, 3

- [ ] **LD-based recent-Ne estimation (GONE-style)** and an **Ne-through-time output.** Creationist
  history posits recent severe bottlenecks (Flood, Babel); PSMC/MSMC/GONE inferences of these are
  actively debated. From known simulated truth, DRIFT can show what these methods get right/wrong.

### 6e. Molecular-clock / dating layer — goals 1, 3

- [ ] **Coalescence → calendar dates under explicit, swappable mutation-rate / generation-time
  assumptions**, with a **sensitivity analysis** showing how the famous dates (~200 kya) move as
  assumptions change. Builds on existing Y-Adam/Mt-Eve generations-back output. Makes the
  timescale dependency quantitative rather than rhetorical.

### 6f. DFE & genetic load — goal 1

- [ ] **Literature-anchored distribution of fitness effects** (gamma/Weibull mixtures — already
  using Weibull), with genetic-load and mutational-meltdown tracking across timescales. Connects
  to the dominance/epistasis items in §1.

### 6g. Haplotype & selection statistics — goal 2

- [ ] **EHH, iHS, XP-EHH, and LD-decay curves.** Standard selection-scan / haplotype-structure
  stats; XP-EHH is cross-population. Useful for realism and for testing selection claims tied to
  the OoA expansion.

### 6h. Credibility backbone (do regardless) — goal 3

- [x] **Validation suite against analytic neutral expectations.** Demonstrates the engine
  reproduces (or, where it deviates, transparently *characterizes*) textbook neutral results
  before any creationist argument is made. **Landed 2026-07-07.**
  - **KEY ARCHITECTURAL FINDING.** DRIFT has two independent genetic-bookkeeping systems, and
    only one supports the θ = 4Nµ equilibrium: the **Chromosomes bitfield** (read by
    SFS/VCF/Fst/joint-SFS/het-distance) holds only FOUNDER standing variation seeded at t=0 and
    receives **no de-novo input** during a run — it only drifts to loss/fixation, never reaching
    mutation-drift equilibrium. De-novo mutations enter the **mutation pool** (`pop.IndMutations`
    + `pop.MutationPool`, gated by `track_mutations`), which has continuous Poisson(mu) input and
    hence the equilibrium the neutral expectations require. §6h therefore computes its SFS/θ/
    Tajima's D/HWE from the **mutation pool**, not the bitfield.
  - **Statistics** (`pkg/analysis/neutral.go`: `ComputeNeutralStats`, `SampleLiving`, expected
    1/i SFS + `SFSShapeChiSquare`, `BuildReport`/`SaveNeutralValidation`). Each distinct mutation
    id is one infinite-sites segregating "site"; derived count = carrying strand copies in the
    sample; het iff exactly one strand carries it. Genome-wide θ_W = S/a_n, θ_π = Σ 2d(2n−d)/
    (2n(2n−1)); Tajima's D / Fu-Li reuse `stats.go`. Unit tests (`neutral_test.go`) verify every
    statistic on a hand-built pool (exact θ_π=1.25, θ_W=3/a₇, D=0.3304, F_IS=0.3143, SFS bins).
  - **End-to-end harness** (`pkg/validation/`, `RunNeutral` + `-validate` CLI). Drives the REAL
    engine (Birth/Mating/Death, mutation kernel, xoshiro RNG) through a strictly neutral scenario
    (`track_mutations=1`, `f_neutral=1` ⇒ every mutation neutral, `track_DNA=0`, `num_demes=1`)
    over replicate seeds; averages to tame Monte-Carlo noise and pools the SFS across replicates.
    θ equilibrates by ~500 y; default burn-in 1500 y, 12 replicates, N=120, mu=1.
  - **WHAT DRIFT ACTUALLY DOES (finding, decided with Rob).** Two textbook expectations hold and
    are asserted tightly: **HWE** (F_IS ≈ 0) and **θ recovery** (θ_W and θ_π both recover a
    stable, reproducible θ ⇒ implied **Ne = θ_W/(2mu)** consistent across seeds; Ne ≈ 0.15–0.20 N
    — emergent because DRIFT is overlapping-generations/monogamous/age-structured, NOT WF). But
    DRIFT does **not** reproduce the WF **Tajima's D ≈ 0**: under strict neutrality it robustly
    yields **D ≈ −0.55 to −0.9** and **θ_π/θ_W ≈ 0.75** (a systematic *excess of rare variants*)
    across every sustainable life history tested — a real consequence of its high reproductive
    variance (monogamy ⇒ many zero-offspring individuals; high-fecundity birth-then-cull ⇒ a
    repeated within-generation expansion/contraction that skews the coalescent). Per Rob's call,
    the harness **characterizes** this reproducible baseline (regression guard) rather than
    faking D ≈ 0; any D/SFS-based inference on DRIFT output must account for the skew.
  - **Runnable modes.** `drift -validate` runs the replicated harness, writes
    `neutral_validation.csv`, and exits non-zero on failure (CI gate). Any normal
    mutation-tracking run with `track_validation=1` also emits a single-realization diagnostic
    (`<model>_validation_*.csv`). Shipped **`Neutral`** base model (`static/basemodels/Neutral/`,
    registered) + a compact `users/smoke/models/Neutraltest` fixture (verified e2e: PASS,
    D=−0.94, F_IS=−0.11, Ne/N=0.17). e2e reproducibility + all checks tested in
    `pkg/validation/validation_test.go`.
  - **Also fixed** a latent determinism bug in `SampleIDs`/`SampleLiving`: the id slice was built
    from randomized map-iteration order and then shuffled, so a subsample differed between runs
    under a fixed seed. Now sorted before shuffling (subsample is a deterministic function of the
    RNG stream). Mutation-pool sums also iterate ids in sorted order (float addition isn't
    associative) so θ_π/F_IS are bit-reproducible.
- [ ] **ABC-readiness.** Once standard summary stats are emitted, DRIFT slots into Approximate
  Bayesian Computation — reframing it from "a simulator" to "an inference engine that can fit any
  model, including ours."

## 7. Emerging frontiers / cutting-edge bets

Where the field is heading (~2025–2026) and where DRIFT's architecture gives it a structural edge.

**The key insight:** DRIFT is simultaneously forward-time, individual-based, pedigree-tracking,
*and* spatially explicit. That combination natively produces — as *known ground truth* — the three
objects the modern field spends enormous effort merely *estimating*: the true ARG (genealogy),
time-stamped + geo-located ancient samples, and true IBD segments. Coalescent tools run backward
and struggle with serial/ancient samples and geography; inference tools estimate these with
uncertainty. **DRIFT knows the answer.** The durable edge is not out-computing inference tools but
being the ground-truth generator that exposes how much their deep-time conclusions depend on their
priors. (SLiM + tree-recording has similar *capability* — the white space is the framing and the
single-origin / young-population question, which no one is pursuing from this angle.)

### Field shifts to be aware of (context, not tasks)

- **Tree sequences / ARGs are the new substrate** (`tskit`, `tsinfer`, `Relate`, `ARG-Needle`,
  `SINGER`) — the field moved from summary stats to inferring the full genealogy. See §6a tskit.
- **Ancient-DNA time-series is rewriting human history** (Reich lab AADR, 10,000+ genomes). The OoA
  story is now an *aDNA* story — models get judged against time-transects, not just modern variation.
- **Simulation-based (neural) inference is eating ABC** — the simulator becomes the training-data
  generator. Raw genotype/tree output matters more than hand-computed statistics.
- **IBD-based inference of recent history** (`hap-IBD`, `iLASH`, Refined IBD) sees the recent
  timescale SFS methods are blind to.
- **Mutation-rate crisis** — direct pedigree rates (~1.2×10⁻⁸) are ~half the phylogenetic rate and
  vary across populations/time; unsettled in the field itself, and it compresses clock dates.
- **Pangenomes / T2T / structural variation** — moving past the SNP to complete genomes and graphs.

### Cutting-edge bets (ranked by impact-per-effort)

- [~] **2. Long-IBD as a falsifiable signature of recent common ancestry. [HIGHEST IMPACT/EFFORT]**
  A young, single-origin population predicts a distinctive *excess of long IBD segments* and recent
  TMRCA — exactly the regime classic SFS methods can't see. Compute *true* IBD from DRIFT's real
  pedigrees + recombining genomes and compare to what `hap-IBD` infers from real human data. The
  cleanest, most underexploited testable prediction the framework makes.
  - [x] Level 1 prototype: `pkg/analysis/ibd.go` (`ExtractIBD`/`SaveIBD`), wired into `drift.go`
    behind the `track_IBD` flag (default off). Measures co-inherited founder-derived segments;
    emits a segment table + log-binned length distribution. Builds clean.
  - [ ] Validate the length distribution against a known run; convert bit-lengths → cM (needs §6a
    recombination map).
  - [ ] Level 2 (true IBD): unique per-founder-haplotype labels threaded through meiosis, or
    tree-sequence recording — removes the over-coarsening noted in `ibd.go`.
  - [x] Fixed a latent bit-order bug found while building IBD: `CountContiguousBlocks`
    (`pkg/utils/genetics.go`) used MSB-first `%064b` strings while bits are set LSB-first, so any
    haplotype block straddling a 64-bit word boundary was over-counted — inflating
    `IndData[NumBlocks]`. Rewritten to scan genome positions with the canonical convention;
    regression test in `pkg/utils/genetics_test.go`. NOTE: previously-saved `NumBlocks` values are
    biased high and not comparable to new runs.
- [ ] **1. Known-truth ARG, stress-tested for assumption-dependence.** Emit a true tree sequence
  (depends on §6a tskit), run `Relate`/`tsinfer`/`SINGER` on DRIFT's genotypes, quantify the
  true-vs-inferred gap, and show how inferred dates/Ne move as clock and generation-time priors
  vary. Novel, publishable, on-trend.
- [ ] **3. Ancient-sample time-series output.** Emit samples at specified generations and map
  locations, formatted to compare against real aDNA transects. Fights the OoA debate on its actual
  current frontier. Builds on forward-time + `Lat`/`Lon`.
- [ ] **4. Simulation-based-inference data generator.** Position DRIFT as the labeled-data engine
  for neural SBI under alternative (young) models. Future-proofs DRIFT as ABC fades; no one in the
  creationist space is doing it. (Relates to §6h ABC-readiness.)
- [ ] **5. Mutation-rate-honest dating.** Extend the §6e clock layer to use direct pedigree rates
  and explicitly propagate rate uncertainty into dates — capitalizing on a discrepancy the
  mainstream itself has not resolved.

> Defensive note: if DRIFT doesn't speak tree-sequence it becomes an island, and raw
> genotype/tree output is becoming more important than reimplementing every statistic. Lean toward
> emitting VCF + tree sequences and letting standard tools compute stats. This is why §6a tskit was
> promoted to near-term.

## 8. Experiment infrastructure

Two features that turn DRIFT from "run a sim" into "run controlled experiments." Both build on the
Phase 0 determinism work.

### 8a. Pluggable phase modules (registry + dropdown) — DONE (2026-06-25)

Let users select the algorithm for each simulation phase (birth, death, mating, seeding, setup)
by name, and add their own, rather than editing the call site. Chosen over the "edit a line to
call `BobsBirth`" alternative because the selected algorithm is **recorded in the model params**
(provenance — an edited call site leaves no trace in the output, undermining Phase 0), is
**GUI-selectable**, and **composes with forks** (§8b).

Mechanism (Go can't enumerate functions by name at runtime, so a name prefix alone can't drive a
dropdown): a **registry** per phase. Each module self-registers in `init()` under a string key;
the dropdown lists the registry keys; the choice is stored as a `*_style` string param and
dispatched at call time. File/function naming (`birth_*.go`) is human convention, not the
discovery mechanism.

- [x] String-parameter support on the model: `core.Model.StringParams` + `StringParam(key, def)`;
  `SetParameter` and `LoadParameters` route non-numeric values there (`Parameters` stays float64).
- [x] Birth/Death registries + `"standard"` default modules + dispatchers + `Available*Styles()`
  (`pkg/simulation/module_registry.go`). Existing `Birth`/`Death` became `birthStandard`/
  `deathStandard`, self-registered via `init()`; the exported `Birth`/`Death` are now dispatchers
  keyed on `birth_style` / `death_style`. Defaults added to `parameter_defaults.csv`. Unit tests
  in `module_registry_test.go`. Verified end-to-end: normal run dispatches; bogus style →
  `unknown birth_style "bogus"; available: [standard]`.
- [x] Validate the selected style at model load — `ValidateStyles(model)` called in `drift.go`
  (fails fast with an "available: [...]" message).
- [x] Unified all phases in a new `pkg/modules` package (imports only `core`, so any package can
  register without a cycle). Birth/Death register from `simulation`, mating from `simulation`,
  seeding from `events`, setup from `config` — each via `init()`. `Birth`/`Death`/`Mating`/`Seed`/
  `InitializePop` became thin dispatchers.
- [x] Legacy **integer** `mating_style` / `seed_style` codes still resolve (0/1/2 → names) via
  `resolveStyle`; new models can use string names, which take precedence. Setup stays
  scenario-driven with an optional `setup_style` override. Verified: legacy-int and string-style
  runs both dispatch with zero fallback warnings.
- [x] Webserver endpoint `GET /api/modules/list` returns `modules.AllStyles()` (available names per
  phase) for the GUI dropdowns. Verified live: returns birth/death/mating/seed/setup lists.
- [x] Frontend: `birth_style` / `death_style` dropdowns in the model editor ("Simulation Modules"
  group), options fetched from `/api/modules/list` (`loadPhaseModules` in `web/static/app.js`).
  String params round-trip through save/load: `handleModelSave` routes strings to `StringParams`,
  `SaveParameters` persists them, `handleModelLoad` merges them back. Verified live: load returns
  the styles, save persists without corrupting the CSV. (The already-working integer
  `mating_style` / `seed_style` dropdowns were left as-is; making them registry-populated too is
  optional polish.)
- [ ] (Future) per-module parameter schemas; optional embedded scripting (Lua/Starlark) for true
  no-recompile user algorithms — deferred (hot-path cost + complexity).

### 8b. Checkpoint / save / resume / fork — CORE DONE (2026-06-25)

Stop a run and serialize full state; resume it exactly, or **fork** it into N divergent
continuations. Enables the core experimental design: run to year T (e.g. a Flood/Babel
bottleneck), checkpoint that population as a **baseline**, then vary one thing (a param, a seed,
or — via §8a — an algorithm) across forks while holding all shared history constant.

- [x] `pkg/checkpoint`: versioned `Checkpoint` (schema v1) serializes `Pop` (all maps) +
  `FreeParameters` counters + RNG state via `encoding/gob`; `Save`/`Load`/`Restore`. Model config
  is reloaded from files on resume, then the checkpoint is overlaid (state only — no `Model.Map` /
  animation state). `ensurePopMaps` guards against nil-map panics on resume. Round-trip +
  exact-stream + fork-divergence unit tests in `checkpoint_test.go`.
- [x] `RNGState() []byte` / `RestoreRNG([]byte)` in `pkg/utils/rng.go`. Required swapping the shared
  generator's source to a small custom xoshiro256** (seeded via splitmix64) because math/rand's
  default source is not serializable in Go 1.23; state is 4×uint64 = 32 bytes. Determinism tests
  still pass.
- [x] CLI: `-checkpoint-out`, `-checkpoint-interval`, `-checkpoint-in`, `-fork-seed`. Run loop
  refactored into a `runOne(run, startYear, pop)` closure so a fresh multi-run and a single resumed
  continuation share one body. Verified end-to-end: exact resume is byte-identical to an
  uninterrupted run (years ≥260 matched); a fork (`-fork-seed=999`) diverges from the same
  checkpoint.
- [x] **Named saves** (`pkg/checkpoint/named.go`): a finished run can save its end state under a
  name into the model's `saves/` directory (`-save-state=<name>`), reload it to start other runs
  (`-load-state=<name>`, combine with `-fork-seed`), and list them (`-list-saves`). API:
  `SaveNamed`/`LoadNamed`/`List`/`SavePath` (names contained to the saves dir). Verified end-to-end:
  normal run saves "baseline"; loading it continues byte-identically to the reference run; forking
  it diverges. Tests in `named_test.go`.
- [x] **GUI save/load/fork.** A "Run State" panel in the model editor (`web/static/index.html` +
  `app.js`): save the end-of-run state under a name, pick a saved state to start from, and set a
  fork seed. Endpoints `GET /api/saves/list` and `POST /api/saves/delete`; `save_state`/`load_state`/
  `fork_seed` thread through `StartSimulationRequest` → `SimulationJob` → `runSimulation` CLI args.
  Saved states are anchored to the model (`users/<user>/models/<model>/saves`) via
  `SavesDir`/`ModelFor` regardless of the run's output dir. Verified live: list/delete round-trip;
  a CLI save appears in the canonical dir the endpoint reads.
- [x] **Graceful stop-and-save of a running sim.** File-based signal (`pkg/checkpoint/stop.go`:
  `RequestStop`/`StopRequested`/`ClearStop`) — chosen over OS signals, which are unreliable on
  Windows. The webserver's new `StopSaveJob` drops a `.stop` (with a save name) into the run's
  output dir; the sim checks it each year, saves state under that name, and finishes the run cleanly
  (results + analyses still written) — no hard kill. GUI: a "Stop & Save" button on the Progress tab
  (`/api/simulation/stop-save`); the forceful `Cancel` (hard kill) remains as a fallback. Added a
  `stopping` job status. Verified end-to-end: sim honors the stop, saves at the current year, clears
  the flag, and the save appears in the model's saves list. Tests in `stop_test.go`. (Also works
  from the CLI: `touch <output-dir>/.stop`.)
- [ ] Resume currently continues a single run; extend to resume/fork within a multi-run batch if
  needed.

### 8c. Batch / experiment runner (auto-run many models & configurations)

Automatically launch many runs from one spec, so a whole experiment is a single action instead of
dozens of manual starts. Ties together replicates (same params, many seeds → a distribution),
parameter sweeps (varying params → a grid), multiple models, and §8b saved states/forks. This is
the applied form of §4's "parameter sweeps / batch experiments" and the natural payoff of the
determinism (§0) + checkpoint/fork (§8b) work.

Requested capabilities:
- [ ] **Run many models in one batch** — a list of models, each started either **from scratch** or
  **from a named saved state** (§8b) as the common baseline.
- [ ] **Replicates** — run the same model+parameters N times with distinct, deterministic seeds
  (`baseSeed + index`) to build a distribution across stochastic realizations.
- [ ] **Parameter sweeps** — expand a grid/list of parameter values (e.g. `mu × max_pop_size`) into
  one run per combination; optionally crossed with replicates (combinations × replicates).
- [ ] **Fork sweeps** — from one saved baseline, fan out runs that vary a single thing (seed, a
  parameter, or a phase module from §8a) while holding all shared history constant.

Design notes:
- An **experiment spec** (CSV/JSON, or a GUI form): models, source (scratch | saved-state name),
  replicate count, and the sweep grid. Expand it into N queued jobs.
- Reuse the existing web **queue** (`pkg/webserver/queue.go`) — a batch is just many `AddJob`s; or a
  CLI batch runner for headless use. Each job gets a unique output dir + recorded seed/params.
- **Aggregate output**: a cross-run summary table (one row per run: params, seed, key final stats)
  plus per-run result CSVs, so a sweep is analyzable without hand-collating.
- **Provenance**: every run records its resolved params + seed (+ source saved-state, + git SHA) —
  see the run-manifests item in §5. Essential for reproducibility of the whole experiment.
- Concurrency/throughput: the queue currently runs jobs one at a time; a batch of hundreds may want
  a worker pool (bounded parallelism) — note the determinism guarantee holds because each run is a
  separate process with its own seed.

---

### Suggested starting order

1. Per-run deterministic RNG (§0) — makes all later results trustworthy and reproducible.
2. Dominance-aware fitness (§0 / §1) — unlocks the founder-effect / bottleneck research questions.
3. Fitness → survivorship (§2) and a test harness (§5) so the above can be validated.

**OoA-testing thread (§6), in order:** VCF export (§6a) → Fst + f-statistics + joint SFS (§6b) →
implement the standard OoA model as a scenario (§6c) → validate against neutral expectations
(§6h). This chain lets DRIFT run the mainstream model and a Babel model in the *same engine*,
score both with the *same* statistics the field uses, on data comparable to *real* human
genomes — the comparison a critic can't dismiss on methodological grounds.

**Cutting-edge thread (§7):** the enabling substrate is tree-sequence (tskit) output (§6a,
now near-term). The sharpest standalone first bet that does *not* require tskit is long-IBD of
recent common ancestry (§7 bet 2) — it rides on true pedigrees + recombining genomes DRIFT
already has, and yields a clean falsifiable prediction. tskit then unlocks bets 1 and 3.
