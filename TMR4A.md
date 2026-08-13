# TMR4A in DRIFT — implementation specification

**Status:** **step 1 of §7 landed** on `rob-TMR4A-1` — W1 + W2 (focal mode) + W5 + the W8 tests that
guard them. Ground-truth TMR-*K*-A now runs end-to-end on ordinary diploid runs. W3, W4, W6, W7 and
W9 are still design only. See §8 for the implementation record.
**Scope:** what must be built to (1) compute *ground-truth* TMR4A inside DRIFT, and (2) calibrate
ARGweaver's *inferred* TMR4A against that truth.
**Roadmap ties:** README items 5a (four-alleles test) and 5b (ARGweaver); TODO §6a (real-data
interop), §7 (tree-sequence substrate).

---

## 0. The assumption this module exists to expose

TMR4A as published (Swamidass, on ARGweaver-inferred ARGs over 1000 Genomes) is the time at which
the marginal genealogy at a locus has exactly four lineages. Reading that as "the earliest date a
founding couple is consistent with" smuggles in two premises:

1. **Universal common ancestry of alleles.** Every allele at a locus descends from one ancestral
   allele; all observed differences accumulated by mutation since. The marginal tree therefore has
   a root, and every pair of alleles has a finite coalescence time.
2. **Four is the ceiling.** Two diploid parents transmit at most four alleles per autosomal locus.

Neither is entailed by a created origin.

**On (1).** If Adam and Eve were created carrying distinct alleles, those alleles are not related by
descent at all. No mutational path connects them and no time exists at which they coalesce. Fitting
a coalescent model to created divergence does not recover a date — it converts *sequence
differences* into *apparent time* through the assumed mutation rate. A created allele pair differing
at *d* sites is read as ~*d*/(2μL) generations old whether or not a single generation elapsed. This
is the mechanism by which TMR4A lands near half a million years, and it is measurable: DRIFT can
generate data whose true TMR4A is known and small, then see what ARGweaver reports.

**On (2).** Four is the ceiling only if each founder's germline is a faithful copy of a single
diploid genome. If the reproductive cells of the created couple carried variation among themselves —
created germline heterogeneity, not somatic diploidy — then the founding pair can transmit **more
than four alleles per locus**, and the "4" in TMR4A is simply the wrong bound to test.

**Consequences for the build, and they are structural, not cosmetic:**

- The module computes **TMR-*K*-A** with *K* a model parameter. *K*=4 is one case, not the design
  centre. Cost of the generalisation: zero — the walker emits the full lineage-count-vs-time
  trajectory and any *K* is a lookup against it.
- The backward walk must **terminate at created founders and say so**, rather than forcing a root.
  Three distinct outcomes must be reportable, never collapsed: *coalesced to K*,
  *censored at founding with >K lineages*, *K never approached*.
- The headline falsifiable output is not a date. It is **the number of distinct *created* alleles
  surviving in the sample at each locus.** If that number is ≤ the founding couple's transmissible
  allele count, a two-person origin is consistent with the data at that locus regardless of what an
  inferred TMR4A says.
- Founder alleles need **identity labels, not just states**. Identity-by-state and
  identity-by-descent diverge the moment created alleles can coincide at a site by chance, and the
  whole argument turns on telling them apart. See W3.
- DRIFT must be able to **represent germline heterogeneity in the founders**, which it currently
  cannot — founders are strictly diploid and meiosis draws from exactly two strands. See W4.

---

## 1. What DRIFT has today

Verified against the tree on `rob-working`.

| Capability | State | Where |
|---|---|---|
| Per-run deterministic RNG | Present (§0) | [pkg/utils/rng.go](pkg/utils/rng.go) |
| Founding-couple scenario | Present, 3 couples, long lifespans — clone as the Eden template | [pkg/methods/setup_pop_flood.go](pkg/methods/setup_pop_flood.go) |
| Founder genome seeding | Per-site independent alt-allele draws at `init_heterozygosity` | [pkg/methods/seeding_population.go:53-80](pkg/methods/seeding_population.go#L53-L80) |
| Recombination, real cM maps | Present (§6a) | [pkg/simulation/birth.go:246-342](pkg/simulation/birth.go#L246-L342), [pkg/core/recomb.go](pkg/core/recomb.go) |
| De-novo mutation, classes, variable rate | Present (§1) | [pkg/utils/mutation.go](pkg/utils/mutation.go), [pkg/utils/mutation_class.go](pkg/utils/mutation_class.go) |
| Pool/bitfield co-segregation | Present under `linkage_model="arm"` | [pkg/utils/mutation.go:55-90](pkg/utils/mutation.go#L55-L90) |
| VCF export | **Founder bitfield only** — de-novo pool is not written | [pkg/analysis/vcf.go:159-190](pkg/analysis/vcf.go#L159-L190) |
| VCF import (real 1000G/HGDP) | Present, loads into the bitfield | [pkg/analysis/vcf_import.go](pkg/analysis/vcf_import.go) |
| Y / mt coalescence | Present | [pkg/analysis/coalescence.go](pkg/analysis/coalescence.go) |
| **Diploid pedigree** | **Absent** | see below |
| **Segment provenance (ARG)** | **Absent** | see below |
| **Allele identity labels** | **Absent** — bitfield is biallelic 0/1 | — |
| **Founder germline heterogeneity** | **Absent** — founders are strictly diploid | — |

### The two gaps that matter most

**No diploid pedigree.** [birth.go:160-174](pkg/simulation/birth.go#L160-L174) writes
`MaleDB[child]` **only when the child is male** (recording the father) and `FemaleDB[child]` **only
when female** (recording the mother). A son's mother and a daughter's father are recorded nowhere.
What exists is the Y-chain and the mt-chain, which is all `FindYAdam`/`FindMtEve` need — and not a
pedigree. The prune rule in [death.go:316-321](pkg/simulation/death.go#L316-L321) (drop a dead node
only when it left no same-sex offspring) is exactly the right retention discipline and should be
generalised rather than replaced. Note `IndData[ind]` is deleted unconditionally at
[death.go:325](pkg/simulation/death.go#L325); birth years for the dead survive only inside the
`Ancestor` records ([pkg/core/types.go:352-355](pkg/core/types.go#L352-L355)).

**No segment provenance.** The meiosis masks are locals at
[birth.go:84-87](pkg/simulation/birth.go#L84-L87), consumed by `meiosis` and `InheritMutations` and
discarded. Which parental strand a given position came from is unrecoverable after the fact.

---

## 1A. Region sizes in the published analysis, and what they imply

Established from the source papers:

- [Rasmussen et al. 2014](https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1004342)
  analysed **54 individuals = 108 haplotypes** (Complete Genomics "69 genomes").
- The genome was partitioned into **~2 Mb blocks, ~1,376 of them**, ARGweaver run on each in
  parallel. This is the *computational* unit — how the job was chunked, not the unit the statistic
  is computed over.
- The **analysis unit is the local tree**: the marginal genealogy over an interval between inferred
  recombination breakpoints. TMR4A is evaluated per local tree and aggregated genome-wide, weighted
  by span. Its size is set by recombination, not chosen by the analyst.

Derived, not quoted — flagged as such wherever it is used downstream. With 108 samples and Ne≈10⁴
the marginal tree's total branch length is ~2×10⁵ generations; at r≈1.2×10⁻⁸/bp/gen that places a
recombination event roughly every **few hundred bp**. So local trees span hundreds of bp to ~1 kb,
and there are **millions** genome-wide, not thousands.

**The power consequence.** Human μ and r are nearly equal (≈1.25×10⁻⁸ vs ≈1.2×10⁻⁸ per bp per
generation), so each non-recombining block carries **on the order of one segregating site**. Per
local tree the data is almost devoid of information: the inferred depth comes from the coalescent
prior with its assumed Ne history, plus information borrowed along the sequence via the SMC. The
genome-wide median is therefore not millions of independent estimates but a heavily correlated,
prior-shaped field. This is the quantitative core of what the calibration study tests.

**Two consequences for this build:**

1. The focal-locus design in W2/W5 is correct as written — a marginal tree at a point is the right
   analogue, and DRIFT computes it exactly rather than inferring it. No change.
2. **A single chromosome arm is a sufficient test object.** ARGweaver's own native unit is a ~2 Mb
   block, so one arm is *larger* than what it was run on. There is no methodological loss in
   scoping the study to one arm rather than a whole genome.

## 1B. The assumed effective population size — a primary object of study, not a nuisance parameter

**Ne is an input to ARGweaver, not an output.** The published human analysis states "we assumed
[N] … generations" — the value is a supplied prior. *(Exact value TBC: it renders as an image in the
PLOS HTML; read it off the PDF and record it here.)*

**The analytic consequence, and it is severe.** Under a constant-*N* coalescent the expected time
until *k* lineages remain is 4*N*(1/*k* − 1/*n*). For *k*=4 and *n*=108 that is **≈0.96 *N*
generations**. So:

> TMR4A ≈ *N* generations, before any data is consulted.

At 25-30 yr/generation and *N*≈10⁴ the prior alone predicts ~250-300 ky — the same order as the
published ~500 ky median. The data adjusts the prior expectation by a factor of order one; it does
not determine it. Combined with §1A's finding that each local tree carries ~1 segregating site, the
statistic is substantially a restatement of its own prior.

**Why constant Ne≈10⁴ is a poor description of recent human history.** It is a deep-time coalescent
abstraction, harmonic-mean-dominated by the smallest sizes in the history, and it is not census
size. Human census size has grown by orders of magnitude in the last ten thousand years, and recent
explosive growth is independently evidenced by the large excess of rare variants. A constant-*N*
model can represent neither recent explosive growth nor a recent founding bottleneck — precisely the
two features that distinguish the models under test. Fitting either history with a constant *N*
does not fail loudly; it silently returns a depth.

**Structure is absorbed as depth.** A structured population of modest census size coalesces more
slowly than a panmictic one of the same size, so it presents as a large panmictic Ne. DRIFT models
structure explicitly and already ships it — demes, migration matrices, barriers, habitat
suitability, geographic Fst (§3) — and measures realised Ne through time (§6d `track_Ne`). So
"known census size + known structure + known recent origin ⇒ what depth does ARGweaver infer?" is
directly testable here, against a true Ne trajectory DRIFT can report rather than assume.

**This promotes a work item (W9) and reorders the plan** — see §7.

**Correction to §1A's derived block sizes.** Those were computed at *N*=10⁴ and inherit its
assumption. Total branch length scales with *N*, so a smaller or younger population produces *fewer*
recombination events per bp — longer local trees, more mutations per block, better-determined
genealogies. The ≤100 bp/bit resolution target in W7 is therefore the conservative bound set by the
mainstream comparison arm; a recent-origin model needs less resolution, not more.

## 2. Work items

Every item is opt-in and must leave a default run byte-identical (the `-validate` gate:
D=−0.6611 / π-W=0.809 / Ne-N=0.313 / SFS-χ² 4.77). Any item that draws RNG on the default path is
wrong by construction. New `.go` files LF on disk.

### W1 — Diploid pedigree retention `track_pedigree` — **LANDED**

Record both parents for every child and retain any node with surviving descendants.

- Add `ParentDB map[int][2]int` (dad, mom) plus a descendant refcount to `core.Pop`
  ([pkg/core/types.go:118-127](pkg/core/types.go#L118-L127)). Exported so gob checkpointing carries
  it ([pkg/checkpoint/](pkg/checkpoint/)).
- Write at [birth.go:160-174](pkg/simulation/birth.go#L160-L174), unconditional on child sex.
- Prune at [death.go:316-325](pkg/simulation/death.go#L316-L325) only at refcount zero, cascading to
  parents. This is the same discipline as the §6f mutation-pool refcount fix — and the same failure
  mode if it is asymmetric, so mirror that structure deliberately.
- Store `BirthYear` alongside, since `IndData` is gone after death.

*Files:* `core/types.go`, `simulation/birth.go`, `simulation/death.go`, checkpoint schema.
*Size:* ~150 lines + tests. *Risk:* unbounded growth on long runs — see §5.

### W2 — Strand provenance `track_arg` — **focal mode LANDED**

Two modes, same consumer:

- **Focal-loci mode** (cheap, do first): for a configured locus set *L*, store per birth per gamete
  which parental strand each *p* ∈ *L* came from — |*L*| bits per gamete. For a handful of loci this
  is ~2 bits per birth.
- **Breakpoint mode** (the real substrate): store the crossover breakpoint list per gamete
  (~23-40 ints), which yields *every* locus. This is the tskit-edge form TODO §7 wants and should be
  the eventual target.

Capture point is [birth.go:84-87](pkg/simulation/birth.go#L84-L87), where the masks already exist.
`createMaskFromMap` ([birth.go:342](pkg/simulation/birth.go#L342)) already knows the breakpoints
before rasterising them to a mask — take them there rather than re-deriving from bits.

*Size:* ~200 lines + tests.

### W3 — Founder allele identity track `track_allele_labels`

The bitfield holds one bit per site: two states, no identity. Two created alleles that happen to
carry the same bit at a site are indistinguishable from two copies of one allele, which is precisely
the IBS/IBD confusion the test cannot afford.

- Per strand, at focal loci only, carry a small integer **founder-allele label** (0..*A*−1 for *A*
  created alleles at that locus) that is inherited through meiosis alongside the bitfield.
- Labels are inert: no fitness, no mutation, no effect on any existing statistic. They exist to let
  the walker answer "same created allele?" exactly.
- De-novo mutations keep their pool IDs as their identity; do **not** key allele identity on
  position, since `GenerateNewMutations` draws `RandIntn(genome_bits)` and positions recur.

*Size:* ~150 lines + tests.

### W4 — Created founder alleles with germline heterogeneity `founder_allele_model`

This is the item that encodes the objection in §0 and has no analogue in the current code.

- **`founder_alleles_per_locus` (*A*)** — number of distinct created alleles at a focal locus.
  *A*=4 reproduces the strict two-diploid-parents case; *A*>4 models created variation among the
  reproductive cells.
- **Distribution across gametes** — how the *A* alleles are apportioned to each founder's gametes.
  Minimum viable: a per-founder pool of alleles with weights, sampled per gamete at birth. This
  breaks the invariant that a parent transmits one of exactly two strands, so it must be a distinct,
  clearly-gated path through meiosis rather than a patch to the existing one.
- **Divergence between created alleles** — the created alleles' sequence difference at the locus,
  set independently of elapsed time. This is the dial that drives apparent-vs-true TMR4A and is the
  single most important parameter in the whole study.
- Restrict the mechanism to founders. Once past generation 0, inheritance is ordinary diploid
  meiosis.

*Files:* new `pkg/methods/seeding_created.go`, new `pkg/methods/setup_pop_eden.go` (clone of
[setup_pop_flood.go](pkg/methods/setup_pop_flood.go), one couple), registered through the §8a module
registry ([pkg/modules/registry.go](pkg/modules/registry.go)).
*Size:* ~300 lines + tests. *This is the largest and least conventional item.*

### W5 — The TMR-K-A walker `analysis/tmrka.go` — **LANDED**

Backward from a sample of 2*N* strands at each focal locus:

- At each step, map (individual, strand) → (parent, parent-strand) via W1 + W2. Two lineages
  arriving at the same (individual, strand) have coalesced.
- Emit the **full lineage-count-vs-year trajectory**; TMR-*K*-A for any *K* is a lookup.
- Terminate at the founder generation and classify each locus:
  - `coalesced` — reached *K* by descent, with the year;
  - `censored_at_founding` — >*K* lineages remained when the walk hit creation;
  - `floor` — the count of distinct **created** alleles (W3 labels) surviving in the sample.
- Also emit TMRCA, the lineage-count at founding, and the created-vs-mutational share of pairwise
  differences at the locus.

No MCMC, no prior, no inference error — this is the true genealogy.

*Size:* ~250 lines + tests.

### W6 — Merged VCF export (**prerequisite for W7**)

[ExportVCF](pkg/analysis/vcf.go#L109) writes the founder bitfield only. ARGweaver needs the
segregating sites *including* de-novo mutations. Overlay the pool onto the bitfield per sampled
strand; they already share one coordinate space (`Mutation.Position` is a genome-bit index) and
co-segregate under `linkage_model="arm"`. Key by mutation ID, not position, for the collision reason
in W3. Optionally emit the W3 labels as an INFO/FORMAT field for downstream truth-vs-inference joins.

*Size:* ~100 lines + tests.

### W7 — ARGweaver bridge

Do **not** reimplement the MCMC (README 5b). Run the real tool externally.

- Converter to ARGweaver `.sites` (+ region `.bed`), from W6 output.
- **Resolution and rate ratio are independent dials in DRIFT.** Multiplying the chromosome CSV's
  arm Start/Length by a scale factor X raises `genome_bits`
  ([pkg/utils/csv.go:205-222](pkg/utils/csv.go#L205-L222)) without touching either rate anchor:
  mutations are a per-individual Poisson count (`RandPoisson(class.Rate)`,
  [pkg/utils/mutation.go:131](pkg/utils/mutation.go#L131)) spread uniformly over `genome_bits`, and
  crossovers per meiosis come from the CSV's separate cM column
  ([pkg/utils/csv.go:174-180](pkg/utils/csv.go#L174-L180)). Scaling by X therefore divides per-bit μ
  **and** per-bit r by X, leaving **θ/ρ invariant**. X buys resolution only.
- **So the ratio is set by two per-generation counts**, not by the genome definition: mutations per
  individual per generation vs crossovers per meiosis. Put both at real human values and the ratio
  is correct at any X. This is the parameter pair to document and defend.
- **Choose X for resolution.** Target enough bits that recombination breakpoints are resolvable and
  each non-recombining block can carry its ~1 segregating site — i.e. bit resolution comfortably
  finer than a local tree span, so on the order of ≤100 bp per bit. Over a single arm of ~50-100 Mb
  that is ~10⁵-10⁶ bits, well inside the range the scale factor reaches. Compute and memory scale
  accordingly; see §5.
- The residual scaling question is therefore narrow and answerable: state the bp-per-bit convention,
  state the two per-generation rate anchors, show θ/ρ matches human values. Not the open-ended
  problem framed in earlier drafts.
- Parser for ARGweaver output → inferred TMR-*K*-A on the same loci, joined against W5 truth.

*Size:* ~250 lines + a written scaling rationale that stands on its own.

### W8 — Validation — **LANDED for W1/W2/W5**

- Walker correctness against a hand-built pedigree with known coalescence times (follow the
  hand-built-panel pattern used for §6b/§6g).
- Byte-identity regression: all flags off ⇒ `-validate` unchanged.
- Determinism: identical seed + different focal-locus sets ⇒ identical genealogy (this is what makes
  locus sweeps safe, see §5).
- A neutral large-*N* run where the true TMR4A is known analytically-ish, as a sanity anchor.

### W9 — Ne-prior sensitivity sweep (**no DRIFT code required**)

The cheapest high-value experiment in the programme, and it can run before any of W1-W8 exists.

- Run ARGweaver on a **fixed** dataset (real 1000G/Complete Genomics, or simulated) under a range of
  assumed effective sizes and size histories, including constant *N* across a wide span, a recent
  explosive-growth history, and a recent-bottleneck history.
- Plot inferred TMR-*K*-A against assumed *N*. The §1B prediction is that it tracks ≈*N* generations.
  If it does, the statistic is largely determined by its prior and that is demonstrable without
  simulating anything.
- Second axis once DRIFT is in play: generate structured populations of known census size and known
  recent origin (§3 machinery), export via W6, and compare ARGweaver's inferred depth against
  DRIFT's true Ne trajectory from §6d `track_Ne`. This separates "large ancient population" from
  "structured recent population", which a constant-*N* prior cannot distinguish.

*Size:* external tooling and run management; a parameter sweep harness, no engine changes.
*Why first:* it is the load-bearing claim of the whole critique, and it is testable immediately.

---

## 3. Parameters to add

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `track_pedigree` | int | 0 | Retain the full diploid pedigree (W1) |
| `track_arg` | int | 0 | 0 off / 1 focal-loci bits / 2 breakpoints (W2) |
| `arg_loci` | string | "" | Focal locus set: bit positions or a range spec |
| `track_allele_labels` | int | 0 | Carry founder-allele identity labels (W3) |
| `founder_allele_model` | string | "diploid" | `diploid` \| `created` (W4) |
| `founder_alleles_per_locus` | int | 4 | *A* — created alleles per locus; >4 = germline heterogeneity |
| `founder_allele_divergence` | float | 0 | Created sequence divergence between founder alleles |
| `tmrka_k` | int | 4 | *K* in TMR-*K*-A |

---

## 4. Output

One row per focal locus per run:

```
locus, k, outcome, tmrka_year, tmrca_year, lineages_at_founding,
created_alleles_surviving, mutational_diffs, created_diffs
```

plus the full lineage-count trajectory as a separate series, and — when W7 has run — the ARGweaver
inferred TMR-*K*-A joined on `locus` for the truth-vs-inference plot.

The headline figure of the study is **inferred TMR4A vs true TMR4A**, with created-allele divergence
as the swept variable.

---

## 5. Scale, compute, and the escape hatch

- The expensive case is the *mainstream comparison arm*, where the model under test posits a large
  ancient population — order 10⁴ diploids over order 10⁴ generations, forward-time. Note this is the
  size the mainstream model *assumes*, not an established fact about human history (§1B); the
  recent-origin arm is far cheaper, because a young structured population needs neither the depth nor
  the census size. The `child_founders/Q` rescaling already in
  [demography.go:275-285](pkg/demography/demography.go#L275-L285) applies: TMR-*K*-A measured in
  generations rescales the same way, so a Q=10 run may suffice. Worth confirming against your actual
  throughput before committing to a full-scale run.
- **Genome memory is the binding constraint at high X, and there is a way out.** Chromosomes are
  dense per-individual bitfields, so at 10⁶ bits an individual costs ~250 KB and a population of
  10⁴ costs ~2.5 GB resident — before the pedigree. But the TMR-*K*-A walker never reads the
  bitfield: it needs only the pedigree (W1) and strand provenance (W2). Sequences are needed *only*
  for the ARGweaver export. So with W2 in breakpoint mode the run can carry no per-individual
  genomes at all, and mutations can be dropped onto the recorded genealogy afterwards to synthesise
  the sample's sequences — the msprime/tskit pattern. This removes the memory ceiling on X and is
  the strongest argument for prioritising breakpoint mode over focal-loci mode once the walker is
  validated.
- W1 pedigree memory grows with everyone who left surviving descendants. Refcount pruning bounds it,
  but measure it early — this is the item most likely to force a design change.
- **Determinism is the escape hatch.** Per-run deterministic RNG means re-running with a different
  focal-locus set reproduces a byte-identical genealogy. Sweep loci in passes instead of holding a
  full ARG in memory. W8 must actually test this, or the whole strategy rests on an assumption.

---

## 6. Risks

1. **Rate anchors (W7)** — reduced from a scaling problem to a documentation problem now that θ/ρ is
   known to be scale-invariant (§1A, W7). What must be defended in writing before any headline run:
   the bp-per-bit convention, mutations per individual per generation, and crossovers per meiosis.
   The related *substantive* point is not a risk to the build but its subject matter: at real human
   rates each non-recombining block carries ~1 segregating site, so ARGweaver's per-block estimates
   are prior-dominated. Quantifying that is the study.
2. **W4 has no precedent in the codebase** and touches meiosis, the most correctness-critical path
   in the engine. Gate it hard; keep the created-gamete path structurally separate.
3. **Pedigree memory (W1)** at mainstream scale.
4. **Recurrent mutation** — positions collide under `RandIntn(genome_bits)`; identity must key on
   mutation ID and founder label throughout, never on coordinate.
5. **Interpretive risk.** Reporting a censored walk as a date would reproduce exactly the error the
   module exists to expose. The three outcome classes in W5 must survive all the way into the
   plots and the prose.

---

## 7. Suggested order

0. **W9 (constant-*N* sweep)** — needs no DRIFT code and tests the load-bearing claim directly: does
   inferred TMR4A track the assumed Ne at ≈*N* generations? Do this first. If it holds, every later
   result lands on prepared ground; if it doesn't, the whole critique needs rethinking before any
   engine work is committed.
1. **W1 + W2 (focal mode) + W5** — ground-truth TMR-*K*-A on ordinary diploid runs. Self-contained,
   publishable on its own, and validates the machinery before any created-allele modelling.
2. **W3 + W4** — created alleles and germline heterogeneity; the distinctively creationist model.
3. **W6 + W7 + W9 second axis** — the ARGweaver calibration, including structured-population runs
   scored against DRIFT's true §6d Ne trajectory.
4. **W2 breakpoint mode** — promote to the tskit substrate TODO §7 wants, unlocking genome-wide
   TMR-*K*-A instead of focal loci. **Promote this earlier if the memory measurement in §5 bites**:
   breakpoint mode plus after-the-fact mutation placement removes per-individual genomes from the
   run entirely, which is what makes a high-resolution single-arm study affordable.

Steps 1 and 3 are the calibration paper. Step 2 is the one that tests whether the metric means what
it is claimed to mean.

---

## 8. Implementation record — step 1 (W1 + W2 focal + W5 + W8)

Landed on `rob-TMR4A-1`. Every flag defaults off; the `-validate` gate is unchanged
(**D=−0.6611 / π-W=0.809 / Ne-N=0.313 / SFS-χ² 4.767**) and no path draws RNG unless opted into.

### What was built

| Item | Files | Notes |
|---|---|---|
| W1 pedigree | [pkg/core/pedigree.go](pkg/core/pedigree.go), `core.Pop.Pedigree`, [birth.go](pkg/simulation/birth.go), [death.go](pkg/simulation/death.go), [initializepop.go](pkg/config/initializepop.go) | `PedNode{Dad,Mom,BirthYear,Sex,Refs,Alive}`; refcount = alive + retained children; release cascades to parents, iteratively (deep lineages would blow a recursive stack). Founder nodes seeded once after setup dispatch, so every setup module is covered. |
| W2 provenance | [pkg/core/arg.go](pkg/core/arg.go), `core.Pop.ARGDB`, `core.Model.ARGLoci` | Focal-loci mode only. One flat `[]uint64` per birth: gamete 0 then gamete 1, one bit per focal locus. Pruned by the W1 cascade, so it cannot outlive the pedigree that indexes it. |
| W5 walker | [pkg/analysis/tmrka.go](pkg/analysis/tmrka.go), [tmrka_output.go](pkg/analysis/tmrka_output.go) | Max-heap over (individual, strand) lineages, popped newest-first. Emits the full trajectory; TMR-*K*-A for any *K* is `TMRKAForK`. |
| W8 tests | `pkg/core/{pedigree,arg}_test.go`, `pkg/simulation/{pedigree,arg,tmrka_sweep}_test.go`, `pkg/analysis/tmrka_test.go` | Hand-built pedigrees with known coalescence times; byte-identity; locus-sweep reproducibility. |

### Parameters

`track_pedigree` (0) · `track_arg` (0) · `arg_loci` ("") · `track_tmrka` (0) · `tmrka_k` (4) ·
`tmrka_sample_size` (0). `arg_loci` grammar: comma/space-separated `P`, `A-B`, or `A-B:S`.
`track_arg` requires `track_pedigree` — drift.go warns and skips capture otherwise.

### Decisions worth knowing before extending this

- **Mask inversion lives in exactly one function.** `meiosis` computes
  `(mask & parent0) | (^mask & parent1)`, so a SET mask bit means parental copy **0**. `ARGCapture`
  is the only place that inverts it, and `TestARGRecordedStrandMatchesMeiosisResult` checks the
  recorded strand against the allele meiosis actually delivered, under both the legacy and the §6a
  map recombination paths. Do not re-derive the convention at a read site.
- **Coalescences are dated to the ancestor's birth year**, and the event list is **sorted by year
  before** being turned into a trajectory. It is not discovered in time order: the walk pops by the
  descendant's birth year, and with DRIFT's long founder lifespans a late birth to an ancient parent
  yields an older event than one found later. An unsorted trajectory attaches the wrong counts to
  the wrong years.
- **The id tie-break in the walk heap is load-bearing.** Correctness rests on never arriving at a
  node already stepped past. That holds because every arrival is strictly older than the node just
  popped; descending-id tie-breaking preserves it even if a model ever allows a parent and child to
  share a birth year, since DRIFT allocates ids monotonically.
- **`FinalLineages` == `FounderStrandsSurviving` by construction** — surviving lineages sit on
  distinct founder strands. They are computed by two different routes and cross-checked, because
  this is the headline falsifiable number (§0) and a silent off-by-one in it would be invisible.
- **`created_alleles_surviving`, `mutational_diffs`, `created_diffs` are emitted as `-1`.** They are
  identity-by-*state* quantities and need W3 labels. Reporting the by-descent count in their place
  would collapse exactly the distinction the study turns on.
- **Censored walks write EMPTY date cells**, not zeros or sentinels — §6 risk 5. An empty cell
  cannot be averaged by accident; the `outcome` column says why in words.

### Smoke model

`users/smoke/models/TMR4ATest` — 200 diploids, 500 years, 13 focal loci at `0-3000:250`.
Representative output (rng_seed 4242): 13 loci → **5 coalesced, 8 censored at founding, 0 with no
coalescence**; mean TMR-4-A depth 297.8 y (11.9 generations); founder haplotypes surviving per
locus min 3 / mean 5.85 / max 11. Note the shape of that result: of 346 sampled lineages, all but
3–11 coalesce, but at most loci the survivors are *more* than four distinct founder haplotypes, so
there is no TMR-4-A to report — which is the censoring the module exists to make visible.

### Not yet done, in the §7 order

- **W9** (constant-*N* ARGweaver sweep) — needs no DRIFT code, but does need ARGweaver runs; it is
  external tooling and run management, not an engine change. Still the highest-value next step.
- **W3 + W4** — founder allele identity labels and created germline heterogeneity. W5 already has
  the fields reserved and the walker already reports the founder-strand count they refine.
- **W6 + W7** — merged VCF export and the ARGweaver bridge.
- **W2 breakpoint mode** — promote when the §5 memory measurement bites.
