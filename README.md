# reality

A Standard Model particle physics simulator in Rust, with a toolkit for proposing changes to the Standard Model and testing them against data without fitting.

It has two halves:

- **Standard Model calculations**: particles, tree-level cross sections from explicit helicity amplitudes, decay tables, Monte Carlo event generation and running couplings.
- **Model amendments**: describe new particles in a small TOML file, check the result for internal consistency, and derive predictions it was never fitted to: the strong coupling α_s(M_Z) from grand unification, and the proton lifetime. You can scan the new particles' masses, or search hundreds of thousands of candidate models automatically.

## Build

```sh
cargo build --release
cargo test
```

The examples below use `cargo run --release --`; the built binary is `target/release/reality`.

## Standard Model commands

```sh
reality particles [quarks|leptons|bosons|all]
reality cross-section ee_mumu --sqrt-s 91.19       # also ee_mumu_ew, ee_qqbar, ee_ww
reality decay z                                    # also w, higgs, top
reality generate ee_mumu_ew --sqrt-s 91.19 --n-events 100 --format table
reality coupling alpha_s --scale 1000              # also alpha_em, sin2w, gf, gw
```

## Unification and model amendments

### `reality unification`

One-loop running of the three gauge couplings, in GUT normalization (α₁ = 5/3 α_em / cos²θ_W), with an optional supersymmetric threshold.

```sh
reality unification                    # the SM couplings miss each other
reality unification --susy-scale 1000  # the MSSM couplings meet near 1.5e16 GeV
```

### `reality amend <model.toml>`

Evaluates a model amendment and prints a report:

1. **Consistency gates**, which need no data: gauge anomaly cancellation (U(1)³, gravity × U(1), SU(3)² × U(1), SU(2)² × U(1), SU(3)³), the Witten SU(2) anomaly, integer charges for colour singlets, and perturbativity up to 10¹⁹ GeV.
2. **Held-out prediction.** Two couplings are inputs. The engine requires all three to meet, and predicts the third at M_Z. Running is two-loop with the top Yukawa by default (`--loops 1` for the analytic one-loop result). The pull includes a theory uncertainty: 2% at two loops, 5% at one.
3. **Proton lifetime** for p → e⁺π⁰ via the X/Y gauge bosons, compared with the Super-K limit and the Hyper-K reach. The report gives a kill condition where the model is testable.
4. **Ledger**: the parameters set by hand, the parameters fitted (always zero), and the predictions made.

```sh
reality amend models/sm.toml              # alpha_s = 0.074, pull +18 sigma
reality amend models/mssm.toml            # alpha_s = 0.121, pull -1.0 sigma
reality amend models/mssm_su5.toml        # adds a GUT-threshold band
reality amend models/lone_higgsino.toml   # rejected by the anomaly gates
```

### `reality scan <model.toml>`

Moves the masses of the added fields across a log grid. The lightest one sweeps from `--from` to `--to`, and the others keep their ratio to it. Reports the best fit and the 1σ and 2σ mass ranges. `--format csv` exports the data for plotting.

```sh
reality scan models/mssm.toml --from 100 --to 1e7 --points 31
```

### `reality search`

Brute-force search for non-supersymmetric additions to the SM that unify. It enumerates every combination of up to `--max-fields` multiplets from a 20-entry menu, at each mass on a grid. The menu has vectorlike fermions, Weyl fermions in real representations, and scalars, so every combination is anomaly-free. Candidates are screened at one loop and refined at two. A model is accepted if:

- |pull| ≤ 1 on α_s
- the proton outlives the Super-K limit
- M_GUT ≤ 10¹⁸ GeV
- all couplings stay perturbative

Results are ranked by parsimony: multiplets, then mass scales, then new states, then |pull|. `--write DIR` saves the winners as model files.

```sh
reality search --max-fields 3 --top 6
reality search --max-fields 3 --masses 1e3,3.16e3,1e4,3.16e4,1e5 --write found/
```

The search runs in parallel. On a 10-core machine, 3 fields at 7 masses (477k candidates) takes about 8 s, and 4 fields (17M candidates) about 4 minutes.

## Model file format

```toml
name = "Vectorlike quark doublet at 2 TeV"
description = "Optional one-line description."
base = "sm"                       # "sm" (default) or "none"
inputs = ["alpha_1", "alpha_2"]   # two of alpha_1, alpha_2, alpha_3; the third is predicted
top_yukawa = "auto"               # "auto" (default) or "none"

[[field]]                         # a particle added to the running
name = "vectorlike Q"
kind = "dirac"                    # weyl, dirac, complex_scalar, real_scalar
su3 = "3"                         # 1, 3, 3bar, 6, 6bar, 8
su2 = 2                           # dimension of the SU(2) representation
y = "1/6"                         # hypercharge, Q = T3 + Y
copies = 1
mass = 2000                       # GeV; active above this scale

[[heavy]]                         # a field at the unification scale (threshold matching)
name = "colour-triplet Higgs"
kind = "chiral"                   # chiral, vector, vector_superfield
su3 = "3"
y = "-1/3"
copies = 2
ratio = [0.1, 10]                 # mass / M_X: fixed (1.0) or a range giving a band
group = "triplets"                # optional: fields in a group share one ratio
```

Notes:

- `weyl` is one left-handed Weyl fermion and can make a model anomalous. `dirac` is a vectorlike pair and is always anomaly-free.
- Two-loop coefficients are derived from the field content. Energy ranges where the content is supersymmetric (one gaugino per gauge group, with every fermion paired with a scalar) are detected automatically, and use the superfield formula and the MSSM top Yukawa.
- `[[heavy]]` fields don't run. Each shifts 1/αᵢ at M_X by bᵢ/(2π) · ln(ratio). Ranges are evaluated at every corner, which gives a band on α_s and M_X.

### Included models

| File | What it is |
|---|---|
| `models/sm.toml` | The Standard Model alone |
| `models/mssm.toml` | MSSM with superpartners at 1 TeV |
| `models/split_susy.toml` | Split SUSY: gauginos and higgsinos at 1 TeV, scalars at 10⁹ GeV |
| `models/mssm_su5.toml` | MSSM plus the GUT-scale fields of minimal SUSY SU(5) |
| `models/vectorlike_quark.toml` | One vectorlike quark doublet: consistent, but predicts α_s badly |
| `models/lone_higgsino.toml` | A deliberately inconsistent model |
| `models/found/*.toml` | The best results from `reality search` |

## Validation

The test suite checks the engine against published results:

- the one- and two-loop SM and MSSM beta-function coefficients, derived from field content
- the Hisano–Murayama–Yanagida relation for SU(5) threshold corrections (exact at one loop)
- the analytic near-threshold cross section for e⁺e⁻ → W⁺W⁻
- the known α_s predictions of the SM (about 0.07) and the MSSM (about 0.13 at two loops with superpartners at M_Z)

## Limitations

- **Sharp thresholds.** Each new particle switches on at its mass, with no matching corrections at superpartner or vectorlike-fermion thresholds.
- **Uncertain proton lifetime.** It's a dimension-6 estimate, uncertain by about a factor of 10. Dimension-5 decays in SUSY models (p → K⁺ν̄) aren't included.
- **Top Yukawa only.** Two-loop running includes only the top Yukawa, with the SM value at M_Z (tan β is ignored).
- **The search menu limits results.** It finds the best models among its 20 multiplets, not every possible model.
