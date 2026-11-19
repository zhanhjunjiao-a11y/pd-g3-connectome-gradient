# Supplementary Software 1

## Scope

Analysis and figure-generation source corresponding to the reported numerical results for:
**Altered whole-brain connectome gradients and associated functional connectivity changes in Parkinson’s disease**

The source corresponds to the 334 PD / 141 HC disease cohort, the same 235 clinical and imaging pairs, and the 359 PD / 79 HC external cohort. A historical 236-record input pool is retained by the local workflow; `paired_analysis_policy.json` selects the frozen 235 active cases. That policy and the individual data are controlled inputs, not supplied in this archive.

This is a **source-code archive**, not a self-contained individual-data rerun bundle. All 26 aggregate result tables accompany the manuscript and are provided in [Supplementary_Data](Supplementary_Data/). In this repository S03 is a [gzip-compressed CSV](Supplementary_Data/Supplementary_Data_S03_wholebrain_edge_FC.csv.gz). Restore its original CSV filename using the command in the [data notes](Supplementary_Data/README.md) before using tools that expect the uncompressed filename. The data dictionary uses the restored CSV name. No participant records, subject identifiers, raw images, clinical registers, individual plotting observations, or internal audit logs are redistributed here. File names containing “anonymous” in the plotting code still refer to controlled individual observations.

## Source map

- `modules/core/run.py`: regional gradient, complete-system, paired and clinical models, including the current eight-measure clinical families.
- `modules/connectivity/run.py`: FC-to-gradient calculations, edgewise and region-to-target FC, external profiles and target decompositions.
- `modules/external/`: PPMI endpoints, adjusted comparisons, restricted wild-bootstrap effect-map tests, cohort characteristics and Table 1.
- `modules/figures/`: four main and three supplementary figures.
- `requirements.txt`: package versions from the local analysis environment.

## Controlled inputs and execution

Running an analysis requires the authorised package-relative `inputs_private` files, frozen participant policy, and corresponding `verification` reference files described in each source module. These must be requested under the article's Data availability terms; the archive does not grant access to individual PPMI data. Atlas/mesh resources also require their respective access and licence conditions.

With authorised inputs in place and a compatible Python environment, module entry points include:

```text
python modules/core/run.py --output-dir runs/core
python modules/connectivity/run.py --output-dir runs/connectivity
python modules/external/run.py --output-dir runs/external
python modules/figures/run.py --output-dir runs/figures
```

Inspect each entry point's arguments before running. Figure 2 and Figure 3A require individual plotting data; their point-level data are not part of the aggregate release. Analysis-generated endpoints, design records and diagnostic outputs may also be individual-level and should remain controlled.

Raw MRI preprocessing is not rerun by these modules. The manuscript describes the recovered preprocessing evidence. This packaging step checked syntax, source hashes and the absence of embedded participant records; it did not rerun the models.

## Changes made for this archive

Eighteen of the 19 analysis/plot source files are copied unchanged. In `modules/external/endpoints.py`, only the `report_role` assignment is updated: `scanner_adjusted_motion` is primary and `age_sex_motion` is a sensitivity model. All model definitions, scanner categories, inputs, estimates and BH families remain unchanged. The portable `guarded_runner.py` retains the previously adapted removal of institution-specific historical absolute paths. It restricts executed scripts to the package and disables network socket connections, but does not claim to be a complete file-access sandbox. No scientific model or statistic was changed, and no models were rerun for this reporting update.

For external PPMI reporting, the primary model adjusts for age, sex, meanFD and categorical scanner model. The age/sex/meanFD base model and shared clinical SITE model remain sensitivity analyses. The ambiguous-age exclusion, screening-failure exclusion and original-age models are base-model sensitivities without scanner covariates; they do not test robustness of the scanner-adjusted primary model. S14, S16 and S23 carry updated reporting roles; the frozen S24/S25 model codes identify the same primary/sensitivity mapping. The single fixed 19-region score has nominal uncorrected P = 0.0388 in the scanner-adjusted primary model and P = 0.0639 in the base-model sensitivity.

`MANIFEST.json` provides source hashes. The authors’ analysis code is available under the [MIT License](LICENSE.txt). Third-party atlas and mesh resources retain their respective licences. Participant-level data remain subject to institutional and PPMI access conditions.

## Final figure presentation

The four main figures supplied with the 19 November manuscript incorporate presentation edits completed on 18 November 2026. These edits use lowercase panel labels, expand the deep cerebellar nuclei label in Figure 1, and remove the construction schematic below Figure 2. Figure 2f retains the connectivity-profile similarity results. The analysis values and plotted observations are preserved.

The archived plotting scripts reproduce the verified numerical plots and their earlier presentation. The final main-figure PNGs are supplied in this repository's [Figures directory](Figures/). This distinction concerns presentation, including labels and the supplementary construction schematic, and is separate from statistical reproducibility.

The complete controlled local package passed numerical verification on 3 November 2026 following restriction to the confirmed 235 paired participants. That check covered all 26 result datasets, Table 1, four main and three supplementary figures, and 38 figure-source comparisons. The current public archive has undergone ZIP integrity, source-hash, syntax and input-boundary checks. Full model fitting requires the authorised controlled inputs described above.
