# PD whole-brain connectome gradients

Analysis code, aggregate results and final main figures for **Altered whole-brain connectome gradients and associated functional connectivity changes in Parkinson’s disease**.

Version 2026.11.19 accompanies the current manuscript. The discovery sample includes 334 patients with Parkinson’s disease (PD) and 141 healthy controls. The acute levodopa analyses use the same 235 patients for clinical and imaging comparisons. The independent Parkinson’s Progression Markers Initiative (PPMI) sample includes 359 patients with PD and 79 controls.

## Repository contents

The [consolidated Excel workbook](Supplementary_Data/Supplementary_Data_S01-S26.xlsx) provides all 26 result tables as sheets S01–S26, together with a contents index, field definitions and analysis notes. It is the single spreadsheet attachment for journal submission. The original CSV files below remain available for analysis use.

- [Supplementary_Data/](Supplementary_Data/) contains all 26 Supplementary Data S01–S26 result tables, with S03 supplied as a gzip-compressed CSV, together with their [data dictionary](Supplementary_Data/Data_dictionary.csv) and [analysis notes](Supplementary_Data/README.md).
- [modules/core/](modules/core/) contains regional-gradient, system, paired and clinical analyses.
- [modules/connectivity/](modules/connectivity/) contains FC and gradient calculations, regional and system connectivity analyses, and connectivity-profile analyses.
- [modules/external/](modules/external/) contains PPMI analyses and descriptive tables.
- [modules/figures/](modules/figures/) contains the figure-generation source.
- [Figures/](Figures/) contains the four final main-figure PNGs.
- [CODE_README.md](CODE_README.md) describes controlled inputs, execution and the source archive.
- [requirements.txt](requirements.txt) records dependencies. [MANIFEST.json](MANIFEST.json) records analysis and plotting source hashes.

## Opening Data S03

[Data S03](Supplementary_Data/Supplementary_Data_S03_wholebrain_edge_FC.csv.gz) contains the complete 39,903-row whole-brain connection table in gzip format. From the repository root, restore its original CSV filename with:

```text
python -m gzip -d Supplementary_Data/Supplementary_Data_S03_wholebrain_edge_FC.csv.gz
```

The data dictionary refers to the decompressed filename. The SHA-256 of the restored CSV is `85455beb0b31759546140177d09543af1bee939b8e8dc42b3f50c7398f95d8a0`. The manuscript’s consolidated Excel attachment contains the same table on sheet S03.

## Analysis definitions

The complete somatomotor network contains 30 regions and the cerebellum contains 31. Their external connectivity profiles use the same 222 remaining regions. The fixed 19-region G3 score uses the discovery cohort’s 10 higher-G3 regions minus its nine lower-G3 regions.

The primary PPMI model adjusts for age, sex, mean framewise displacement and categorical scanner model (`scanner_adjusted_motion`). The age/sex/motion model and the shared clinical SITE model are sensitivity analyses. The ambiguous-age exclusion, screening-failure exclusion and original-age checks use the base covariates. The data notes specify the sample sizes and multiple-testing families for each analysis.

## Running the code

Use a Python environment compatible with [requirements.txt](requirements.txt). Analysis execution requires the authorised `inputs_private/` files, frozen participant-selection policy and reference files in `verification/`, as documented in the source modules and [CODE_README.md](CODE_README.md). Obtain access under the manuscript’s Data availability terms through the corresponding authors. PPMI data access remains subject to PPMI’s access conditions.

With the required environment and authorised inputs in place, the module entry points are:

```text
python modules/core/run.py --output-dir runs/core
python modules/connectivity/run.py --output-dir runs/connectivity
python modules/external/run.py --output-dir runs/external
python modules/figures/run.py --output-dir runs/figures
```

The repository supplies aggregate results and source code. Participant-level clinical records, FC or gradient arrays, and individual plotting observations are controlled-access inputs. Raw MRI preprocessing is described in the manuscript.

## Verification and figures

The complete controlled local package passed numerical verification on 3 November 2026 after the paired sample was fixed at 235 patients. Verification covered all 26 result datasets, Table 1, four main figures, three supplementary figures and 38 figure-source comparisons. The 19 November public package was checked for source syntax, source hashes, CSV completeness, data-dictionary coverage and the public release boundary.

The [final figures](Figures/) incorporate presentation edits completed on 18 November, including lowercase panel labels, the expanded deep cerebellar nuclei label and removal of the construction schematic below Figure 2. The archived scripts reproduce the verified statistical plots in their earlier presentation. Figure 2f retains the connectivity-profile similarity results.

## Citation

[CITATION.cff](CITATION.cff) provides the authors, study title and version for citing this repository.

## Code licence

The authors’ analysis code uses the [MIT License](LICENSE.txt). Third-party resources and controlled participant-level data retain their respective licence and access conditions.
