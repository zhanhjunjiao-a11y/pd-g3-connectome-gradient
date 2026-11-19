# Supplementary statistical data

**Altered whole-brain connectome gradients and associated functional connectivity changes in Parkinson’s disease**

This folder contains complete aggregate results supporting the article and Supplementary Information, including tests that did not reach the significance threshold. Data S01–S26 comprise 26 comma-separated-value (CSV) tables. S03 is supplied as a gzip-compressed CSV, and the other 25 tables are supplied as plain CSV files. [Data_dictionary.csv](Data_dictionary.csv) defines each field by file. The files contain aggregate statistical, regional or connection-level results, rather than participant-level records or participant identifiers. Region identifiers refer to atlas parcels, not participants.

## Opening S03

The file `Supplementary_Data_S03_wholebrain_edge_FC.csv.gz` preserves the complete original CSV. From the repository root, run:

```text
python -m gzip -d Supplementary_Data/Supplementary_Data_S03_wholebrain_edge_FC.csv.gz
```

This restores `Supplementary_Data_S03_wholebrain_edge_FC.csv`, the filename used in the data dictionary. Its SHA-256 is `85455beb0b31759546140177d09543af1bee939b8e8dc42b3f50c7398f95d8a0`. Row counts, field definitions and numerical results apply to this decompressed table.

The consolidated [Excel workbook](Supplementary_Data_S01-S26.xlsx) contains sheets S01–S26, a contents index, the data dictionary and analysis notes. Original CSV files remain available for analysis use (S03 is gzip-compressed). The workbook preserves source row and column order; numeric columns containing more than 15 significant digits are stored as text to avoid Excel precision loss.

## File index

Row counts exclude the header. Data-file numbers and Supplementary Table numbers are separate numbering systems.

| Data | Contents | Rows |
|---|---|---:|
| [S01](Supplementary_Data_S01_atlas_and_networks.csv) | Atlas regions and system assignments | 283 |
| [S02](Supplementary_Data_S02_disease_gradients.csv) | Regional disease-related gradients | 849 |
| [S03](Supplementary_Data_S03_wholebrain_edge_FC.csv.gz) | Whole-brain edgewise functional connectivity | 39,903 |
| [S04](Supplementary_Data_S04_G3_selected19_edge_display.csv) | Connection-display results for the 19 selected G3 regions | 5,358 |
| [S05](Supplementary_Data_S05_ROI_to_network_mean_FC.csv) | Mean FC from each region to nine target categories | 2,547 |
| [S06](Supplementary_Data_S06_complete_system_disease_family.csv) | Complete-system disease effects | 6 |
| [S07](Supplementary_Data_S07_paired_regional_gradients.csv) | Paired regional gradient effects | 849 |
| [S08](Supplementary_Data_S08_complete_system_paired_family.csv) | Complete-system paired effects | 24 |
| [S09](Supplementary_Data_S09_motor_severity.csv) | OFF-state motor severity associations | 8 |
| [S10](Supplementary_Data_S10_motor_improvement.csv) | Acute motor improvement associations | 8 |
| [S11](Supplementary_Data_S11_edge_and_ROI_mean_sensitivity.csv) | Edgewise and regional-FC sensitivity analyses | 15,468 |
| [S12](Supplementary_Data_S12_complete_system_auxiliary_family.csv) | Auxiliary complete-system disease measures | 7 |
| [S13](Supplementary_Data_S13_G3_motion_sensitivity.csv) | G3 motion sensitivity analyses | 1,132 |
| [S14](Supplementary_Data_S14_PPMI_complete_system_external.csv) | PPMI complete-system FC | 24 |
| [S15](Supplementary_Data_S15_model_diagnostics.csv) | Model diagnostics | 857 |
| [S16](Supplementary_Data_S16_PPMI_fixed19_score_external.csv) | PPMI fixed 19-region G3 score | 6 |
| [S17](Supplementary_Data_S17_profile_target_networks.csv) | Profile decomposition by target network | 35 |
| [S18](Supplementary_Data_S18_profile_target_regions.csv) | Profile decomposition by target region | 222 |
| [S19](Supplementary_Data_S19_profile_geometry_correspondence.csv) | Associations between profile distance and gradient geometry | 6 |
| [S20](Supplementary_Data_S20_profile_distance_mean_dispersion.csv) | Descriptive mean–dispersion decomposition | 7 |
| [S21](Supplementary_Data_S21_complete_system_geometry_paired.csv) | Paired complete-system gradient geometry | 3 |
| [S22](Supplementary_Data_S22_PPMI_characteristics.csv) | Aggregate PPMI characteristics and acquisition metadata | 92 |
| [S23](Supplementary_Data_S23_PPMI_system_geometry.csv) | PPMI complete-system gradient geometry | 18 |
| [S24](Supplementary_Data_S24_PPMI_effect_map_agreement.csv) | PPMI–discovery gradient disease-effect map correspondence | 9 |
| [S25](Supplementary_Data_S25_PPMI_regional_gradients.csv) | PPMI whole-brain regional gradient tests | 2,547 |
| [S26](Supplementary_Data_S26_complete_system_motion_sensitivity.csv) | Complete-system disease comparisons: head-motion sensitivities | 39 |

## Definitions and interpretation

PD denotes Parkinson’s disease; HC, healthy controls; FC, functional connectivity; SMN, somatomotor network; CB, cerebellum; DMN, default mode network; meanFD, mean framewise displacement; and ΔmeanFD, ON−OFF meanFD. HC3 denotes heteroscedasticity-consistent robust standard errors. BH denotes Benjamini–Hochberg multiple-comparison correction. The number following BH identifies the size of a correction family, not the number of significant results.

The discovery cohort comprised 334 PD and 141 HC. All paired clinical, imaging and clinical-change analyses included the same 235 patients with complete paired data, all of whom belonged to the baseline PD sample. Disease effects are PD−HC, imaging changes are ON−OFF, and motor improvement is OFF−ON. For motion-adjusted paired models, the reported estimate is the intercept at zero motion change in `ON−OFF change ~ 1 + uncentred ΔmeanFD`, rather than the raw mean change.

All hypothesis tests are two-sided. Confidence intervals (CIs) are pointwise, rather than simultaneous. The significance threshold is BH q < 0.05 within the stated family. Correction families include all tested outcomes, including those not displayed in the main figures. Failure to reach the significance threshold does not establish equivalence. Selected-region FC follow-up analyses are exploratory and use the same underlying FC data as the gradients.

Mean FC is expressed in Fisher-z units. External-profile angular similarity is dimensionless and equals `1 − arccos(cosine similarity)/π`. Here, “external” refers to targets outside both systems in a pair. Target-distance contributions sum to `1 − cosine similarity`, not angular similarity. Gradient measures retain aligned coordinate units. Clinical standardised β is the unstandardised coefficient multiplied by `sample SD(predictor)/sample SD(outcome)`.

The primary PPMI model (`scanner_adjusted_motion`) includes 359 PD and 79 HC and adjusts for age, sex, meanFD and categorical scanner model using the frozen scanner categories. The base model (`age_sex_motion`), adjusting for age, sex and meanFD, is retained as a sensitivity analysis. Clinical SITE sensitivity (`site_shared_motion`) uses the 27 clinical sites containing both groups (304 PD, 78 HC), adjusting for age, sex, meanFD and clinical SITE; SITE codes do not establish physical MRI acquisition facilities. Other scalar-outcome base-model sensitivities exclude nine ambiguous-date participants, exclude one control with a later screening-failure label, or retain all original ages. These three checks adjust for age, sex and meanFD without scanner covariates; they do not assess robustness of the scanner-adjusted primary model. In file names, raw means before ComBat, not before Fisher transformation. All 438 FC inputs were verified against saved processed time series; no second Fisher transformation was applied. Upstream workflows and derivative metadata were inspected, and BOLD-to-time-series-to-FC reconstruction was checked in three PD and three HC examples. MRI-day medication status is not uniformly documented. Clinical assessment labels must not be used to infer MRI-day state.

## Analysis-specific notes

**S01 — Atlas and system assignments.** The 283 atlas regions and the fixed network labels used in the analyses are provided. SMN contains 30 regions and the cerebellum 31 regions. Coordinates are not included because coordinate-column labels in the source file were ambiguous.

**S02 — Regional disease-related gradients.** All 849 G1–G3 tests are included. Ordinary least-squares models used conventional standard errors, adjusting for age, sex and meanFD, with df = 470. BH283 correction was applied separately to each gradient. CIs were not available in the archived output.

**S03 — Whole-brain edgewise FC.** All 39,903 unique connections are included. HC3 models adjusted for age, sex and meanFD, with df = 470 and BH39903 correction.

**S04 — Connection display for the 19 selected G3 regions.** The 19 × 282 rows correspond to 5,187 unique connections. Connections between two selected regions appear twice. The q values are inherited from S03 and were not recalculated within this subset.

**S05 — Mean FC from regions to nine target categories.** All 2,547 region–target tests are included, excluding self-connections. Models used HC3 standard errors and adjusted for age, sex and meanFD, with BH2547 correction. The file includes the results used for the G2 and G3 follow-up displays.

**S06 — Complete-system disease effects.** All six tests in the disease family are provided, including the G3 mean difference and the within- minus between-system FC difference. Models used HC3 standard errors and adjusted for age, sex and meanFD, with BH6 correction.

**S07 — Paired regional gradients.** All 849 G1–G3 tests in 235 paired observations are included. The model was `ON−OFF change ~ 1 + uncentred ΔmeanFD`, using HC3 standard errors and df = 233. Separate columns retain the primary BH283 correction and sensitivity BH849 correction.

**S08 — Complete-system paired effects.** Two models each contain 12 tests: the motion-adjusted analysis and the unadjusted analysis of the same 235 paired observations. BH12 correction was applied separately to each model, not jointly across all 24 rows. G2 system-pair differences are the second system minus the first. Shared external targets exclude both systems in the pair: 207 targets for DMN–SMN, 206 for DMN–CB and 222 for SMN–CB.

**S09 — OFF-state motor severity.** Eight separate models, each containing one imaging predictor, were fitted in 334 PD, adjusting for age, sex and meanFD. The outcome was OFF-state MDS-UPDRS-III. Models used HC3 standard errors, df = 329 and BH8 correction.

**S10 — Acute motor improvement.** Eight separate models, each containing one imaging-change predictor, were fitted in 235 paired observations. The outcome was the OFF−ON MDS-UPDRS-III difference, and the predictor was ON−OFF imaging change. Covariates were age, sex, ΔmeanFD, OFF-state motor score and the corresponding OFF-state imaging measure. Models used HC3 standard errors, df = 228 and BH8 correction. The correction family is separate from S09.

**S11 — Edgewise and regional-FC sensitivity analyses.** Analyses additionally adjusted for whole-brain mean FC or excluded one control flagged during an earlier record-linkage review because of conflicting duplicate clinical entries. Independently traceable imaging, covariate and motion records supported retention in the primary analysis; exclusion was not selected using an influence-statistic threshold. The displayed edgewise results retain BH39903 correction, and regional means retain BH2547 correction, separately for each model. These analyses did not test complete-system profile similarity.

**S12 — Auxiliary complete-system disease measures.** Seven outcomes include 3D/10D geometry, centred-profile correlation, profile root-mean-square difference and whole-brain mean FC. Models used HC3 standard errors and adjusted for age, sex and meanFD, with BH7 correction.

**S13 — G3 motion sensitivity analyses.** Four G3 models are included: the reference model with linear meanFD, a model additionally including meanFD², restriction to participants with meanFD < 0.2 mm, and HC3 inference. Each model used separate BH283 correction. These checks concern regional effect maps only.

**S14 — PPMI system FC.** Four outcomes in six models (24 rows). Primary (`scanner_adjusted_motion`): age, sex, meanFD and categorical scanner model (n = 438); base-model sensitivity (`age_sex_motion`): age, sex and meanFD (n = 438); shared-SITE sensitivity (`site_shared_motion`): age, sex, meanFD and clinical SITE (n = 382). Other base-model sensitivities exclude nine ambiguous-date cases (n = 429), exclude a control with a later screening-failure label (n = 437), or use original ages (n = 438). These three models do not include scanner covariates and do not assess robustness of the scanner-adjusted primary model. HC3 standard errors, residual-df t inference and BH4 correction separately within each model. All results, including nonsignificant outcomes, are retained.

**S15 — Model diagnostics.** Diagnostics cover regional paired-gradient models and four FC clinical models; they are not a complete diagnostic set for clinical analyses of all eight imaging measures. Diagnostic P values are descriptive. Participant-level residuals and identifiers are not included.

**S16 — PPMI fixed 19-region G3 score.** Six models using the same covariate and sample rules as S14. Discovery regions, signs and gradient reference were fixed. The single composite outcome is reported with uncorrected P values. The scanner-adjusted primary model (`scanner_adjusted_motion`) is nominally significant (uncorrected P = 0.0388); the base-model sensitivity (`age_sex_motion`) is not (P = 0.0639). This single-score result is uncorrected for multiple comparisons and varies across models. The ambiguous-age, screening-failure and original-age checks are base-model sensitivities without scanner covariates.

**S17 — Target-network decomposition.** The 35 outcomes belong to four separate BH families: seven cosine-distance contributions, seven centred contributions, 14 raw FC outcomes and seven signed SMN−CB differences. Models used HC3 standard errors and adjusted for age, sex and meanFD, with df = 470.

**S18 — Target-region decomposition.** The 222 rows jointly contain 222 regional distance tests and 444 seed–target raw FC tests. BH222 and BH444 are separate correction families. All targets are retained.

**S19 — Profile distance and gradient geometry.** Six within-group associations examine cosine distance in relation to the signed G3 mean difference and 3D/10D geometry, separately in HC and PD. Models adjusted for age, sex and meanFD and used HC3 standard errors, with BH6 correction. Predictors and outcomes were scaled using within-group sample SDs.

**S20 — Descriptive mean–dispersion decomposition.** Seven target categories are included. The unadjusted mean individual distance equals one-half the squared group mean of the normalised profile difference plus one-half the between-participant variance of that difference, using variance denominator n. Terms are summed within each category. No tests or CIs are calculated, and the adjusted group effect is not decomposed. See the article’s Methods for the definitions.

**S21 — Paired complete-system geometry.** Three motion-adjusted ON−OFF intercept tests in 235 paired observations used HC3 standard errors, df = 233 and BH3 correction. Raw paired mean changes and t-based CIs are reported separately as descriptive results, using df = 234.

**S22 — PPMI characteristics.** Aggregate summaries, variable-specific denominators and documented age rules describe the 359 PD and 79 HC. All have derivative-based meanFD; 437 have a clinical SITE code. Twenty-four ages were corrected from MRI-source records; nine unresolved-date cases retained their original recorded ages. The other recorded ages remained unchanged. Clinical summaries retain the prior visit-and-scan-month matching; the 291 untreated/OFF examinations comprise 284 labelled untreated and seven labelled OFF, which does not determine MRI-day medication status. State-specific clinical subsets can overlap and must not be summed. Saved-series length and nominal TR × frame count are not evidence for original acquisition duration or initial-volume removal. Participant identifiers, dates and local source paths are omitted.

**S23 — PPMI system gradient geometry.** Three outcomes in six models (18 rows), using the same covariate/sample rules as S14: `scanner_adjusted_motion` is primary; `age_sex_motion` and `site_shared_motion` are sensitivity models. The remaining three checks are base-model sensitivities without scanner covariates. Outcomes are the signed SMN-minus-cerebellum G3 mean difference and 3D/10D centroid distances, using the complete 30-region SMN and 31-region cerebellum. Coordinates retain the original aligned scale. HC3 standard errors and residual-df t inference; BH3 correction separately for each model.

**S24 — PPMI–discovery gradient disease-effect map correspondence.** Three gradients in the scanner-adjusted primary model (`scanner_adjusted_motion`), base-model sensitivity (`age_sex_motion`) and shared-SITE sensitivity (`site_shared_motion`) (nine rows). All 283 regional group coefficients are used without significance selection. Discovery maps remain fixed. Two-sided P values use 9,999 restricted HC2 wild-bootstrap samples (seed 20261030) with BH3 within each model. The null is no regional group effects in PPMI conditional on the design and fixed discovery map. Covariate-only residuals are divided by sqrt(1 − h), and participant-level Rademacher multipliers are shared across all regions and gradients. The full group model is refitted and its Spearman correlation with the fixed discovery map calculated. P = (1 + number of absolute null correlations at least as extreme as observed)/10,000. This is not a spatial-rotation specificity test and does not propagate uncertainty in the discovery map. Confidence intervals from a separate unrestricted bootstrap are retained as internal diagnostics, not combined here with this test.

**S25 — PPMI regional gradient tests.** Three models × three gradients × 283 regions give 2,547 rows. The scanner-adjusted primary model (`scanner_adjusted_motion`) and base-model sensitivity (`age_sex_motion`) use 359 PD and 79 HC; the shared-SITE sensitivity (`site_shared_motion`) uses 304 PD and 78 HC. The primary model adjusts for age, sex, meanFD and categorical scanner model; the base model adjusts for age, sex and meanFD, and the shared-SITE model adds clinical SITE to those base covariates. PD−HC coefficients use HC3 standard errors, residual-df t inference and BH283 separately for each gradient and model. These families are distinct from the map-correspondence BH3 families in S24.

**S26 — Complete-system head-motion sensitivity analyses.** The reference model (334 PD, 141 HC), a low-motion model restricted to meanFD < 0.2 mm (200 PD, 105 HC), and an additional meanFD² model (334 PD, 141 HC) each include all 13 outcomes. Existing participant-level summaries for the complete 30-region SMN and 31-region cerebellum were reused without re-embedding or re-selecting regions. All models adjust for age, sex and linear meanFD; the quadratic model also includes uncentred meanFD². PD OFF−HC effects use HC3 standard errors and two-sided residual-df t inference (df = 470, 300 and 469, respectively). The original six-outcome disease family and seven-outcome auxiliary family retain separate BH correction within each model, including every family member even when only seven core outcomes are displayed in Table S15. The six-outcome family comprises the signed G3 mean difference, angular profile similarity, SMN-internal FC, cerebellar-internal FC, between-system FC and within-minus-between FC. The seven-outcome family comprises the two system G3 means, profile Pearson-z, profile RMS difference, 3D/10D centroid distances and whole-brain mean FC. The reference results reproduce S06/S12. These sensitivities concern disease-group comparisons only, not paired medication-state changes or clinical associations. They do not use frame censoring and cannot establish the absence of residual motion effects.

## Release notes for 19 November 2026

The 26 result CSV files and data dictionary are the frozen files used for the 15 November supplement update. Their numeric results, model identifiers and reporting roles are preserved. In archived metadata, the term “blinded clinical SITE” denotes anonymised clinical-site coding. The manuscript and explanatory text use “clinical SITE”.
