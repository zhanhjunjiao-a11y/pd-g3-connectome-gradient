# PD cortico-cerebellar connectome-gradient study — analysis code

Code to reproduce the results, figures, and tables of the study of third-connectome-gradient (G3)
compression in Parkinson's disease. Data are archived separately on Zenodo (DOI: **[add at
submission]**); subject-level imaging data are available on request under a data-use agreement
(see the paper's Data-availability statement).

## Repository layout
```
03_analysis/          the 9-stage pipeline
  a_gradient/          connectome gradients (diffusion-map embedding + Procrustes)
  b_cpm_prediction/    ridge / leave-one-out clinical prediction (exploratory; null)
  c_pet_neurotransmitter/  spatial correlation vs 20 neuromaps PET maps
  d_gene_expression/   AHBA/abagen expression → ROI matrices
  e_gsea_pathway/      genome-wide pre-ranked GSEA (KEGG/GO/Reactome/WikiPathways)
  f_cell_type/         cortical cell-type marker enrichment
  g_yeo7_network/      Yeo-7 network / compartment enrichment
  h_spin_null/         variogram (BrainSMASH) spatial-null construction
  i_reliability/       test–retest reliability
04_figures/           main + supplementary figure scripts (shared style: _nc_style.py)
results/              per-analysis robustness/sensitivity runs (run_*.py) incl.
  reviewer_tissue_confound/   tissue-class (cortex/subcortex/cerebellum) confound sensitivity
  item_gmv_structural_covariate/  gray-matter-volume-adjusted replication
```

## Environment
Python 3.11. Install: `pip install -r requirements.txt`.
Preprocessing (upstream, not in this repo): fMRIPrep 23.2.0a3, FSL 6.0.5.

`requirements.txt` is pinned to the **primary-analysis versions** reported in the paper
(brainspace 0.1.10, brainsmash 0.11, gseapy 1.0.6). The reviewer-response tissue-class
sensitivity re-runs used gseapy 1.1.13 (direction/significance unchanged); see the note in
`requirements.txt`.

## Reproducing
Paths in scripts are absolute (`E:\预测模型\...`); set a project root and adjust, or place the
Zenodo `zenodo_data/` archive alongside. Group-level inputs (t-maps, spatial-null surrogates,
atlas assets, expression tables) are on Zenodo; subject-level inputs are controlled-access.
Each `run_*.py` writes a `_RESULTS_SUMMARY.md` next to its outputs. See Supplementary §S20 for
version pins and step order.
