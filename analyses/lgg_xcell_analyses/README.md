
## Author: Komal S. Rathi, Adam Kraya
 
### Purpose

Preliminary analysis for immune profiles of low-grade gliomas. Total number of RNA-seq samples that match our criteria of `short_histology == "LGAT"`, `broad_histology == "Low-grade astrocytic tumor"` and `cancer_group == "low-grade glioma/astrocytoma"`: **536**.

 ### Data version

Latest data release: [Release 20230826](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230826_release)

Latest imaging clusters: `ImagingClusterAssignment_Aug2023.xlsx`

### 0. xCell-derived immune scores

`00-run_xcell.R`: Run `xCell` deconvolution on LGG samples of interest. 

#### Run script

```
Rscript 00-run_xcell.R --help

Options:
  --mat=MAT
    input matrix for immune profiling (.rds)

  --histology=HISTOLOGY
    histology file for all OpenPedCan samples (.tsv)

  --short_hist=SHORT_HIST
    short histology filter

  --broad_hist=BROAD_HIST
    broad histology filter

  --cg=CG
    cancer group filter

  --output_dir=OUTPUT_DIR
    output directory

  -h, --help
    Show this help message and exit

```

#### Output

Output files in results directory:

```
results/xcell_output
└── xcell_scores.rds
```

### 1. xCell-derived clustering

`01-xcell_clustering.R`: Perform clustering on xCell-derived scores using [ClusTarIDseq](https://github.com/d3b-center/ClusTarIDseq). 
 
#### Run script

```
Rscript 01-xcell_clustering.R --help

Options:
  --mat=MAT
    input matrix (.rds)

  --histology=HISTOLOGY
    histology file (.tsv)

  --output_dir=OUTPUT_DIR
    output directory for clustering analyses

  --xcell_output_dir=XCELL_OUTPUT_DIR
    output directory for xcell analyses

  --vtest_output_dir=VTEST_OUTPUT_DIR
    output directory for v-test analyses

  -h, --help
    Show this help message and exit
```

Following key parameters were used:

```
filter_expr = FALSE 
protein_coding_only = FALSE
feature_selection = "variance"
var_prop = 100
transformation_type = "none" 
max_k = 10
```

#### Output

- ClusTarIDseq output files in the results directory:

```
results/clustering
├── ccp_output
│ ├── ccp_optimal_clusters.tsv
│ ├── {hc, km, pam}_{binary, canberra, euclidean, maximum, pearson, spearman}_100.pdf
│ └── {hc, km, pam}_{binary, canberra, euclidean, maximum, pearson, spearman}_100.rds
├── dbscan_output
│ ├── KNNdistplot_output.pdf
│ ├── cluster_plot.pdf
│ └── dbscan_optimal_clusters.tsv
├── final_score
│ └── final_clustering_output.tsv # tsv file with ranks assigned to each evaluated method
└── intnmf_output
    ├── intnmf_consensus_plot.pdf
    ├── intnmf_fit.rds
    ├── intnmf_optimal_clusters.tsv
    ├── intnmf_optk.rds
    └── intnmf_silhouette_plot.pdf
```

- xCell-derived clusters are mapped back to xCell scores + clinical variables like `molecular_subtype` for downstream analyses and written to the following file:

```
results/xcell_output
└── xcell_score_cluster.tsv
```

- Additionally, a v-test is applied to each gene to identify its contribution to individual clusters. These v-test scores are written as follows:

```
results/vtest_analysis
└── vtest_scores_all.tsv
```

### 2. xCell-derived score heatmap

`02-xcell_heatmap.R`: For the top clustering method chosen in the above script, a heatmap is generated with annotations like `Cluster assignment` and `Molecular subtype`.

#### Run script

```
Rscript 02-xcell_heatmap.R --help

Options:
  --xcell_file=XCELL_FILE
    xcell score + cluster file (.tsv)

  --plots_dir=PLOTS_DIR
    output directory to save plots

  -h, --help
    Show this help message and exit
```

#### Output

Heatmap of xCell-derived scores annotated with clusters and molecular subtypes:

```
plots/xcell_cluster_analysis
└── xcell_heatmap.pdf
```

### 3. xCell-derived cluster analysis

`03-xcell_cluster_analysis.R`: This script generates stats and plots for the comparison between xCell-derived clusters and clinical variables like `TMB` and `molecular_subtype`.

#### Run script

```
Rscript 03-xcell_cluster_analysis.R --help

Options:
  --xcell_file=XCELL_FILE
    xcell score + cluster file (.tsv)

  --histology=HISTOLOGY
    histology file (.tsv)

  --output_dir=OUTPUT_DIR
    output directory to write files

  --plots_dir=PLOTS_DIR
    output directory to save plots

  -h, --help
    Show this help message and exit
```

#### Output

- Chi-square test of association between xCell-derived clusters and molecular subtypes:

```
results/xcell_cluster_analysis
└── xcell_clusters_vs_subtype_chisquare.txt
```

- Correlation and Balloon plots for xCell-derived clusters and molecular subtypes. Only subtypes with >= 5 samples were plotted. 
- Violin plot of xCell-derived clusters and TMB (shows kruskal-wallis p-value):

```
plots/xcell_cluster_analysis
├── xcell_clusters_vs_subtype_balloonplot.pdf
├── xcell_clusters_vs_subtype_corrplot.pdf
└── xcell_clusters_vs_tmb.pdf 
```

### 4. xCell-derived cluster analysis

`04-vtest_plots.R`: Generates radar plots of top 5 and bottom 5 immune cell types per cluster.

#### Run script

```
Rscript 04-vtest_plots.R --help

Options:
  --vtest_scores=VTEST_SCORES
    v-test scores file (.tsv)

  --plots_dir=PLOTS_DIR
    output directory to save plots

  -h, --help
    Show this help message and exit
```

#### Output

- Radar plots for top 5 and bottom 5 features per cluster.
- Circular barplots for top 5 and bottom 5 features per cluster.
- Bubble plots for top 5 and bottom 5 features per cluster. 

```
plots/vtest_analysis
├── cluster1_vtest_down.pdf
├── cluster1_vtest_up.pdf
├── cluster2_vtest_down.pdf
├── cluster2_vtest_up.pdf
├── cluster3_vtest_down.pdf
├── cluster3_vtest_up.pdf
├── vtest_bubble_plot.pdf
├── vtest_down_per_cluster.pdf
└── vtest_up_per_cluster.pdf
```

### 5. TIS analysis 

`05-tis_analysis.R`: Compute TIS score from gene expression (count) data and generate violin plot and heatmap of TIS score per cluster. 

#### Run script

```
Rscript 05-tis_analysis.R --help

Options:
  --count_file=COUNT_FILE
    count matrix file (.rds) 

  --xcell_file=XCELL_FILE
    xcell score + cluster file (.tsv)

  --tis_genes=TIS_GENES
    Tumor inflammation signature genes file (.tsv)

  --output_dir=OUTPUT_DIR
    output directory to write files

  --plots_dir=PLOTS_DIR
    output directory to save plots

  -h, --help
    Show this help message and exit
```

#### Output

TSV file with TIS scores for each sample:

```
results/tis_analysis
└── TIS_scores.tsv
```

- Violin plot of average and total TIS scores per cluster.
- Heatmap of TIS score annotated by cluster and `molecular_subtype`. Two normalization methods were used for representation: 1) quantile normalization and 2) z-score normalization. Only subtypes with >= 5 samples were plotted. 

```
plots/tis_analysis
├── TIS_cluster_avg_violinplot.pdf
├── TIS_cluster_sum_violinplot.pdf
├── TIS_heatmap_quantile_norm.pdf
└── TIS_heatmap_zscore.pdf
```

### 6. I) Survival analysis

`06-survival_analysis.R`: Kaplan-Meier and Cox regression Survival analysis of OS and EFS for xCell-derived clusters.

#### Run script

```
Rscript --vanilla 06-survival_analysis.R
```

#### Output

Summary of survival analysis:

```
results/survival_analysis
├── coxph_EFS_model_summary.tsv
└── coxph_OS_model_summary.tsv
```

 - Kaplan Meier curves of OS and EFS for xCell-derived clusters
 - Cox regression (forest plots) of OS and EFS for xCell-derived clusters

```
plots/survival_analysis
├── KM_EFS_forestplot_lgg_clusters.pdf
├── KM_EFS_plot_lgg_clusters.pdf
├── KM_OS_forestplot_lgg_clusters.pdf
└── KM_OS_plot_lgg_clusters.pdf
```

### 6. II) Survival analysis with TIS

`06-survival_analysis_tis.R`: Kaplan-Meier and Cox regression Survival analysis of OS and EFS for TIS groups. TIS scores are divided into two groups: high and low  based on the median value.

#### Run script

```
Rscript --vanilla 06-survival_analysis_tis.R
```

#### Output

Summary of survival analysis:

```
results/survival_analysis
├── coxph_EFS_model_summary_tis.tsv
├── coxph_OS_model_summary_tis.tsv
└── lgat_tis_groups.txt
```

 - Kaplan Meier curves of OS and EFS for TIS groups
 - Cox regression (forest plots) of OS and EFS for TIS groups

```
plots/survival_analysis
├── KM_EFS_forestplot_lgg_tis.pdf
├── KM_EFS_plot_lgg_tis.pdf
├── KM_OS_forestplot_lgg_tis.pdf
└── KM_OS_plot_lgg_tis.pdf
```

### 7. Tumor Purity analysis

`07-tumor_purity_scatter.R`: Scatter plots for tumor purity and immune score.

#### Run script

```
Rscript 07-tumor_purity_scatter.R --help

Options:
  --mat=MAT
    input matrix (.rds)

  --histology=HISTOLOGY
    histology file (.tsv)

  --plots_dir=PLOTS_DIR
    output directory for plots

  -h, --help
    Show this help message and exit
```

#### Output

Scatter plots showing correlation between xCell-derived immune scores and 
- tumor_fraction
- tumor_fraction_RFpurify_ABSOLUTE
- tumor_fraction_RFpurify_ESTIMATE
- tumor_fraction_LUMP

```
plots/tumor_purity_analysis

├── immunescore_RFpurify_ABSOLUTE_correlation_plot.pdf
├── immunescore_RFpurify_ESTIMATE_correlation_plot.pdf
├── immunescore_tumor_fraction_LUMP_correlation_plot.pdf
└── immunescore_tumorfraction_correlation_plot.pdf
```

### 8. Immune Signature analysis

`08-immuneSig_db_gsva.R`: Using TPM gene expression data and `MSigDB C7: immunologic` signature gene sets, compute differential regulated pathways between each pair of xCell-derived clusters.

#### Run script

```
Rscript 08-immuneSig_db_gsva.R --help

Options:
  --mat=MAT
    gene expression matrix preferably TPM (.rds)

  --xcell_file=XCELL_FILE
    xcell score + cluster file (.tsv)

  --output_dir=OUTPUT_DIR
    output directory for files

  --plots_dir=PLOTS_DIR
    output directory for plots

  -h, --help
    Show this help message and exit
```

#### Output

GSVA and limma output for differentially regulated pathways:

```
results/gsva_analysis
├── gsva_igg_cluster_result.RData
└── limma_test_results.txt
```

Heatmap for the top 50 differentially expressed pathways:

```
plots/gsva_analysis
└── gsva_heatmap.pdf
```
