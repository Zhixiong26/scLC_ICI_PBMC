# scLC_ICI_PBMC

本仓库汇总 scLC ICI PBMC 项目的分析脚本与 PNG 结果图。三个分析子项目按照统一结构组织，每个子项目根目录仅保留最新日期，历史日期统一存放在 `Archive/` 中。

## 目录约定

每个日期目录均包含以下三个同级目录：

- `Results/`：仅保存 PNG 结果图。
- `Scripts/`：保存 Python、R、Shell 等分析脚本。
- `Supplementary_materials/`：预留补充材料；空目录使用 `.gitkeep` 保留。

当前日期组织如下：

- **Methscan**：最新批次 `20260815`；归档批次 `20260716`、`20260718`。
- **MethylVI**：最新批次 `20260816`；归档批次 `20260810`、`20260813`。
- **Scanpy**：最新批次 `20260815`；归档批次 `20260810`。

仓库当前包含 177 个分析脚本和 223 张 PNG 结果图。三个原始子仓库不在本仓库的跟踪范围内，其远程仓库配置保持不变。

## 完整文件树

```text
.
├── .gitignore
├── README.md
├── Methscan/
│   ├── 20260815/
│   │   ├── Results/
│   │   │   ├── 01_Upstream/
│   │   │   │   ├── DMRtypeMean_Zscore/
│   │   │   │   │   ├── DMRtypeMean_Zscore__IR01.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__IR02.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__IR03.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__IR04.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__IR05.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__NR01.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__NR02.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__NR03.png
│   │   │   │   │   ├── DMRtypeMean_Zscore__NR04.png
│   │   │   │   │   └── DMRtypeMean_Zscore__NR05.png
│   │   │   │   ├── DMRwise_Zscore/
│   │   │   │   │   ├── DMRwise_Zscore__IR01.png
│   │   │   │   │   ├── DMRwise_Zscore__IR02.png
│   │   │   │   │   ├── DMRwise_Zscore__IR03.png
│   │   │   │   │   ├── DMRwise_Zscore__IR04.png
│   │   │   │   │   ├── DMRwise_Zscore__IR05.png
│   │   │   │   │   ├── DMRwise_Zscore__NR01.png
│   │   │   │   │   ├── DMRwise_Zscore__NR02.png
│   │   │   │   │   ├── DMRwise_Zscore__NR03.png
│   │   │   │   │   ├── DMRwise_Zscore__NR04.png
│   │   │   │   │   └── DMRwise_Zscore__NR05.png
│   │   │   │   ├── DMRwise_Zscore_ColorClip1/
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__IR01.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__IR02.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__IR03.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__IR04.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__IR05.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__NR01.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__NR02.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__NR03.png
│   │   │   │   │   ├── DMRwise_Zscore_ColorClip1__NR04.png
│   │   │   │   │   └── DMRwise_Zscore_ColorClip1__NR05.png
│   │   │   │   └── MeanRatio/
│   │   │   │       ├── MeanRatio__IR01.png
│   │   │   │       ├── MeanRatio__IR02.png
│   │   │   │       ├── MeanRatio__IR03.png
│   │   │   │       ├── MeanRatio__IR04.png
│   │   │   │       ├── MeanRatio__IR05.png
│   │   │   │       ├── MeanRatio__NR01.png
│   │   │   │       ├── MeanRatio__NR02.png
│   │   │   │       ├── MeanRatio__NR03.png
│   │   │   │       ├── MeanRatio__NR04.png
│   │   │   │       └── MeanRatio__NR05.png
│   │   │   └── Archive/
│   │   │       └── root_workflows/
│   │   │           └── Meth diff/
│   │   │               └── result/
│   │   │                   └── clean_celltype_IR_vs_NR/
│   │   │                       ├── top1/
│   │   │                       │   ├── top1_PCA.png
│   │   │                       │   ├── top1_UMAP_by_cell_type.png
│   │   │                       │   ├── top1_UMAP_by_leiden.png
│   │   │                       │   ├── top1_UMAP_by_response.png
│   │   │                       │   ├── top1_UMAP_by_sample.png
│   │   │                       │   └── top1_UMAP_response_by_cell_type.png
│   │   │                       ├── top2/
│   │   │                       │   ├── top2_PCA.png
│   │   │                       │   ├── top2_UMAP_by_cell_type.png
│   │   │                       │   ├── top2_UMAP_by_leiden.png
│   │   │                       │   ├── top2_UMAP_by_response.png
│   │   │                       │   ├── top2_UMAP_by_sample.png
│   │   │                       │   └── top2_UMAP_response_by_cell_type.png
│   │   │                       └── top5/
│   │   │                           ├── top5_PCA.png
│   │   │                           ├── top5_UMAP_by_cell_type.png
│   │   │                           ├── top5_UMAP_by_leiden.png
│   │   │                           ├── top5_UMAP_by_response.png
│   │   │                           ├── top5_UMAP_by_sample.png
│   │   │                           └── top5_UMAP_response_by_cell_type.png
│   │   ├── Scripts/
│   │   │   ├── 01_Upstream/
│   │   │   │   ├── 01_check_cov_duplicates.sh
│   │   │   │   ├── 02_deduplicate_cov_by_probability.sh
│   │   │   │   ├── 03_run_upstream_pipeline.sh
│   │   │   │   ├── 04_run_all_samples_to_smooth.sh
│   │   │   │   ├── 05_run_all_samples_dmr.sh
│   │   │   │   ├── 06_select_top200_dmrs.sh
│   │   │   │   ├── 07_compute_top200_dmr_matrix.sh
│   │   │   │   ├── 08_plot_all_top200_heatmaps.sh
│   │   │   │   ├── 09_rerun_rawp_no_null_fdr.sh
│   │   │   │   └── lib/
│   │   │   │       ├── build_scanpy_clean_cell_list.py
│   │   │   │       ├── check_cov_duplicates_one_sample.sh
│   │   │   │       ├── deduplicate_cov_by_probability_one_sample.sh
│   │   │   │       ├── workflow_common.sh
│   │   │   │       └── methdiff/
│   │   │   │           ├── run_single_sample_dmr.sh
│   │   │   │           └── python/
│   │   │   │               ├── 02_merge_sample_dmrs.py
│   │   │   │               ├── 04_plot_single_cell_dmr_heatmaps.py
│   │   │   │               ├── 05_extract_celltype_hypo_dmrs_top1500.py
│   │   │   │               └── 06_compute_dmr_mean_of_cpg_ratios.py
│   │   │   └── Archive/
│   │   │       ├── 01_Upstream/
│   │   │       │   ├── 10_prepare_merged_response_input.sh
│   │   │       │   ├── 11_run_merged_celltype_ir_vs_nr_dmr.sh
│   │   │       │   ├── 12_select_merged_ir_nr_candidate_dmrs.sh
│   │   │       │   ├── 13_run_merged_response_scan_matrix.sh
│   │   │       │   └── archive_legacy/
│   │   │       │       ├── 04_run_ir01_single_sample_dmr_entry_legacy.sh
│   │   │       │       ├── Methscan_10samples_merged_batch_correct_legacy.sh
│   │   │       │       ├── Methscan_full_legacy.sh
│   │   │       │       ├── run_upstream_pipeline_legacy.sh
│   │   │       │       ├── filter_threshold_audit_legacy/
│   │   │       │       │   └── 90_summarize_filter_thresholds.py
│   │   │       │       └── plot_wrappers_before_unified_20260812/
│   │   │       │           ├── 08_plot_top200_dmr_heatmap.sh
│   │   │       │           ├── 09_plot_top200_dmr_zscore_heatmap.sh
│   │   │       │           ├── 09b_plot_top200_dmr_type_mean_zscore_heatmap.sh
│   │   │       │           ├── 09c_plot_top200_dmr_zscore_maxabs_heatmap.sh
│   │   │       │           ├── 09d_plot_top200_dmr_zscore_clip1_heatmap.sh
│   │   │       │           ├── 09e_plot_top200_dmr_type_mean_zscore_maxabs_heatmap.sh
│   │   │       │           └── 09f_plot_top200_dmr_type_mean_zscore_clip1_heatmap.sh
│   │   │       ├── 02_Methdiff/
│   │   │       │   ├── archive_merged_workflow/
│   │   │       │   │   └── run_methdiff_pipeline_merged_legacy.sh
│   │   │       │   └── Result/
│   │   │       │       ├── 01_extract_celltype_hypo_dmrs.py
│   │   │       │       ├── 03_compute_dmr_mean_cpg_ratio.py
│   │   │       │       ├── 07_audit_sparse_value_origin.py
│   │   │       │       ├── 08_trace_sparse_values_to_cov.py
│   │   │       │       └── run_ir01_top200_heatmap.sh
│   │   │       └── root_workflows/
│   │   │           ├── 03_MethExprBubble/
│   │   │           │   ├── 01_audit_promoter_dmr_expression_inputs.py
│   │   │           │   ├── 02_map_ir_hypo_dmrs_to_promoters.py
│   │   │           │   ├── 03_compute_promoter_methylation_pseudobulk.py
│   │   │           │   ├── 04_compute_rna_expression_pseudobulk.py
│   │   │           │   ├── 05_correlate_promoter_methylation_expression.py
│   │   │           │   ├── 06_plot_promoter_effect_quadrant.py
│   │   │           │   ├── run_promoter_dmr_expression_workflow.sh
│   │   │           │   └── workflow_config.py
│   │   │           ├── 04_Hamming_distance/
│   │   │           │   ├── hamming_scwgbs.py
│   │   │           │   └── run_hamming_pipeline.sh
│   │   │           ├── 05_VMR_clustering/
│   │   │           │   ├── run_vmr_clustering.sh
│   │   │           │   └── vmr_clustering.R
│   │   │           ├── 06_PromoterDMR_eQTM/
│   │   │           │   ├── 01_map_bidirectional_dmrs_to_promoters.py
│   │   │           │   ├── 02_compute_promoter_methylation_pseudobulk.py
│   │   │           │   ├── 03_compute_rna_expression_pseudobulk.py
│   │   │           │   ├── 04_compute_pseudobulk_correlations.py
│   │   │           │   ├── run_pseudobulk_correlation.sh
│   │   │           │   └── workflow_config.py
│   │   │           └── DMRs/
│   │   │               └── scripts/
│   │   │                   ├── common/
│   │   │                   │   └── make_promoter_bed_2kb_2kb.py
│   │   │                   ├── q001_2kb_2kb/
│   │   │                   │   ├── 00_make_dna_cell_metadata.py
│   │   │                   │   ├── 01_merge_subtract_dmr.sh
│   │   │                   │   ├── 02_subset_vmr_matrix.py
│   │   │                   │   ├── 03_run_subset_vmr_matrix.sh
│   │   │                   │   ├── 04_summarize_ir_vs_nr_dmr.sh
│   │   │                   │   ├── 05_filter_ir_vs_nr_dmr.sh
│   │   │                   │   ├── 06_make_promoter_bed_2kb_2kb.py
│   │   │                   │   ├── 07_intersect_dmr_promoter.sh
│   │   │                   │   ├── 08_format_promoter_dmr_gene_table.py
│   │   │                   │   ├── 09_filter_protein_coding_genes.py
│   │   │                   │   ├── 10_integrate_pseudobulk_expression.py
│   │   │                   │   ├── 11_run_pseudobulk_expression.sh
│   │   │                   │   └── 12_correlate_dna_methylation_rna_expression.py
│   │   │                   └── q005_2kb_2kb/
│   │   │                       ├── 00_make_dna_cell_metadata.py
│   │   │                       ├── 01_merge_subtract_dmr.sh
│   │   │                       ├── 02_subset_vmr_matrix.py
│   │   │                       ├── 03_run_subset_vmr_matrix.sh
│   │   │                       ├── 04_summarize_ir_vs_nr_dmr.sh
│   │   │                       ├── 05_filter_ir_vs_nr_dmr.sh
│   │   │                       ├── 06_make_promoter_bed_2kb_2kb.py
│   │   │                       ├── 07_intersect_dmr_promoter.sh
│   │   │                       ├── 08_format_promoter_dmr_gene_table.py
│   │   │                       ├── 09_filter_protein_coding_genes.py
│   │   │                       ├── 10_integrate_pseudobulk_expression.py
│   │   │                       ├── 11_run_pseudobulk_expression.sh
│   │   │                       └── 12_correlate_dna_methylation_rna_expression.py
│   │   └── Supplementary_materials/
│   │       └── .gitkeep
│   └── Archive/
│       ├── 20260716/
│       │   ├── Results/
│       │   │   ├── .gitkeep
│       │   │   └── Annotation/
│       │   │       └── result/
│       │   │           └── plots/
│       │   │               ├── ALL_PCA_200k.png
│       │   │               ├── ALL_umap_plot_by_cell_type_200k.png
│       │   │               ├── ALL_umap_plot_by_leiden_200k.png
│       │   │               ├── ALL_umap_plot_by_response_200k.png
│       │   │               ├── ALL_umap_plot_by_sample_200k.png
│       │   │               └── ALL_umap_plot_response_by_cell_type_200k.png
│       │   ├── Scripts/
│       │   │   ├── Annotation/
│       │   │   │   └── scripts/
│       │   │   │       ├── 01_run_All_200k.sh
│       │   │   │       └── 02_All_200k_analysis.R
│       │   │   └── Meth_diff/
│       │   │       ├── 03_meth_diff_celltype_sample_pairwise_200k.sh
│       │   │       ├── 04_generate_celltype_sample_pairwise_groups.R
│       │   │       ├── 08_merge_celltype_sample_pairwise_dmrs.sh
│       │   │       ├── 10_prepare_individual_effect_mask.sh
│       │   │       ├── 11_run_threshold_vmrs_remove_individual.sh
│       │   │       ├── 12_build_threshold_before_vmr_matrix.sh
│       │   │       ├── 13_recluster_threshold_clean_vmrs.R
│       │   │       ├── 13_run_threshold_before_vmr_reclustering.sh
│       │   │       ├── 13_run_threshold_clean_vmr_reclustering.sh
│       │   │       ├── 14_collect_threshold_metrics.sh
│       │   │       ├── validate_threshold_matrix.py
│       │   │       └── archive_old_workflow/
│       │   │           ├── 01_meth_diff_pairwise_200k.sh
│       │   │           ├── 02_generate_cell_groups.R
│       │   │           ├── 05_meth_diff_sample_pairwise_200k.sh
│       │   │           ├── 06_generate_sample_pairwise_groups.R
│       │   │           ├── 07_merge_celltype_ir_vs_nr_dmrs.sh
│       │   │           ├── 09_subtract_within_group_sample_dmrs_from_ir_nr_dmrs.sh
│       │   │           ├── 10_map_clean_dmrs_to_all_vmrs.sh
│       │   │           ├── 11_select_top_clean_vmrs.sh
│       │   │           ├── 13_run_top_clean_vmr_reclustering.sh
│       │   │           ├── run_steps_07_09.sh
│       │   │           ├── run_steps_10_11.sh
│       │   │           ├── run_steps_12_13.sh
│       │   │           └── legacy/
│       │   │               ├── 12_run_subset_top_clean_vmr_matrices.sh
│       │   │               └── 12_subset_top_clean_vmr_matrices.py
│       │   └── Supplementary_materials/
│       │       └── .gitkeep
│       └── 20260718/
│           ├── Results/
│           │   └── Meth_diff/
│           │       └── result/
│           │           └── UMAP_before_after_all/
│           │               ├── cell_type/
│           │               │   ├── IR_threshold001_after_UMAP_by_cell_type.png
│           │               │   ├── IR_threshold001_before_UMAP_by_cell_type.png
│           │               │   ├── IR_threshold002_after_UMAP_by_cell_type.png
│           │               │   ├── IR_threshold002_before_UMAP_by_cell_type.png
│           │               │   ├── IR_threshold005_after_UMAP_by_cell_type.png
│           │               │   ├── IR_threshold005_before_UMAP_by_cell_type.png
│           │               │   ├── NR_threshold001_after_UMAP_by_cell_type.png
│           │               │   ├── NR_threshold001_before_UMAP_by_cell_type.png
│           │               │   ├── NR_threshold002_after_UMAP_by_cell_type.png
│           │               │   ├── NR_threshold002_before_UMAP_by_cell_type.png
│           │               │   ├── NR_threshold005_after_UMAP_by_cell_type.png
│           │               │   └── NR_threshold005_before_UMAP_by_cell_type.png
│           │               └── sample/
│           │                   ├── IR_threshold001_after_UMAP_by_sample.png
│           │                   ├── IR_threshold001_before_UMAP_by_sample.png
│           │                   ├── IR_threshold002_after_UMAP_by_sample.png
│           │                   ├── IR_threshold002_before_UMAP_by_sample.png
│           │                   ├── IR_threshold005_after_UMAP_by_sample.png
│           │                   ├── IR_threshold005_before_UMAP_by_sample.png
│           │                   ├── NR_threshold001_after_UMAP_by_sample.png
│           │                   ├── NR_threshold001_before_UMAP_by_sample.png
│           │                   ├── NR_threshold002_after_UMAP_by_sample.png
│           │                   ├── NR_threshold002_before_UMAP_by_sample.png
│           │                   ├── NR_threshold005_after_UMAP_by_sample.png
│           │                   └── NR_threshold005_before_UMAP_by_sample.png
│           ├── Scripts/
│           │   └── Meth_diff/
│           │       ├── 01_prepare_response_data.sh
│           │       ├── 02_scan_response_vmrs.sh
│           │       ├── 03_prepare_individual_effect_dmr_union.sh
│           │       ├── 04_build_response_clean_vmr_matrix.sh
│           │       ├── 04b_build_response_before_vmr_matrix.sh
│           │       ├── 05_recluster_response_clean_vmrs.R
│           │       ├── 05_run_response_reclustering.sh
│           │       ├── 06_collect_response_metrics.sh
│           │       ├── 06_compare_before_after_metrics.py
│           │       ├── 07_filter_celltype_dmrs.py
│           │       ├── 07_prepare_filtered_celltype_dmrs.sh
│           │       ├── 08_build_filtered_celltype_dmr_matrix.sh
│           │       ├── 09_recluster_filtered_celltype_dmrs.R
│           │       ├── 09_run_filtered_celltype_dmr_reclustering.sh
│           │       ├── 10_prepare_matched_individual_effect_clean_dmrs.sh
│           │       ├── validate_response_matrix.py
│           │       └── archive_old_workflow/
│           │           ├── 01_meth_diff_pairwise_200k.sh
│           │           ├── 02_generate_cell_groups.R
│           │           ├── 03_meth_diff_celltype_sample_pairwise_200k.sh
│           │           ├── 04_generate_celltype_sample_pairwise_groups.R
│           │           ├── 05_meth_diff_sample_pairwise_200k.sh
│           │           ├── 06_generate_sample_pairwise_groups.R
│           │           ├── 07_merge_celltype_ir_vs_nr_dmrs.sh
│           │           ├── 08_merge_celltype_sample_pairwise_dmrs.sh
│           │           ├── 09_subtract_within_group_sample_dmrs_from_ir_nr_dmrs.sh
│           │           ├── 10_map_clean_dmrs_to_all_vmrs.sh
│           │           ├── 10_prepare_individual_effect_mask.sh
│           │           ├── 11_run_threshold_vmrs_remove_individual.sh
│           │           ├── 11_select_top_clean_vmrs.sh
│           │           ├── 13_recluster_threshold_clean_vmrs.R
│           │           ├── 13_run_threshold_clean_vmr_reclustering.sh
│           │           ├── 13_run_top_clean_vmr_reclustering.sh
│           │           ├── 14_collect_threshold_metrics.sh
│           │           ├── run_steps_07_09.sh
│           │           ├── run_steps_10_11.sh
│           │           ├── run_steps_12_13.sh
│           │           ├── validate_threshold_matrix.py
│           │           └── legacy/
│           │               ├── 12_run_subset_top_clean_vmr_matrices.sh
│           │               └── 12_subset_top_clean_vmr_matrices.py
│           └── Supplementary_materials/
│               └── .gitkeep
├── MethylVI/
│   ├── 20260816/
│   │   ├── Results/
│   │   │   ├── blacklist_f0p2/
│   │   │   │   └── 01_before_methylvi/
│   │   │   │       ├── allcools_5kb_tsne_L1.png
│   │   │   │       └── allcools_5kb_umap_L1.png
│   │   │   ├── blacklist_f0p2_100k/
│   │   │   │   └── 01_before_methylvi/
│   │   │   │       ├── allcools_5kb_tsne_L1.png
│   │   │   │       └── allcools_5kb_umap_L1.png
│   │   │   ├── blacklist_f0p2_4819/
│   │   │   │   ├── 01_before_methylvi/
│   │   │   │   │   ├── allcools_original_embedding_cell_type.png
│   │   │   │   │   ├── allcools_original_embedding_condition.png
│   │   │   │   │   └── allcools_original_embedding_sample_id.png
│   │   │   │   ├── 02_after_methylvi/
│   │   │   │   │   ├── methylvi_umap_cell_type.png
│   │   │   │   │   ├── methylvi_umap_condition.png
│   │   │   │   │   └── methylvi_umap_sample_id.png
│   │   │   │   └── 03_supervised_umap/
│   │   │   │       ├── target_weight_0p2/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       ├── target_weight_0p5/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       ├── target_weight_0p7/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       └── target_weight_0p9/
│   │   │   │           ├── methylvi_supervised_umap_cell_type.png
│   │   │   │           ├── methylvi_supervised_umap_condition.png
│   │   │   │           ├── methylvi_supervised_umap_sample_id.png
│   │   │   │           ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │           └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   ├── blacklist_f0p2_4819_latent15/
│   │   │   │   ├── 01_before_methylvi/
│   │   │   │   │   ├── allcools_original_embedding_cell_type.png
│   │   │   │   │   ├── allcools_original_embedding_condition.png
│   │   │   │   │   └── allcools_original_embedding_sample_id.png
│   │   │   │   ├── 02_after_methylvi/
│   │   │   │   │   ├── methylvi_umap_cell_type.png
│   │   │   │   │   ├── methylvi_umap_condition.png
│   │   │   │   │   └── methylvi_umap_sample_id.png
│   │   │   │   └── 03_supervised_umap/
│   │   │   │       ├── target_weight_0p2/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       ├── target_weight_0p5/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       ├── target_weight_0p7/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       └── target_weight_0p9/
│   │   │   │           ├── methylvi_supervised_umap_cell_type.png
│   │   │   │           ├── methylvi_supervised_umap_condition.png
│   │   │   │           ├── methylvi_supervised_umap_sample_id.png
│   │   │   │           ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │           └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   ├── blacklist_f0p2_4819_latent25/
│   │   │   │   ├── 01_before_methylvi/
│   │   │   │   │   ├── allcools_original_embedding_cell_type.png
│   │   │   │   │   ├── allcools_original_embedding_condition.png
│   │   │   │   │   └── allcools_original_embedding_sample_id.png
│   │   │   │   ├── 02_after_methylvi/
│   │   │   │   │   ├── methylvi_umap_cell_type.png
│   │   │   │   │   ├── methylvi_umap_condition.png
│   │   │   │   │   └── methylvi_umap_sample_id.png
│   │   │   │   └── 03_supervised_umap/
│   │   │   │       ├── target_weight_0p2/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       ├── target_weight_0p5/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       ├── target_weight_0p7/
│   │   │   │       │   ├── methylvi_supervised_umap_cell_type.png
│   │   │   │       │   ├── methylvi_supervised_umap_condition.png
│   │   │   │       │   ├── methylvi_supervised_umap_sample_id.png
│   │   │   │       │   ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │       │   └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   │       └── target_weight_0p9/
│   │   │   │           ├── methylvi_supervised_umap_cell_type.png
│   │   │   │           ├── methylvi_supervised_umap_condition.png
│   │   │   │           ├── methylvi_supervised_umap_sample_id.png
│   │   │   │           ├── methylvi_supervised_umap_sequencing_depth_absolute.png
│   │   │   │           └── methylvi_supervised_umap_sequencing_depth.png
│   │   │   ├── blacklist_f0p2_50k/
│   │   │   │   └── 01_before_methylvi/
│   │   │   │       ├── allcools_5kb_tsne_L1.png
│   │   │   │       └── allcools_5kb_umap_L1.png
│   │   │   └── qc_comparison_6199_vs_4819/
│   │   │       ├── target_weight_0p2/
│   │   │       │   ├── qc_retained_vs_removed.png
│   │   │       │   ├── sequencing_depth_absolute_removed_cells_overlay.png
│   │   │       │   └── sequencing_depth_removed_cells_overlay.png
│   │   │       ├── target_weight_0p5/
│   │   │       │   ├── qc_retained_vs_removed.png
│   │   │       │   ├── sequencing_depth_absolute_removed_cells_overlay.png
│   │   │       │   └── sequencing_depth_removed_cells_overlay.png
│   │   │       ├── target_weight_0p7/
│   │   │       │   ├── qc_retained_vs_removed.png
│   │   │       │   ├── sequencing_depth_absolute_removed_cells_overlay.png
│   │   │       │   └── sequencing_depth_removed_cells_overlay.png
│   │   │       └── target_weight_0p9/
│   │   │           ├── qc_retained_vs_removed.png
│   │   │           ├── sequencing_depth_absolute_removed_cells_overlay.png
│   │   │           └── sequencing_depth_removed_cells_overlay.png
│   │   ├── Scripts/
│   │   │   └── .gitkeep
│   │   └── Supplementary_materials/
│   │       └── .gitkeep
│   └── Archive/
│       ├── 20260810/
│       │   ├── Results/
│       │   │   └── .gitkeep
│       │   ├── Scripts/
│       │   │   ├── 00_config.sh
│       │   │   ├── 02_prepare_allcools.sh
│       │   │   ├── 03_cluster_allcools.py
│       │   │   ├── 04_verify_inputs.py
│       │   │   ├── 05_build_methylvi_input.py
│       │   │   ├── 06_train_methylvi.py
│       │   │   ├── 07_plot_embeddings.py
│       │   │   ├── 08_plot_supervised_umap.py
│       │   │   ├── 09_run_pipeline.sh
│       │   │   ├── 10_plot_sequencing_depth.py
│       │   │   ├── 11_compare_qc_cell_sets.py
│       │   │   ├── 12_plot_cpg_sites.py
│       │   │   ├── mvi_utils.py
│       │   │   ├── tests/
│       │   │   │   ├── test_methylvi_smoke.py
│       │   │   │   └── test_mvi_utils.py
│       │   │   └── yuanpei/
│       │   │       ├── cluster_5kbin.py
│       │   │       └── reproducible_methylVI_pipeline/
│       │   │           ├── 00_verify_inputs.py
│       │   │           ├── config.sh
│       │   │           └── run_pipeline.sh
│       │   └── Supplementary_materials/
│       │       └── .gitkeep
│       └── 20260813/
│           ├── Results/
│           │   ├── blacklist_f0p2_100k/
│           │   │   └── 01_before_methylvi/
│           │   │       ├── allcools_5kb_tsne_L1.png
│           │   │       └── allcools_5kb_umap_L1.png
│           │   └── blacklist_f0p2_50k/
│           │       └── 01_before_methylvi/
│           │           ├── allcools_5kb_tsne_L1.png
│           │           └── allcools_5kb_umap_L1.png
│           ├── Scripts/
│           │   └── .gitkeep
│           └── Supplementary_materials/
│               └── .gitkeep
└── Scanpy/
    ├── 20260815/
    │   ├── Results/
    │   │   ├── 12_umap_final_cell_type_IR01.png
    │   │   ├── 12_umap_final_cell_type_IR02.png
    │   │   ├── 12_umap_final_cell_type_IR03.png
    │   │   ├── 12_umap_final_cell_type_IR04.png
    │   │   ├── 12_umap_final_cell_type_IR05.png
    │   │   ├── 12_umap_final_cell_type_NR01.png
    │   │   ├── 12_umap_final_cell_type_NR02.png
    │   │   ├── 12_umap_final_cell_type_NR03.png
    │   │   ├── 12_umap_final_cell_type_NR04.png
    │   │   ├── 12_umap_final_cell_type_NR05.png
    │   │   └── 13_umap_clean_cells_final_annotation.png
    │   ├── Scripts/
    │   │   └── .gitkeep
    │   └── Supplementary_materials/
    │       └── .gitkeep
    └── Archive/
        └── 20260810/
            ├── Results/
            │   ├── current/
            │   │   ├── 05_umap_by_leiden_integrated_副本.png
            │   │   ├── 06_umap_by_final_cell_type.png
            │   │   ├── 12_umap_final_cell_type_IR01.png
            │   │   ├── 12_umap_final_cell_type_IR02.png
            │   │   ├── 12_umap_final_cell_type_IR03.png
            │   │   ├── 12_umap_final_cell_type_IR04.png
            │   │   ├── 12_umap_final_cell_type_IR05.png
            │   │   ├── 12_umap_final_cell_type_NR01.png
            │   │   ├── 12_umap_final_cell_type_NR02.png
            │   │   ├── 12_umap_final_cell_type_NR03.png
            │   │   ├── 12_umap_final_cell_type_NR04.png
            │   │   ├── 12_umap_final_cell_type_NR05.png
            │   │   └── 13_umap_clean_cells_final_annotation.png
            │   └── Result0810/
            │       ├── 12_umap_final_cell_type_IR01.png
            │       ├── 12_umap_final_cell_type_IR02.png
            │       ├── 12_umap_final_cell_type_IR03.png
            │       ├── 12_umap_final_cell_type_IR04.png
            │       ├── 12_umap_final_cell_type_IR05.png
            │       ├── 12_umap_final_cell_type_NR01.png
            │       ├── 12_umap_final_cell_type_NR02.png
            │       ├── 12_umap_final_cell_type_NR03.png
            │       ├── 12_umap_final_cell_type_NR04.png
            │       ├── 12_umap_final_cell_type_NR05.png
            │       └── 13_umap_clean_cells_final_annotation.png
            ├── Scripts/
            │   ├── 01_integration.py
            │   ├── 02_annotation_config.py
            │   ├── 03_annotation.py
            │   ├── 04_export_figures.py
            │   ├── 05_run_integration.sh
            │   ├── 06_run_annotation.sh
            │   └── 07_run_export_figures.sh
            └── Supplementary_materials/
                └── .gitkeep
```
