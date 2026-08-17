# scLC_ICI_PBMC

本仓库汇总 scLC ICI PBMC 项目的分析脚本与 PNG 结果图。三个分析子项目按照统一结构组织，每个子项目根目录仅保留最新日期，历史日期统一存放在 `Archive/` 中。

## 目录约定

每个日期目录均包含以下三个同级目录：

- `Results/`：Git 仅跟踪 PNG 结果图；服务器运行时可在其中生成被忽略的中间文件。
- `Scripts/`：保存 Python、R、Shell 等分析脚本。
- `Supplementary_materials/`：保存运行所需的小型补充材料；空目录使用 `.gitkeep` 保留。

当前日期组织如下：

- **Methscan**：最新批次 `20260815`；归档批次 `20260716`、`20260718`。
- **MethylVI**：最新批次 `20260816`；归档批次 `20260810`、`20260813`。
- **Scanpy**：最新批次 `20260815`；归档批次 `20260810`。

仓库当前包含 175 个分析脚本和 223 张 PNG 结果图。三个原始子仓库不在本仓库的跟踪范围内，其远程仓库配置保持不变。

### Methscan 当前流程顺序

`Methscan/20260815/Scripts/01_Upstream` 使用 `00–08` 表示执行顺序：`00` 为公共配置，`01–03` 完成 cov 检查、概率去重和上游处理，`04` 运行逐样本 DMR，`05–07` 完成 Top200 DMR 选择、矩阵计算和绘图，`08` 处理无 null DMR 时的 raw-p 回退。带 `a/b` 后缀的文件是同一步骤调用的辅助程序，不需要单独调整执行顺序。`01`、`02` 和 `04` 均通过单一入口同时提供批处理与单样本模式，不再保留单独的单样本 Shell 实现。原 `04_run_all_samples_to_smooth.sh` 已合并为 `03_run_upstream_pipeline.sh smooth-all [sample_jobs]` 和 `smooth-status`。

## 服务器部署与路径配置

服务器上的标准部署目录为：

```text
/share/home/rzli/scLC_ICI_PBMC
```

仓库根目录的 `project_config.sh` 是三个当前流程共用的唯一路径入口。它会根据自身位置解析绝对项目根目录，并统一导出数据、参考文件、Conda、三个项目日期以及各项目 `Scripts/`、`Results/` 路径。默认外部路径如下：

- 数据：`/share/LCZX_Data/data`
- 参考文件：`/share/LCZX_Data/ref`
- Conda：`/share/home/rzli/miniconda3`

部署后可检查配置：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
source project_config.sh
printf '%s\n' "$SCLC_PROJECT_ROOT" "$SCLC_SCANPY_RESULTS"
```

若以后移动仓库或数据目录，只需在加载配置前导出对应的 `SCLC_*` 环境变量。当前 Methscan、MethylVI 和 Scanpy 启动脚本会自动加载该配置；`Archive/` 中的历史脚本保留原始路径，不作为当前服务器流程入口。

## 文件树

为便于阅读，`Results/` 目录仅显示目录节点，不展开其中的 PNG 文件和子目录。

```text
.
├── .gitignore
├── README.md
├── project_config.sh
├── Methscan/
│   ├── 20260815/
│   │   ├── Results/
│   │   ├── Scripts/
│   │   │   ├── README.md
│   │   │   ├── 01_Upstream/
│   │   │   │   ├── 00_workflow_common.sh
│   │   │   │   ├── 01_check_cov_duplicates.sh
│   │   │   │   ├── 02_deduplicate_cov_by_probability.sh
│   │   │   │   ├── 03_run_upstream_pipeline.sh
│   │   │   │   ├── 03a_build_scanpy_clean_cell_list.py
│   │   │   │   ├── 04_run_celltype_dmr.sh
│   │   │   │   ├── 05_select_top200_dmrs.sh
│   │   │   │   ├── 05a_extract_celltype_hypo_dmrs.py
│   │   │   │   ├── 05b_merge_sample_dmrs.py
│   │   │   │   ├── 06_compute_top200_dmr_matrix.sh
│   │   │   │   ├── 06a_compute_dmr_mean_of_cpg_ratios.py
│   │   │   │   ├── 07_plot_all_top200_heatmaps.sh
│   │   │   │   ├── 07a_plot_single_cell_dmr_heatmaps.py
│   │   │   │   └── 08_rerun_rawp_no_null_fdr.sh
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
│   │       ├── README.md
│   │       └── .gitkeep
│   └── Archive/
│       ├── 20260716/
│       │   ├── Results/
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
│   │   ├── Scripts/
│   │   │   ├── README.md
│   │   │   ├── 00_config.sh
│   │   │   ├── 02_prepare_allcools.sh
│   │   │   ├── 03_cluster_allcools.py
│   │   │   ├── 04_verify_inputs.py
│   │   │   ├── 05_build_methylvi_input.py
│   │   │   ├── 06_train_methylvi.py
│   │   │   ├── 07_plot_embeddings.py
│   │   │   ├── 08_plot_supervised_umap.py
│   │   │   ├── 09_run_pipeline.sh
│   │   │   ├── 10_plot_sequencing_depth.py
│   │   │   ├── 11_compare_qc_cell_sets.py
│   │   │   ├── 12_plot_cpg_sites.py
│   │   │   ├── mvi_utils.py
│   │   │   └── tests/
│   │   │       ├── test_methylvi_smoke.py
│   │   │       └── test_mvi_utils.py
│   │   └── Supplementary_materials/
│   │       ├── README.md
│   │       ├── 01_sample_metadata.tsv
│   │       ├── ENCFF356LFX_GRCh38_blacklist.bed.gz
│   │       └── hg38.canonical.chrom.sizes
│   └── Archive/
│       ├── 20260810/
│       │   ├── Results/
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
│           ├── Scripts/
│           │   └── .gitkeep
│           └── Supplementary_materials/
│               └── .gitkeep
└── Scanpy/
    ├── 20260815/
    │   ├── Results/
    │   ├── Scripts/
    │   │   ├── README.md
    │   │   ├── 00_config.sh
    │   │   ├── 01_integration.py
    │   │   ├── 02_annotation_config.py
    │   │   ├── 03_annotation.py
    │   │   ├── 04_export_figures.py
    │   │   ├── 05_run_integration.sh
    │   │   ├── 06_run_annotation.sh
    │   │   └── 07_run_export_figures.sh
    │   └── Supplementary_materials/
    │       ├── README.md
    │       └── .gitkeep
    └── Archive/
        └── 20260810/
            ├── Results/
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
