# stniche/__init__.py

from .stniche import (
    compute_groupwise_adjacency_matrix, #计算邻接spot的差异
    plot_adj_difference_heatmap, #绘制邻接spot的差异
    extract_and_group_3rd_structures_from_2nd_with_ratio, #从二级到三级结构
    filter_unique_significant_structures, #筛选三级结构
    run_iterative_analysis_over_df, #识别niche
    plot_structure_hex, #绘制niche结构
    highlight_niche_on_spatial, #在切片上绘制niche
    run_niche_differential_and_enrichment, #识别niche功能
    plot_spatial_communication_sankey #识别niche信号通路
)
