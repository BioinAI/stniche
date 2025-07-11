
# stniche

**stniche** is a Python toolkit for identifying and analyzing spatial niches in spatial transcriptomics data. It integrates spatial connectivity, geometric motif extraction, statistical significance, differential gene expression, pathway enrichment, and spatial communication visualization.

---

## 🧬 Features

- Compute neighborhood connectivity differences across groups
- Identify enriched 2nd/3rd/4th-order spatial motifs
- Filter structures by statistical significance, fold change, and coverage
- Iteratively expand niche structures to discover higher-order spatial architecture
- Highlight niche points on Visium histology images
- Perform differential expression and enrichment analysis of niche spots
- Visualize niche communication using Sankey diagrams

---

## 📦 Installation

```bash
pip install .
```

---

## 🚀 Usage Overview

### 1. Compute Adjacency Differences

```python
from stniche import compute_groupwise_adjacency_matrix

adj_df, fdr_filtered_2nd = compute_groupwise_adjacency_matrix(
    adata,
    row_key='array_row',
    col_key='array_col',
    sample_key='sample_id',
    cluster_key='leiden_0.7',
    group_key='class',
    groups=('basal cell carcinoma', 'healthy'),
    focus_group='basal cell carcinoma',
    enrichment_fold=2,
    P_value=0.05
)
```

---

### 2. Visualize Adjacency Difference Heatmap

```python
from stniche import plot_adj_difference_heatmap

plot_adj_difference_heatmap(
    adj_df,
    save_as='bcc_vs_healthy.pdf'
)
```

---

### 3. Extract and Group 3rd-Order Structures

```python
from stniche import extract_and_group_3rd_structures_from_2nd_with_ratio

grouped_structures, sample_wise_structures, df_counts = extract_and_group_3rd_structures_from_2nd_with_ratio(
    adata,
    fdr_filtered_2nd,
    row_key='array_row',
    col_key='array_col',
    sample_key='sample_id',
    group_key='class',
    focus_group='basal cell carcinoma',
    leiden_key='leiden_0.7',
    coverage_threshold=0.8,
    fc_threshold=4.0,
    p_threshold=0.05
)
```

---

### 4. Filter Unique Significant Structures

```python
from stniche import filter_unique_significant_structures

df_unique_filtered = filter_unique_significant_structures(
    df_counts,
    grouped_structures,
    sample_wise_structures,
    adjusted_p_threshold=0.05
)
```

---

### 5. Iteratively Expand and Analyze Structures

```python
from stniche import run_iterative_analysis_over_df

all_iterative_results = run_iterative_analysis_over_df(
    df_unique_filtered,
    grouped_structures,
    sample_wise_structures,
    adata,
    coverage_threshold=0.8,
    sample_key='sample_id',
    row_key='array_row',
    col_key='array_col',
    leiden_key='leiden_0.7',
    focus_group='basal cell carcinoma',
    group_key='class',
    fc_threshold=4.0,
    p_threshold=0.05
)
```

---

### 6. Plot Structure Geometry as Hex View

```python
from stniche import plot_structure_hex

plot_structure_hex(
    all_iterative_results,
    group=3,
    title='Niche Structure',
    hex_size=1,
    save_as='structure_hex.svg'
)
```

---

### 7. Highlight Niche on Visium Image

```python
from stniche import highlight_niche_on_spatial

highlight_niche_on_spatial(
    iter_result=all_iterative_results[3],
    sample_id='Sample123',
    sample_path='/path/to/visium/sample',
    save_path='highlighted_niche.pdf',
    show=True
)
```

---

### 8. Run Niche Differential Expression + Enrichment

```python
from stniche import run_niche_differential_and_enrichment

run_niche_differential_and_enrichment(
    adata,
    all_iterative_results,
    group_index=3,
    padj_thr=0.05,
    lfc_thr=1,
    enrichr_gene_set='KEGG_2021_Human',
    organism='Human',
    top_n=20,
    output_prefix='niche_analysis'
)
```

---

### 9. Visualize Spatial Communication (Sankey)

```python
from stniche import plot_spatial_communication_sankey

plot_spatial_communication_sankey(
    adata,
    all_iterative_results,
    group_index=3,
    direction='niche',
    top_n=50,
    output_dir='output_dir',
    show_plot=True
)
```

---

## 📄 License

MIT License

---

## ✉️ Author

- **Mintian Cui**
- Contact: [1308318910@qq.com](mailto:1308318910@qq.com)
- GitHub: [https://github.com/mintian-cui/stniche](https://github.com/mintian-cui/stniche)
