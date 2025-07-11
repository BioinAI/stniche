
# stNiche

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
pip install stniche
```

---

## 🚀 Usage Overview

### 1. Compute Adjacency Differences

```python
import stNiche

adj_df, fdr_filtered_2nd = stNiche.compute_groupwise_adjacency_matrix(
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
stNiche.plot_adj_difference_heatmap(
    adj_df,
    cmap='RdBu',
    save_as=None,    
    show=True
)
```

---

### 3. Extract and Group 3rd-Order Structures

```python
grouped_structures, sample_wise_structures, df_counts = stNiche.extract_and_group_3rd_structures_from_2nd_with_ratio(
    adata,
    fdr_filtered_2nd,
    row_key='array_row',
    col_key='array_col',
    sample_key='sample_id',
    group_key='class1',
    focus_group=None,
    leiden_key='leiden_0.7',
    coverage_threshold = 0.8,    
    fc_threshold=4.0,
    p_threshold=0.05
)
```

---

### 4. Filter Unique Significant Structures

```python

df_unique_filtered = stNiche.filter_unique_significant_structures(
    df_counts,
    grouped_structures,
    sample_wise_structures,
    adjusted_p_threshold=0.05
)
```

---

### 5. Iteratively Expand and Analyze Structures

```python
all_iterative_results = stNiche.run_iterative_analysis_over_df(
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
stNiche.plot_structure_hex(
    all_iterative_results,
    group = 0,
    title="Structure (Hex View)",
    hex_size=1,
    save_as=None,
    show=True
)
```

---

### 7. Highlight Niche on Visium Image

```python
stNiche.highlight_niche_on_spatial(
    iter_result=all_iterative_results[3],
    sample_id='Sample123',
    sample_path='/path/to/visium/sample',
    save_path='highlighted_niche.pdf',
    show=True,
    spot_size=30,
    niche_color='red',
    background_color='gray',
    image_resolution='hires'
)
```

---

### 8. Run Niche Differential Expression + Enrichment

```python
stNiche.run_niche_differential_and_enrichment(
    adata,
    all_iterative_results,
    group_index=0,
    padj_thr=0.05,
    lfc_thr=1,
    enrichr_gene_set='KEGG_2021_Human',
    organism='Human',
    top_n=20,
    output_prefix=None,
    show_plot=True
)
```

---

### 9. Visualize Spatial Communication (Sankey)

```python
stNiche.plot_spatial_communication_sankey(
    adata,
    all_iterative_results,
    group_index=0,
    direction='niche',  # or 'non_niche'
    top_n=50,
    output_dir=None,
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
- GitHub: [https://github.com/BioinAI/stniche](https://github.com/BioinAI/stniche)
