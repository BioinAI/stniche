import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from tqdm import tqdm
from collections import Counter
from itertools import combinations
from matplotlib.patches import RegularPolygon
import scanpy as sc
import os
import re
import gseapy as gp
import squidpy as sq
import plotly.graph_objects as go
import plotly.colors as pc


def compute_groupwise_adjacency_matrix(
    adata,
    row_key='array_row',
    col_key='array_col',
    sample_key='sample_id',
    cluster_key='leiden_0.7',
    group_key='class',
    groups=('GroupA', 'GroupB'),
    focus_group=None,  # 新增：指定关注的 group
    enrichment_fold=4,   # 可选：fold change 倍数阈值
    P_value = 0.05
):
    """
    Compute group-wise adjacency matrices from spatial transcriptomics data
    and identify significantly enriched cell-cell interactions between groups.

    This function calculates cell-type (Leiden cluster) adjacency matrices 
    for each sample, computes average adjacency strengths for two groups (e.g. tumor vs normal), 
    performs Mann-Whitney U tests for differences, and applies FDR correction.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (from Scanpy), with spatial coordinates and cluster assignments.
    row_key : str, default='array_row'
        Column name in `adata.obs` representing row positions of spots.
    col_key : str, default='array_col'
        Column name in `adata.obs` representing column positions of spots.
    sample_key : str, default='sample_id'
        Column name identifying sample IDs in `adata.obs`.
    cluster_key : str, default='leiden_0.7'
        Column name representing cell-type or cluster labels (e.g., Leiden clustering).
    group_key : str, default='class'
        Column in `adata.obs` defining the biological group for comparison (e.g., 'tumor', 'healthy').
    groups : tuple of str, default=('GroupA', 'GroupB')
        A pair of group labels to compare (e.g., ('tumor', 'normal')).
    focus_group : str or None, optional
        Group to focus on when filtering for enriched interactions.
        Used to determine whether a cell-cell interaction is enriched in this group.
    enrichment_fold : float, default=4
        Fold-change threshold to define enrichment of adjacency strength.
        Only used when `focus_group` is provided.
    P_value : float, default=0.05
        FDR-corrected p-value threshold to define statistical significance.

    Returns
    -------
    merged_df : pandas.DataFrame
        Group-wise average adjacency matrices, indexed by group and source cluster,
        with columns as target clusters.
    fdr_filtered : pandas.DataFrame or None
        Table of significantly different adjacency interactions (non-diagonal),
        filtered by FDR < `P_value` and optionally by `enrichment_fold`.
        Returns None if fewer than two groups are provided.

    Raises
    ------
    ValueError
        If `focus_group` is not in the provided `groups`.

    Notes
    -----
    - Assumes 7-neighbor hexagonal spatial topology for adjacency.
    - Each sample contributes one adjacency matrix.
    - Fold-enrichment is calculated only for the `focus_group`.
    """
    spot_data = adata.obs[[row_key, col_key, sample_key, cluster_key]].copy()
    sample_adj_counts = {}

    # 每个样本构建邻接矩阵
    for sample in spot_data[sample_key].unique():
        sample_spots = spot_data[spot_data[sample_key] == sample]
        spot_dict = {
            (row, col): cl for row, col, cl in zip(
                sample_spots[row_key], sample_spots[col_key], sample_spots[cluster_key]
            )
        }

        adj_counts = {}
        for (row, col), cl in spot_dict.items():
            if cl not in adj_counts:
                adj_counts[cl] = {}
            neighbors = [
                (row - 1, col), (row + 1, col),
                (row, col - 2), (row, col + 2),
                (row - 1, col - 1), (row + 1, col - 1),
                (row - 1, col + 1), (row + 1, col + 1)
            ]
            for nb in neighbors:
                if nb in spot_dict:
                    nb_cl = spot_dict[nb]
                    if nb_cl not in adj_counts[cl]:
                        adj_counts[cl][nb_cl] = 0
                    adj_counts[cl][nb_cl] += 1

        adj_df = pd.DataFrame.from_dict(adj_counts, orient='index').fillna(0).astype(float)
        leiden_counts = sample_spots[cluster_key].value_counts()
        for cl in adj_df.index:
            if cl in leiden_counts:
                adj_df.loc[cl] /= leiden_counts[cl]

        sample_adj_counts[sample] = adj_df

    # 合并所有样本邻接矩阵
    all_samples_adj = pd.concat(sample_adj_counts, names=["library_id", "leiden"])
    all_leiden_classes = sorted(all_samples_adj.columns.union(all_samples_adj.index.levels[1]))

    # 获取每组邻接矩阵（保留为 NumPy 用于统计分析）
    group_adj_matrices = {}
    for group in groups:
        samples = adata.obs[adata.obs[group_key] == group][sample_key].unique().tolist()
        group_matrices = []
        for s in samples:
            try:
                df = all_samples_adj.xs(s, level="library_id")
                df = df.fillna(0).reindex(index=all_leiden_classes, columns=all_leiden_classes, fill_value=0)
                group_matrices.append(df.values)
            except KeyError:
                continue
        if group_matrices:
            group_adj_matrices[group] = np.array(group_matrices)

    # 构建合并均值矩阵 DataFrame（用于展示）
    group_avg_flattened = {}
    for group in group_adj_matrices:
        mean_matrix = np.mean(group_adj_matrices[group], axis=0)
        group_avg_flattened[group] = pd.DataFrame(mean_matrix, index=all_leiden_classes, columns=all_leiden_classes)
    merged_df = pd.concat(group_avg_flattened, names=["Group_samll"])

    # 仅支持两个组比较
    if len(groups) != 2:
        print("Only two groups comparison is supported for statistical testing.")
        return merged_df, None

    group1, group2 = groups
    data1 = group_adj_matrices[group1]
    data2 = group_adj_matrices[group2]

    # Mann-Whitney U test + FDR
    p_values, comparisons, group1_means, group2_means = [], [], [], []

    for i, l1 in enumerate(all_leiden_classes):
        for j, l2 in enumerate(all_leiden_classes):
            v1 = data1[:, i, j]
            v2 = data2[:, i, j]

            if np.any(v1) and np.any(v2):
                stat, p = mannwhitneyu(v1, v2, alternative='two-sided')
            else:
                p = 1.0

            p_values.append(p)
            comparisons.append((l1, l2))
            group1_means.append(np.mean(v1))
            group2_means.append(np.mean(v2))

    _, fdrs, _, _ = multipletests(p_values, method='fdr_bh')

    fdr_results = pd.DataFrame({
        "Leiden1": [c[0] for c in comparisons],
        "Leiden2": [c[1] for c in comparisons],
        f"{group1}_mean": group1_means,
        f"{group2}_mean": group2_means,
        "P_value": p_values,
        "FDR_corrected": fdrs
    })

    # 筛选显著差异且非对角（即不同类之间）
    fdr_filtered = fdr_results[
        (fdr_results["FDR_corrected"] < P_value) &
        (fdr_results["Leiden1"] != fdr_results["Leiden2"])
    ]
    
    #fdr_filtered[
    #(fdr_filtered[ f"{group1}_mean"] > 4 * fdr_filtered[f"{group2}_mean"]) &
    #(fdr_filtered[ f"{group1}_mean"] > fdr_filtered[ f"{group1}_mean"].mean())
    #]

    # 若指定了 focus_group，则进一步筛选显著增强的邻接结构
    focus_filtered = None
    if focus_group in groups:
        other_group = [g for g in groups if g != focus_group][0]
        focus_filtered = fdr_filtered[
            (fdr_filtered[f"{focus_group}_mean"] > enrichment_fold * fdr_filtered[f"{other_group}_mean"]) &
            (fdr_filtered[f"{focus_group}_mean"] > fdr_filtered[f"{focus_group}_mean"].mean())
        ]
    else:
        raise ValueError(f"指定组 '{focus_group}' 不在提供的 groups {groups} 中")


    return merged_df, focus_filtered 

def plot_adj_difference_heatmap(
    adj_df,
    cmap='RdBu',
    save_as=None,    # 新增：保存路径 + 格式，例如 'figure.svg'
    show=True
):
    """
    Plot a heatmap of adjacency differences between two groups of cell clusters.

    This function takes a multi-indexed adjacency matrix (group x cluster x cluster),
    computes the difference between the two groups, and visualizes the result as a heatmap.

    Parameters
    ----------
    adj_df : pandas.DataFrame
        Multi-indexed DataFrame where the first level is group label and the 
        rest is a square matrix (e.g., cell-type adjacency averages).
        Typically the output from `compute_groupwise_adjacency_matrix`.

    cmap : str, default='RdBu'
        Colormap used for the heatmap to visualize difference values.
        Choose diverging colormaps (e.g., 'RdBu', 'coolwarm') for better contrast.

    save_as : str or None, optional
        If provided, saves the figure to the specified path (e.g., 'output.svg' or 'figures/heatmap.pdf').
        Format is inferred from the file extension.

    show : bool, default=True
        Whether to display the plot immediately using `plt.show()`.

    Returns
    -------
    None
        Displays and/or saves a matplotlib heatmap showing differences between two group adjacency matrices.

    Notes
    -----
    - Automatically subtracts the second group from the first (e.g., `Group1 - Group2`).
    - Assumes that `adj_df` has exactly two groups in the first-level index.
    - Diagonal entries (self-to-self cluster connections) are set to 0 for clarity.
    """

    # 自动获取 group 名称
    first_group = adj_df.index.levels[0][0]
    second_group = adj_df.index.levels[0][1]

    # 提取 group 的矩阵
    disease_mat = adj_df.loc[first_group]
    healthy_mat = adj_df.loc[second_group]

    # 计算差值矩阵
    diff_mat = disease_mat - healthy_mat

    # 排序 & 数值化
    diff_mat.columns = diff_mat.columns.astype(int)
    diff_mat = diff_mat.reindex(sorted(diff_mat.columns), axis=1)
    diff_mat.index = diff_mat.index.astype(int)
    diff_mat = diff_mat.sort_index(axis=0)

    # 设置对角线为 0
    np.fill_diagonal(diff_mat.values, 0)

    # 绘图
    plt.figure(figsize=(8, 7))
    sns.heatmap(
        diff_mat,
        cmap=cmap,
        linewidths=0.5,
        cbar_kws={'label': f'{first_group} - {second_group}'}
    )
    plt.xlabel('cell cluster')
    plt.ylabel('cell cluster')
    plt.title('Immediate neighborhood')
    plt.tight_layout()

    # 保存图像（如果指定了路径）
    if save_as:
        plt.savefig(save_as, format=save_as.split('.')[-1], dpi=300)

    # 显示或关闭图像
    if show:
        plt.show()
    else:
        plt.close()

# 计算结构的相对几何关系
def geometric_signature(structure):
    coordinates = structure[['x', 'y']].values
    center = np.mean(coordinates, axis=0)
    vectors = coordinates - center
    distances = np.linalg.norm(vectors[:, np.newaxis, :] - vectors[np.newaxis, :, :], axis=-1)

    # 创建距离矩阵
    distance_matrix = pd.DataFrame(distances, index=structure['name'], columns=structure['name'])

    # 关键修正：使用 stack() 但先去掉索引名称，防止与列名冲突
    distance_matrix.index.name = None
    distance_matrix.columns.name = None

    # 转换为 DataFrame
    distance_flat = distance_matrix.stack().reset_index()

    # 确保列名正确
    distance_flat.columns = ['From', 'To', 'Distance']

    # 只保留 From < To，去掉自环(A-A, B-B) 以及重复项(A-B 和 B-A)
    distance_flat = distance_flat[distance_flat['From'] < distance_flat['To']]

    # 按距离排序
    sorted_distances = distance_flat.sort_values(by='Distance').reset_index(drop=True)

    # 获取 label 组合
    labels_dict = structure.set_index('name')['label'].to_dict()
    sorted_label_combinations = [(labels_dict[f], labels_dict[t]) for f, t in zip(sorted_distances['From'], sorted_distances['To'])]

    return sorted_distances, sorted_label_combinations

# 比较签名（分组比较，组内排序后标签一致）
def compare_structures(sig1, sig2, labels1, labels2, tolerance=1e-6):


    # 提取距离列，转换为 NumPy 数组
    distances1 = np.round(sig1['Distance'].to_numpy(), decimals=6)
    distances2 = np.round(sig2['Distance'].to_numpy(), decimals=6)
    
    # 确保两个结构的距离数组长度相同
    if len(distances1) != len(distances2):
        return False

    # 归一化标签顺序：确保每个 (a, b) 组合都是 (min(a,b), max(a,b))
    normalized_labels1 = [tuple(sorted(label)) for label in labels1]
    normalized_labels2 = [tuple(sorted(label)) for label in labels2]
    #print(normalized_labels1)
    #print(normalized_labels2)        
    
    # **按距离分组**
    from collections import defaultdict
    group1 = defaultdict(list)
    group2 = defaultdict(list)

    for d, label in zip(distances1, normalized_labels1):
        group1[d].append(label)

    for d, label in zip(distances2, normalized_labels2):
        group2[d].append(label)
   
    #print(group1)
    #print(group2)    
    # **检查每个分组是否匹配**
    for d in sorted(group1.keys()):  # 确保比较时按距离排序
        if d not in group2:
            return False  # 有不同的距离，结构不匹配
        
        # **排序后比较标签**
        if sorted(group1[d]) != sorted(group2[d]):
            return False  # 组内标签不匹配

    return True

#从二级结构到三级结构
def extract_and_group_3rd_structures_from_2nd_with_ratio(
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
):
    """
    Identify and group third-order niche structures from second-order interactions,
    then perform enrichment filtering based on frequency, coverage, and statistical tests.

    This function builds third-order structures (3-spot motifs) by extending
    significantly enriched 2-spot interactions, and clusters similar structures
    across samples using geometric signatures. It then performs per-structure
    statistical analysis to identify those significantly enriched in a focus group.

    Parameters
    ----------
    adata : AnnData
        Annotated data object containing spatial transcriptomics data,
        including spatial coordinates and clustering information.

    fdr_filtered_2nd : pandas.DataFrame
        DataFrame containing significantly enriched 2-spot cluster interactions
        (typically from `compute_groupwise_adjacency_matrix`), with 'Leiden1' and 'Leiden2' columns.

    row_key : str, default='array_row'
        Column in `adata.obs` representing the row coordinate of each spot.

    col_key : str, default='array_col'
        Column in `adata.obs` representing the column coordinate of each spot.

    sample_key : str, default='sample_id'
        Column in `adata.obs` identifying the sample each spot belongs to.

    group_key : str, default='class1'
        Column in `adata.obs` that defines biological groups (e.g., 'tumor', 'healthy').

    focus_group : str, optional
        The group label (from `group_key`) to be considered the target group for enrichment.

    leiden_key : str, default='leiden_0.7'
        Column in `adata.obs` containing clustering labels (e.g., Leiden clusters).

    coverage_threshold : float, default=0.8
        Minimum proportion of samples in which a structure must appear to be considered valid.

    fc_threshold : float, default=4.0
        Minimum fold change in structure frequency compared to background to be considered enriched.

    p_threshold : float, default=0.05
        FDR-corrected p-value threshold for statistical significance.

    Returns
    -------
    grouped_structures : dict
        A dictionary mapping group IDs to lists of (sample, structure_name) tuples
        representing similar third-order structures across samples.

    sample_wise_structures : dict
        A nested dictionary mapping each sample to its 3rd-order structures and
        their spatial composition.

    df_counts : pandas.DataFrame
        Summary statistics per structure group, including occurrence count, coverage,
        fold change, p-value, and FDR-adjusted p-value.

    Notes
    -----
    - Structures are compared and grouped using geometric signatures (based on relative distances and labels).
    - Mann–Whitney U test is used for enrichment comparison between structure and background.
    - Only structures matching the enrichment criteria are returned in `df_counts`.

    See Also
    --------
    compute_groupwise_adjacency_matrix : For computing significant 2-spot interactions.
    filter_unique_significant_structures : For deduplicating structure groups by label signature.
    """
   
    neighbors_offset = [
        (-1, 0), (1, 0),
        (0, -2), (0, 2),
        (-1, -1), (1, -1),
        (-1, 1), (1, 1)
    ]

    focus_sample_ids = adata.obs[adata.obs[group_key] == focus_group][sample_key].unique()
    sample_wise_structures = {}
    total_sample_count = len(focus_sample_ids)

    
    for sample in focus_sample_ids:
        sample_spots = adata.obs[adata.obs[sample_key] == sample][[row_key, col_key, leiden_key]]
        spot_dict = {(r, c): str(l) for r, c, l in zip(sample_spots[row_key], sample_spots[col_key], sample_spots[leiden_key])}

        visited_2nd = set()
        structures_2nd = {}
        sid_2nd = 0
        for (row, col), spot_leiden in spot_dict.items():
            for dr, dc in neighbors_offset:
                nr, nc = row + dr, col + dc
                if (nr, nc) in spot_dict:
                    pair = {spot_leiden, spot_dict[(nr, nc)]}
                    if tuple(sorted(pair)) in set(tuple(sorted((str(a), str(b)))) for a, b in fdr_filtered_2nd[['Leiden1', 'Leiden2']].values):
                        if frozenset([(row, col), (nr, nc)]) not in visited_2nd:
                            structures_2nd[f'structure_{sid_2nd}'] = [((row, col), spot_leiden), ((nr, nc), spot_dict[(nr, nc)])]
                            visited_2nd.add(frozenset([(row, col), (nr, nc)]))
                            sid_2nd += 1

        visited_3rd = set()
        structures_3rd = {}
        sid_3rd = 0
        for _, pts in structures_2nd.items():
            (p1, l1), (p2, l2) = pts
            for pr, pc in [p1, p2]:
                for dr, dc in neighbors_offset:
                    nr, nc = pr + dr, pc + dc
                    if (nr, nc) in spot_dict and (nr, nc) not in [p1, p2]:
                        key = frozenset([(p1, p2), (nr, nc)])
                        if key not in visited_3rd:
                            structures_3rd[f'structure_{sid_3rd}'] = [(p1, l1), (p2, l2), ((nr, nc), spot_dict[(nr, nc)])]
                            visited_3rd.add(key)
                            sid_3rd += 1

        sample_wise_structures[sample] = {"3rd": structures_3rd}

    structure_signatures = {}
    for sample, data in tqdm(sample_wise_structures.items(), desc="📦 提取结构特征"):
        for structure_name, nodes in data["3rd"].items():
            df = pd.DataFrame(nodes, columns=["(x, y)", "label"])
            df["name"] = df["(x, y)"]
            df["x"] = df["(x, y)"].apply(lambda p: p[0])
            df["y"] = df["(x, y)"].apply(lambda p: p[1])
            df.drop(columns=["(x, y)"], inplace=True)
            sig, labels = geometric_signature(df)
            structure_signatures[(sample, structure_name)] = (sig, labels)

    grouped_structures = {}
    for (sample, sid), (sig, labels) in tqdm(structure_signatures.items(), desc="🧩 匹配结构"):
        matched = None
        for gid, items in grouped_structures.items():
            ref_sample, ref_sid = items[0]
            ref_sig, ref_labels = structure_signatures[(ref_sample, ref_sid)]
            if compare_structures(sig, ref_sig, labels, ref_labels):
                matched = gid
                break
        if matched is not None:
            grouped_structures[matched].append((sample, sid))
        else:
            new_gid = len(grouped_structures)
            grouped_structures[new_gid] = [(sample, sid)]

    print("\n📊 正在统计结构出现频率并比较组间差异...")
    sample_structure_data = []
    all_samples = list(focus_sample_ids)
    sample_total_spots = {s: np.sum(adata.obs[sample_key] == s) for s in all_samples}
    group_ratios = defaultdict(lambda: np.zeros(len(all_samples)))
    sample_index = {s: i for i, s in enumerate(all_samples)}

    for gid, items in grouped_structures.items():
        for sample, _ in items:
            idx = sample_index[sample]
            group_ratios[gid][idx] += 1

    for gid, ratios in group_ratios.items():
        ratio_values = []
        sample_count = 0
        for sid in all_samples:
            idx = sample_index[sid]
            count = ratios[idx]
            total = sample_total_spots[sid]
            ratio = count / total if total > 0 else 0
            ratio_values.append(ratio)
            if count > 0:
                sample_count += 1
        mean_ratio = np.mean(ratio_values)
        sample_structure_data.append([
            gid, sample_count, np.mean(ratios), mean_ratio, ratio_values
        ])

    df_counts = pd.DataFrame(sample_structure_data, columns=[
        "Group_ID", "Count_Sample", "Count_Mean", "Ratio_Mean", "Ratios"
    ])
    
    # 添加 Count_Mean_Other、Fold_Change、Coverage、Significant
    df_counts["Count_Mean_Other"] = [
        np.mean([df_counts.loc[j, "Count_Mean"] for j in range(len(df_counts)) if j != i]) 
        if len(df_counts) > 1 else 1e-9
        for i in range(len(df_counts))
    ]


   # df_counts["Fold_Change"] = df_counts["Count_Mean"] / df_counts["Count_Mean_Other"].replace(0, np.nan)
    df_counts["Fold_Change"] = df_counts["Count_Mean"] / df_counts["Count_Mean_Other"].replace(0, np.nan).fillna(1e-9)

    df_counts["Coverage"] = df_counts["Count_Sample"] / total_sample_count
    
    p_values = []
    for i, row in df_counts.iterrows():
        current_ratios = row["Ratios"]
        other_ratios = np.mean(
            [r for j, r in enumerate(df_counts["Ratios"]) if j != i], axis=0
        )
        try:
            p = mannwhitneyu(current_ratios, other_ratios, alternative="greater").pvalue
        except ValueError:
            p = 1
        p_values.append(p)

    df_counts["P_Value"] = p_values
    df_counts["Adjusted_P"] = multipletests(p_values, method="fdr_bh")[1]
    
    # 新的筛选逻辑
    df_counts["Significant"] = (
        (df_counts["Adjusted_P"] < p_threshold) &
        (df_counts["Fold_Change"] > fc_threshold) &
        (df_counts["Coverage"] >= coverage_threshold)
    )    
    
    
    
    df_counts = df_counts.sort_values(by="Adjusted_P")

    return grouped_structures, sample_wise_structures, df_counts 

#筛选显著的三级结构
def filter_unique_significant_structures(
    df_counts,
    grouped_structures,
    sample_wise_structures,
    adjusted_p_threshold=0.05
):
    """
    Filter and deduplicate significant third-order structure groups based on their label signatures.

    This function identifies unique structure groups by comparing their label compositions
    (converted to hashable signatures), and retains only the most statistically significant
    instance per unique label set. It then filters by adjusted p-value threshold.

    Parameters
    ----------
    df_counts : pandas.DataFrame
        DataFrame containing statistical summaries for each structure group,
        including 'Group_ID', 'P_Value', 'Adjusted_P', etc.

    grouped_structures : dict
        Dictionary mapping group_id to a list of (sample, structure_name) tuples.
        Represents grouped third-order structures across samples.

    sample_wise_structures : dict
        Dictionary mapping sample_id to a nested dictionary with keys like "3rd",
        which stores structure_name → list of (coordinate, label) tuples.

    adjusted_p_threshold : float, default=0.05
        Threshold for filtering structures based on FDR-adjusted p-values.

    Returns
    -------
    df_unique_filtered : pandas.DataFrame
        Filtered DataFrame containing only the most significant structure per
        unique label signature, and passing the p-value threshold.

    Notes
    -----
    - Signature of a structure is defined as a sorted tuple of label counts (from `collections.Counter`).
    - If multiple structure groups share the same label signature, only the one with the lowest p-value is retained.
    """


    def get_structure_label_signature(group_id):
        """获取结构组中的 label signature（Counter，排序并转为 tuple 以便哈希）"""
        structures = grouped_structures.get(group_id, [])
        if not structures:
            return None
        sample, structure_id = structures[0]
        structure = sample_wise_structures.get(sample, {}).get("3rd", {}).get(structure_id, None)
        if not structure:
            return None
        labels = [label for _, label in structure]
        return tuple(sorted(Counter(labels).items()))

    # 收集每个结构标签 signature 对应的最显著一行
    signature_to_best_idx = {}

    for idx, row in df_counts.iterrows():
        group_id = row["Group_ID"]
        signature = get_structure_label_signature(group_id)
        if signature is None:
            continue
        if (
            signature not in signature_to_best_idx or
            row["P_Value"] < df_counts.loc[signature_to_best_idx[signature], "P_Value"]
        ):
            signature_to_best_idx[signature] = idx

    # 获取唯一结构的索引列表
    unique_indices = list(signature_to_best_idx.values())

    # 提取唯一显著结构组并过滤 Adjusted_P
    df_unique_filtered = df_counts.loc[unique_indices].copy()
    df_unique_filtered = df_unique_filtered.sort_values(by="Adjusted_P").reset_index(drop=True)
    df_unique_filtered = df_unique_filtered[df_unique_filtered["Adjusted_P"] <= adjusted_p_threshold]

    return df_unique_filtered

#
# **按标签集合进行初步分组**
def group_by_label_sets(structures_dict):
    label_groups = defaultdict(list)
    
    for struct_name, nodes in structures_dict.items():
        label_set = frozenset([node[1] for node in nodes])  # 只存标签集合
        label_groups[label_set].append(struct_name)

    return label_groups

# **根据几何特征进行最终分组**
def group_structures_by_geometry(structures_dict):
    structure_signatures = {}
    grouped_structures = {}

    # **计算所有结构的几何签名**
    for structure_name, nodes in structures_dict.items():
        df = pd.DataFrame(nodes, columns=["(x, y)", "label"])
        df["name"] = df["(x, y)"]
        df["x"] = df["(x, y)"].apply(lambda p: p[0])
        df["y"] = df["(x, y)"].apply(lambda p: p[1])
        df.drop(columns=["(x, y)"], inplace=True)

        # **计算几何签名**
        signature, labels = geometric_signature(df)

        # **存储签名**
        structure_signatures[structure_name] = (signature, labels)

    # **先按标签集合进行初步分组**
    label_based_groups = group_by_label_sets(structures_dict)

    # **在同一标签组内进行详细几何比较**
    for label_set, structure_names in label_based_groups.items():
        for struct_name in structure_names:
            sig, labels = structure_signatures[struct_name]
            matched_group = None

            for group_key, group_items in grouped_structures.items():
                ref_struct = group_items[0]
                ref_sig, ref_labels = structure_signatures[ref_struct]

                if compare_structures(sig, ref_sig, labels, ref_labels):
                    matched_group = group_key
                    break  # 找到匹配的组就停止

            # **如果匹配到已有分组，就加入该分组**
            if matched_group is not None:
                grouped_structures[matched_group].append(struct_name)
            else:
                # **如果没有匹配的，就创建新分组**
                new_group_key = len(grouped_structures)
                grouped_structures[new_group_key] = [struct_name]

    return grouped_structures

def analyze_structure_groups(grouped_structures, adata, group_key= None, focus_group = None, sample_key=None, fc_threshold=None, p_threshold=None, coverage_threshold=None):
    """
    计算结构组间的差异（不区分 Healthy/Disease），使用每组与其它组进行比较。

    参数：
    - grouped_structures: dict，键为 group_id，值为 (sample, structure_id) 对
    - adata: AnnData 对象
    - sample_key: str，用于标识样本的列名（默认为 'library_id'）
    - fc_threshold: float，fold change 阈值
    - p_threshold: float，p 值显著性阈值
    - coverage_threshold: float，结构覆盖的样本比例阈值

    返回：
    - df_counts: DataFrame，结构组统计结果
    """
    focus_sample_ids = adata.obs[adata.obs[group_key] == focus_group][sample_key].unique()
    #sample_wise_structures = {}
    #total_sample_count = len(focus_sample_ids)
    #all_samples = adata.obs[sample_key].unique()
    all_samples = list(focus_sample_ids)
    
    sample_total_spots = {s: np.sum(adata.obs[sample_key] == s) for s in all_samples}
    total_samples = len(all_samples)    
    #sample_total_spots = {s: np.sum(adata.obs[sample_key] == s) for s in all_samples}
    sample_index = {s: i for i, s in enumerate(all_samples)}

    # 初始化每组结构的每个样本中出现次数
    group_counts = defaultdict(lambda: np.zeros(total_samples))

    for gid, items in grouped_structures.items():
        for sample, _ in items:
            idx = sample_index[sample]
            group_counts[gid][idx] += 1

    # 计算每组在每个样本中的比例
    sample_structure_data = []
    for gid, counts in group_counts.items():
        ratios = []
        sample_count = 0
        for sid in all_samples:
            idx = sample_index[sid]
            total = sample_total_spots[sid]
            ratio = counts[idx] / total if total > 0 else 0
            ratios.append(ratio)
            if counts[idx] > 0:
                sample_count += 1
        mean_ratio = np.mean(ratios)
        mean_count = np.mean(counts)
        sample_structure_data.append([gid, sample_count, mean_count, mean_ratio, ratios])

    df_counts = pd.DataFrame(sample_structure_data, columns=[
        "Group_ID", "Count_Sample", "Count_Mean", "Ratio_Mean", "Ratios"
    ])

    # 添加 Count_Mean_Other、Fold_Change、Coverage
    df_counts["Count_Mean_Other"] = [
        np.mean([df_counts.loc[j, "Count_Mean"] for j in range(len(df_counts)) if j != i]) if len(df_counts) > 1 else 1e-9
        for i in range(len(df_counts))
    ]
    df_counts["Fold_Change"] = df_counts["Count_Mean"] / df_counts["Count_Mean_Other"].replace(0, np.nan).fillna(1e-9)
    df_counts["Coverage"] = df_counts["Count_Sample"] / total_samples

    # 结构组 vs 其它组的 Mann-Whitney U 检验
    p_values = []
    for i, row in df_counts.iterrows():
        current_ratios = row["Ratios"]
        other_ratios = np.mean(
            [df_counts.loc[j, "Ratios"] for j in range(len(df_counts)) if j != i], axis=0
        )
        try:
            p = mannwhitneyu(current_ratios, other_ratios, alternative="greater").pvalue
        except ValueError:
            p = 1.0
        p_values.append(p)

    df_counts["P_Value"] = p_values
    df_counts["Adjusted_P"] = multipletests(p_values, method="fdr_bh")[1]

    # 筛选显著结构
    df_counts["Significant"] = (
        (df_counts["Adjusted_P"] < p_threshold) &
        (df_counts["Fold_Change"] > fc_threshold) &
        (df_counts["Coverage"] >= coverage_threshold)
    )

    return df_counts.sort_values(by="Adjusted_P")

def expand_structure(structure_dict, adata, sample_key=None, row_key=None,
                     col_key=None, leiden_key =None):
    """
    扩展给定的结构到更高一级（在当前结构的基础上增加一个新节点）。
    
    参数：
    - structure_dict: dict，存储当前级别结构的字典，键为 (sample, structure_name)，值为节点列表
    - adata: AnnData 对象，包含样本的空间坐标和类别信息

    返回：
    - expanded_structures: dict，扩展后的结构
    """
    if not structure_dict:
        print("❌ 结构字典为空，无法扩展！")
        return {}
        # **定义邻接关系（六边形网格）**
        
    neighbors_offset = [
        (-1, 0), (1, 0),  # 上下
        (0, -2), (0, 2),  # 左右
        (-1, -1), (1, -1),  # 斜向
        (-1, 1), (1, 1)    # 斜向
    ]
    
    expanded_structures = {}
    visited_expanded = set()
    structure_id = 0

    for (sample, structure), points in structure_dict.items():
        # **获取当前结构的所有点**
        existing_points = {p[0] for p in points}  # 提取 (row, col)
        existing_leidens = {p[1] for p in points}  # 提取 leiden 类型
        current_size = len(points)  # 当前结构的大小

        # **获取该样本的 spot 数据**
        spot_dict = {  # 该样本的 (row, col) -> leiden 映射
            (row, col): str(leiden)
            for row, col, leiden in zip(
                adata.obs.loc[adata.obs[sample_key] == sample, row_key],
                adata.obs.loc[adata.obs[sample_key] == sample, col_key],
                adata.obs.loc[adata.obs[sample_key] == sample, leiden_key]
            )
        }

        # **查找邻居点**
        potential_new_nodes = set()
        for pr, pc in existing_points:
            for dr, dc in neighbors_offset:
                nr, nc = pr + dr, pc + dc

                # **确保新点属于同一 sample，不在已有结构中**
                if (nr, nc) in spot_dict and (nr, nc) not in existing_points:
                    leiden_new = spot_dict[(nr, nc)]
                    potential_new_nodes.add(((nr, nc), leiden_new))

        # **确保至少有一个新点来扩展**
        if not potential_new_nodes:
            continue  # 跳过无法扩展的结构

        # **生成所有可能的新结构**
        for new_node in potential_new_nodes:
            new_structure_key = (sample, f"structure_{structure_id}")

            # **存储扩展后的结构**
            expanded_structures[new_structure_key] = points + [new_node]

            # **防止重复存储**
            visited_expanded.add((sample, tuple(sorted(existing_points | {new_node[0]}))))
            structure_id += 1

    return expanded_structures

def iterative_structure_analysis(initial_structures, adata, reference_row, coverage_threshold=None, 
                                 sample_key=None, row_key=None, col_key=None, 
                                 leiden_key=None,focus_group= None, group_key= None,
                                 fc_threshold= None,p_threshold= None):
    """
    迭代扩展结构，每轮挑出最显著、满足覆盖率的结构，并记录“上一次有效结果”；
    如果下一轮覆盖率不足，则回退到上一次有效结果并结束。
    """
    extracted_structures = initial_structures
    last_valid_structures = initial_structures  # ← 记录上一次通过阈值的结构
    iteration = 0

    while True:
        iteration += 1
        print(f"开始第 {iteration} 轮结构扩展...")

        # 1. 扩展结构
        new_structures = expand_structure(
            extracted_structures, adata, 
            sample_key=sample_key, row_key=row_key,
            col_key=col_key, leiden_key=leiden_key
        )
        if not new_structures:
            print("❌ 无法扩展更多结构。终止迭代。")
            break

        # 2. 分组 & 3. 差异分析
        grouped = group_structures_by_geometry(new_structures)
        df_results = analyze_structure_groups(
            grouped, adata,
            group_key=group_key, focus_group= focus_group,
            sample_key=sample_key,
            fc_threshold=fc_threshold, p_threshold=p_threshold,
            coverage_threshold=coverage_threshold
        )

        # 4. 取出 fold_change >4 的结果
        df_filtered = df_results.query(f"Fold_Change > {fc_threshold}").copy()
        if df_filtered.empty:
            print("没有满足显著性的结构组，终止迭代。")
            break

        # 5. 检查覆盖率
        max_cov = df_filtered["Coverage"].max()
        print(f"本轮最大覆盖率：{max_cov:.3f}")
        if max_cov < coverage_threshold:
            print("覆盖率不足，回退到上一次有效结构，终止迭代。")
            break

        # —— 到这里说明本轮通过了覆盖率阈值 —— 
        # 6. 挑出最优组，并更新 last_valid_structures 和 extracted_structures
        #    （只有在本轮有效时才更新这两个变量）
        # 找到最小 P-value 的候选组
        best_idx = df_filtered["P_Value"].idxmin()
        best_group_id = int(df_filtered.loc[best_idx, "Group_ID"])
        members = grouped[best_group_id]

        # 构建下一轮要扩展的结构集
        next_structures = {
            (sample, struct): new_structures[(sample, struct)]
            for sample, struct in members
        }

        # 更新“上一次有效结果” 和 继续迭代的 extracted_structures
        last_valid_structures = next_structures
        extracted_structures = next_structures

    print(f"✅ 迭代完成，共 {iteration} 轮，返回上一次有效结构集。")
    return last_valid_structures

def run_iterative_analysis_over_df(
    df_filtered_new,
    grouped_structures,
    sample_wise_structures,
    adata,
    coverage_threshold=0.6,
    sample_key='sample_id',
    row_key='array_row',
    col_key='array_col',
    leiden_key='leiden_0.7',
    focus_group=None,
    group_key='class1',
    fc_threshold=4,
    p_threshold=0.05
):
    """
    Run iterative niche analysis for each structure group in the filtered DataFrame.

    For each row in `df_filtered_new` (representing a structure group), this function retrieves
    the corresponding structure instances across samples, performs spatial expansion (iterative search),
    and assesses statistical enrichment. It aggregates the results into a dictionary.

    Parameters
    ----------
    df_filtered_new : pandas.DataFrame
        Filtered DataFrame containing significant structure groups, typically from
        `filter_unique_significant_structures`. Must contain 'Group_ID'.

    grouped_structures : dict
        Mapping from group ID to a list of (sample, structure_name) tuples,
        representing grouped third-order structures.

    sample_wise_structures : dict
        Nested dictionary mapping each sample to a dictionary of 3rd-order structures,
        where each structure is a list of (coordinate, label) tuples.

    adata : AnnData
        Annotated data object containing spot coordinates, cluster assignments,
        and group labels for all samples.

    coverage_threshold : float, default=0.6
        Minimum fraction of samples in which a structure must occur to be considered
        robust during expansion.

    sample_key : str, default='sample_id'
        Column in `adata.obs` identifying each sample.

    row_key : str, default='array_row'
        Column in `adata.obs` representing spot row positions.

    col_key : str, default='array_col'
        Column in `adata.obs` representing spot column positions.

    leiden_key : str, default='leiden_0.7'
        Column in `adata.obs` indicating clustering (e.g., Leiden labels).

    focus_group : str, optional
        Name of the biological group (from `group_key`) to use as the target
        in enrichment comparison.

    group_key : str, default='class1'
        Column in `adata.obs` defining the biological group label (e.g., disease class).

    fc_threshold : float, default=4
        Minimum fold change for enrichment during iterative analysis.

    p_threshold : float, default=0.05
        Adjusted p-value threshold for statistical significance.

    Returns
    -------
    result_dict : dict
        Dictionary where each key is an index from `df_filtered_new` and each value
        is the result from running iterative analysis on the corresponding structure group.

    Notes
    -----
    This function is typically run after identifying significant structures using:
    `compute_groupwise_adjacency_matrix`, `extract_and_group_3rd_structures_from_2nd_with_ratio`,
    and `filter_unique_significant_structures`.

    The output `result_dict` can be used with downstream functions such as:
    `highlight_niche_on_spatial`, `run_niche_differential_and_enrichment`, etc.
    """

    result_dict = {}

    for idx, row in df_filtered_new.iterrows():
        group_id = row["Group_ID"]
        print(idx)
        if group_id not in grouped_structures:
            print(f"⚠️ Group {group_id} 不在 grouped_structures 中，跳过")
            continue

        structure_details = {}
        for sample, structure in grouped_structures[group_id]:
            if sample in sample_wise_structures:
                structures_3rd = sample_wise_structures[sample].get("3rd", {})
                if structure in structures_3rd:
                    structure_details[(sample, structure)] = structures_3rd[structure]
                else:
                    print(f"⚠️ {sample} 中未找到结构 {structure}")
            else:
                print(f"⚠️ 未找到样本 {sample}")

        if not structure_details:
            print(f"⚠️ Group {group_id} 没有有效的结构，跳过")
            continue

        try:
            result = iterative_structure_analysis(structure_details, adata, row, coverage_threshold=coverage_threshold, sample_key=sample_key, row_key=row_key,
                                                  col_key=col_key, leiden_key=leiden_key,focus_group = focus_group, group_key=group_key,
                                                  fc_threshold= fc_threshold,p_threshold= p_threshold)
            result_dict[idx] = result
        except Exception as e:
            print(f"❌ Group {group_id} 出错: {e}")

    return result_dict

def plot_structure_hex(
    all_iterative_results,
    group = 0,
    title="Structure (Hex View)",
    hex_size=1,
    save_as=None,
    show=True
):
    """
    Visualize a spatial niche structure using a hexagonal layout.

    This function renders the spatial organization of a given niche (structure group)
    using axial hexagonal coordinates. Each node is represented as a hexagon and
    colored according to its cluster label.

    Parameters
    ----------
    all_iterative_results : dict
        Dictionary containing iterative structure analysis results, such as the output
        from `run_iterative_analysis_over_df`.

    group : int or str
        The key/index of the structure group to visualize, corresponding to one entry
        in `all_iterative_results`.

    title : str, optional
        Title to display on the plot. Default is "Structure (Hex View)".

    hex_size : float, optional
        Size (radius) of each hexagon in the plot. Default is 1.

    save_as : str or None, optional
        File path to save the plot (including extension such as `.pdf`, `.svg`, etc.).
        If None, the figure is not saved.

    show : bool, optional
        Whether to display the plot interactively. Set to False to suppress display.

    Returns
    -------
    None

    Notes
    -----
    - This function assumes the structure to be visualized is stored in axial hex coordinates
      (q, r) format.
    - Labels are typically cluster IDs used for coloring the hexes.
    - Designed for small-to-medium sized structures; large ones may appear crowded.

    Examples
    --------
    >>> plot_structure_hex(all_iterative_results, group=0, save_as="structure0.svg")
    """

    first_group_key = list(all_iterative_results.keys())[group]
    first_structure_key = list(all_iterative_results[first_group_key].keys())[0]
    structure = all_iterative_results[first_group_key][first_structure_key]
    
    # 六边形的宽和高
    width = np.sqrt(3) * hex_size
    height = 2 * hex_size
    vert_dist = 3/4 * height
    horiz_dist = width

    coords = [pt[0] for pt in structure]
    labels = [int(pt[1]) for pt in structure]

    # 转换轴向坐标为像素坐标
    new_coords = []
    for q, r in coords:
        x = q * horiz_dist
        y = r * vert_dist
        new_coords.append((x, y))

    # 设置图形
    fig, ax = plt.subplots(figsize=(6, 6))

    # 画边（邻接线）
    for i, (x1, y1) in enumerate(new_coords):
        for j, (x2, y2) in enumerate(new_coords):
            if i < j:
                # 邻接判断：六边形邻居有6种位移方向
                dq = abs(coords[i][0] - coords[j][0])
                dr = abs(coords[i][1] - coords[j][1])
                if (dq, dr) in [(1, 0), (0, 2), (1, 1)]:
                    ax.plot([x1, x2], [y1, y2], color='gray', linewidth=1.2, zorder=1)

    # 画六边形节点
    for (x, y), label in zip(new_coords, labels):
        hexagon = RegularPolygon(
            (x, y), numVertices=6, radius=hex_size,
            orientation=np.radians(30), facecolor=plt.cm.tab10(label % 10),
            edgecolor='white', linewidth=1.2, zorder=2
        )
        ax.add_patch(hexagon)
        ax.text(x, y, str(label), ha='center', va='center', color='white', fontsize=12, weight='bold', zorder=3)

    # 设置边界
    xs, ys = zip(*new_coords)
    padding = 2 * hex_size
    ax.set_xlim(min(xs) - padding, max(xs) + padding)
    ax.set_ylim(min(ys) - padding, max(ys) + padding)

    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=14)
    plt.tight_layout()
    
    # 保存图像（如果指定了路径）
    if save_as:
        plt.savefig(save_as, format=save_as.split('.')[-1], dpi=300)

    # 显示或关闭图像
    if show:
        plt.show()
    else:
        plt.close()

def highlight_niche_on_spatial(
    iter_result,
    sample_id,
    sample_path,
    save_path=None,
    show=True,
    spot_size=30,
    niche_color='red',
    background_color='gray',
    image_resolution='hires'
):
    """
    Highlight niche spots on the spatial transcriptomics tissue image.

    This function overlays the inferred niche structure on a Visium tissue image,
    marking niche spots in a specified color while rendering all other spots
    in a muted background.

    Parameters
    ----------
    iter_result : dict
        A single result item from `all_iterative_results`, typically a dictionary
        containing a 'structure' field with coordinates.

    sample_id : str
        Identifier of the sample (e.g., "WSSKNKCLsp12140271"), used to match the tissue image.

    sample_path : str
        Path to the Visium sample folder (should contain the "spatial" subdirectory).

    save_path : str or None, optional
        Full file path to save the output image (e.g., "/path/to/output.pdf").
        If None, the image will not be saved.

    show : bool, default=True
        Whether to display the image interactively.

    spot_size : int, optional
        Size of each spot (dot) plotted. Default is 30.

    niche_color : str, optional
        Color used to highlight niche spots. Default is 'red'.

    background_color : str, optional
        Color used for background (non-niche) spots. Default is 'gray'.

    image_resolution : str, {'hires', 'lowres'}, default='hires'
        Resolution of the tissue image used for background.

    Returns
    -------
    None

    Notes
    -----
    - The function expects the Visium image and coordinates to follow 10x Genomics' format.
    - `iter_result['structure']` must contain a list of spatial (row, col) coordinates.
    - Matplotlib is used for rendering and image export.
    """

    # 提取当前样本的结构信息
    sample_structures = {
        structure: coords_labels
        for (sid, structure), coords_labels in iter_result.items()
        if sid == sample_id
    }

    # 读取样本
    adata_sample = sc.read_visium(sample_path)

    # 构建坐标集合
    niche_all_coords = set()
    for structure, coord_list in sample_structures.items():
        niche_all_coords.update((int(x), int(y)) for (x, y), _ in coord_list)

    # 添加 niche/background 信息
    adata_sample.obs['niche_status'] = adata_sample.obs.apply(
        lambda row: 'niche' if (row['array_row'], row['array_col']) in niche_all_coords else 'background',
        axis=1
    )
    adata_sample.obs['niche_status'] = adata_sample.obs['niche_status'].astype('category')
    adata_sample.obs['niche_status'].cat.set_categories(['background', 'niche'], inplace=True)

    # 设置颜色
    color_palette = [background_color, niche_color]

    # 如果保存路径提供，就设置 figdir 并提取文件名
    save_arg = None
    if save_path:
        figdir, filename = os.path.split(save_path)
        sc.settings.figdir = figdir
        save_arg = filename  # scanpy 会自动拼接 figdir + filename

    # 画图
    sc.pl.spatial(
        adata_sample,
        color='niche_status',
        img_key=image_resolution,
        bw=True,
        palette=color_palette,
        title=f'Niche in {sample_id}',
        spot_size=spot_size,
        show=show,
        save=save_arg
    )

def run_niche_differential_and_enrichment(
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
):
    """
    Perform differential gene expression and functional enrichment analysis for a specific niche group.

    This function compares gene expression between cells within a given niche structure and all
    remaining cells, identifies significantly differentially expressed genes, and performs pathway
    enrichment using Enrichr gene sets. It also optionally saves and plots the results.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix (Scanpy AnnData object), containing gene expression and metadata.

    all_iterative_results : dict
        A dictionary containing the results of iterative niche structure detection.
        Typically generated by `run_iterative_analysis_over_df()`.

    group_index : int, default=0
        Index of the niche structure group to analyze (i.e., the key in `all_iterative_results`).

    padj_thr : float, default=0.05
        Adjusted p-value threshold for calling genes significantly differentially expressed.

    lfc_thr : float, default=1
        Log fold-change threshold for selecting enriched genes (absolute value).

    enrichr_gene_set : str, default='KEGG_2021_Human'
        Name of the gene set collection used for enrichment analysis (via gseapy).

    organism : str, {'Human', 'Mouse'}, default='Human'
        Organism type for enrichment analysis.

    top_n : int, default=20
        Number of top enriched pathways to show and export.

    output_prefix : str or None, optional
        Prefix for output files (e.g., "results/niche_group3"). Files like DEGs, enrichment
        results, and plots will be saved using this prefix. If None, results are not saved.

    show_plot : bool, default=True
        Whether to display enrichment and volcano plots interactively.

    Returns
    -------
    None

    Notes
    -----
    - Requires `scanpy` and `gseapy` to be installed.
    - Output files include: volcano plot, top enriched pathways barplot, DEG list, enrichment results.
    - Designed to help interpret spatial niche structures biologically.
    """

    # 1. 收集 niche spots
    nich_index = set()
    selected_result = all_iterative_results[group_index]
    
    for (sample, _), nodes in selected_result.items():
        for (r, c), _ in nodes:
            hits = adata.obs.reset_index().query(
                "sample_id == @sample and array_row == @r and array_col == @c"
            )['index'].tolist()
            nich_index.update(hits)

    # 2. 找出涉及的所有样本
    samples = {sample for sample, _ in selected_result.keys()}

    # 3. 建立标签
    obs = adata.obs.copy()
    obs['in_nich'] = False
    obs.loc[list(nich_index), 'in_nich'] = True
    mask = obs['sample_id'].isin(samples)
    adata_subset = adata[mask].copy()
    adata_subset.obs['in_nich'] = obs.loc[mask, 'in_nich']

    # 4. 生成标签列
    adata_subset.obs['niche_label'] = adata_subset.obs['in_nich'].map({
        True: 'niche', False: 'background'
    }).astype('category')
    adata_subset.obs['niche_label'] = adata_subset.obs['niche_label'].cat.reorder_categories(
        ['niche', 'background']
    )

    # 5. 差异分析
    sc.tl.rank_genes_groups(
        adata_subset,
        groupby='niche_label',
        groups=['niche'],
        reference='background',
        method='wilcoxon'
    )

    df_de = sc.get.rank_genes_groups_df(adata_subset, group='niche')

    
    # 6. 基因筛选
    sig_genes = df_de[
        (df_de['pvals_adj'] < padj_thr) &
        (df_de['logfoldchanges'] > lfc_thr)
    ].copy().sort_values('logfoldchanges', ascending=False)

    # 7. 富集分析
    gene_list = sig_genes['names'].tolist()
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=enrichr_gene_set,
        organism=organism,
        outdir=None,
        cutoff=0.5
    )

    enrichment_df = enr.results if enr and not enr.results.empty else pd.DataFrame()
    
    if enrichment_df.empty:
        print("❗ 没有富集通路结果，请检查输入参数或基因集是否合理。")
        return sig_genes, None

    # 8. 保存 CSV（如指定了前缀）
    if output_prefix:
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
        df_de.to_csv(f"{output_prefix}/de_results.csv", index=False)
        sig_genes.to_csv(f"{output_prefix}/sig_genes.csv", index=False)
        if not enrichment_df.empty:
            enrichment_df.to_csv(f"{output_prefix}/enrichment.csv", index=False)
    
    
    # 9. 可视化
    kegg_res = enr.results.sort_values('Adjusted P-value').head(top_n).copy()
    kegg_res['Shortened Term'] = kegg_res['Term'].apply(lambda x: re.sub(r"\s*\(.*?\)", "", x))
    kegg_res['-log10(Adjusted P-value)'] = -np.log10(kegg_res['Adjusted P-value'])

    plt.figure(figsize=(12, 8))
    sns.barplot(x='-log10(Adjusted P-value)', y='Shortened Term', data=kegg_res, palette='viridis')
    ax = plt.gca()
    ax.spines['top'].set_color('none')
    ax.spines['right'].set_color('none')
    plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
    plt.ylabel('Pathway', fontsize=12)
    plt.tight_layout()

    if output_prefix:
        pdf_path = f"{output_prefix}/enrichment.pdf"
        plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
        print(f"✅ 图像保存为: {pdf_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()

    return sig_genes, kegg_res

def plot_spatial_communication_sankey(
    adata,
    all_iterative_results,
    group_index=0,
    direction='niche',  # or 'non_niche'
    top_n=50,
    output_dir=None,
    show_plot=True
):
    """
    Visualize spatial communication pathways as a Sankey diagram.

    This function identifies cell-cell communication pairs within or outside niche structures
    and displays the top interactions as a Sankey plot.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix containing spatial transcriptomics data (e.g., from Visium).

    all_iterative_results : dict
        Dictionary containing niche structure detection results, as produced by
        `run_iterative_analysis_over_df()`.

    group_index : int, default=0
        Index of the niche group to analyze (key in `all_iterative_results`).

    direction : {'niche', 'non_niche'}, default='niche'
        Specify whether to extract communications within niche structures ('niche') or
        outside of them ('non_niche').

    top_n : int, default=50
        Number of top communication pairs to include in the Sankey plot (ranked by communication score).

    output_dir : str or None, optional
        Directory to save the Sankey plot and communication table. If provided, files will be saved as:
        - `{output_dir}/sankey.pdf`
        - `{output_dir}/top_communication.csv`

    show_plot : bool, default=True
        Whether to display the Sankey plot interactively.

    Returns
    -------
    None

    Notes
    -----
    - Requires `plotly` and a pre-processed AnnData object with communication scores.
    - Useful for visualizing key ligand-receptor interactions within spatially defined niches.
    """

    # Step 1: 收集所有 niche spot 索引
    selected_result = all_iterative_results[group_index]
    nich_index = set()
    for (sample, _), nodes in selected_result.items():
        for (r, c), _ in nodes:
            hits = adata.obs.reset_index().query(
                "sample_id == @sample and array_row == @r and array_col == @c"
            )['index'].tolist()
            nich_index.update(hits)

    # Step 2: 构建 adata_combined_all，仅包含相关样本
    samples = {sample for sample, _ in selected_result.keys()}
    obs = adata.obs.copy()
    obs['in_nich'] = False
    obs.loc[list(nich_index), 'in_nich'] = True
    mask = obs['sample_id'].isin(samples)
    adata_combined_all = adata[mask].copy()
    adata_combined_all.obs['in_nich'] = obs.loc[mask, 'in_nich']
    adata_combined_all.obsm['spatial'] = adata_combined_all.obs[['array_row', 'array_col']].to_numpy(dtype=float)
    adata_combined_all.obs['niche_group'] = adata_combined_all.obs['in_nich'].map(lambda x: 'niche' if x else 'non_niche').astype('category')

    # Step 3: 建立空间邻接图
    sq.gr.spatial_neighbors(adata_combined_all, coord_type="grid")

    # Step 4: 空间通讯计算
    niche_result = sq.gr.ligrec(
        adata_combined_all,
        cluster_key='niche_group',
        n_perms=100,
        copy=True
    )

    # Step 5: 提取 means
    means = niche_result['means']
    if direction == 'niche':
        comm_scores = means[('niche', 'non_niche')]
    elif direction == 'non_niche':
        comm_scores = means[('non_niche', 'niche')]
    else:
        raise ValueError("❌ direction 参数必须是 'niche' 或 'non_niche'")

    # Step 6: 排序并取前 N
    top_pairs = comm_scores.sort_values(ascending=False).head(top_n)
    top_df = top_pairs.reset_index()
    top_df.columns = ['Ligand', 'Receptor', 'Score']

    # Step 7: 构造 Sankey 图
    ligands = top_df['Ligand'].tolist()
    receptors = top_df['Receptor'].tolist()
    scores = top_df['Score'].tolist()

    nodes = [{"name": name} for name in set(ligands) | set(receptors)]
    all_nodes = [node['name'] for node in nodes]
    links = [{
        "source": all_nodes.index(lig),
        "target": all_nodes.index(rec),
        "value": round(score, 4)
    } for lig, rec, score in zip(ligands, receptors, scores)]

    # 颜色渐变（可选）
    colorscale = pc.sample_colorscale('Reds', [i/len(scores) for i in range(len(scores))])

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=["#FF6666" if n in ligands else "#66B2FF" for n in all_nodes]
        ),
        link=dict(
            source=[link['source'] for link in links],
            target=[link['target'] for link in links],
            value=[link['value'] for link in links],
            color=colorscale
        )
    )])

    fig.update_layout(
        title_text=f"Top {top_n} Communication: {direction.replace('_', ' ')}",
        font_size=12
    )

    # Step 8: 保存 PDF 和 CSV
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        pdf_path = os.path.join(output_dir, "sankey.pdf")
        csv_path = os.path.join(output_dir, "top_communication.csv")
        fig.write_image(pdf_path, format='pdf', scale=2)
        top_df.to_csv(csv_path, index=False)
        print(f"✅ Sankey 图已保存: {pdf_path}")
        print(f"✅ 通讯得分表已保存: {csv_path}")

    if show_plot:
        fig.show()

    return top_df




