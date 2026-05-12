#!/usr/bin/env python3
"""
RF Gini Importance Heatmap with ComplexHeatmap-style annotations.
Loads rf_gini data from .RData files and creates a multi-tissue heatmap.

step 1:
  run /cellfile/cellnet/MutationModel/scripts/05_analysis/01_WESdata/01b_modelEvaluationWEX_RF_chrWise_TissueWiseGI_forRebuttal.R
  
step 2: run this script:

python3 ./scripts/05_analysis/01_WESdata/01b_modelEvaluationWEX_RF_chrWise_TissueWiseGI_forRebuttal.py 
--data-dir fig/modelEvaluation/WES_predictor_cors_csv/ --output fig/modelEvaluation/rf_gini_heatmap_chrWise_tissueWise_forRebuttal.png
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from pathlib import Path
import re

# =============================================================================
# CONSTANTS
# =============================================================================

TISSUES = ['brain', 'breast', 'colon', 'esophagus', 'kidney', 
           'liver', 'lung', 'ovary', 'prostate', 'skin']

TISSUE_NAMES = {
    'brain': 'Brain', 'breast': 'Breast', 'colon': 'Colon',
    'esophagus': 'Esophagus', 'kidney': 'Kidney', 'liver': 'Liver',
    'lung': 'Lung', 'ovary': 'Ovary', 'prostate': 'Prostate', 'skin': 'Skin'
}

# Genomic scales (small to big)
SCALE_ORDER = ['1bp', '10bp', '100bp', '1kb', '10kb', '100kb', '1Mb']
SCALE_VALUES = {s: i for i, s in enumerate(SCALE_ORDER)}

# Reversed: 1bp = yellow, 1Mb = purple
SCALE_CMAP = LinearSegmentedColormap.from_list(
    'genomic_scale',
    ['#fde725', '#a0da39', '#5ec962', '#21918c', '#3b528b', '#482878', '#440154'],
    N=len(SCALE_ORDER)
)
NO_SCALE_COLOR = '#d3d3d3'

# Sub-feature colors (detailed coloring for each feature type)
SUB_FEATURE_COLORS = {
    # GC Content
    'GCcontent': '#1f77b4',
    
    # Distance
    'Distance to centromere': '#ff7f0e',
    'Distance to telomere': '#d62728',
    
    # Replication
    'Replication Valleys': '#006400',
    'Replication Peaks': '#000000',
    'Replication WaveSignal': '#2ca02c',
    
    # HiC
    'HiC PCA compartments': '#d62728',
    'HiC interactions': '#ff9896',
    
    # DNA Accessibility - DNAse-seq (purple tones)
    'DNAse-seq signal': '#9467bd',
    'DNAse-seq peaks': '#000000',
    
    # DNA Accessibility - ATAC-seq (lighter purple/pink tones)
    'ATACseq signal': '#c5b0d5',
    'ATACseq peaks': '#000000',
    
    # DNA Methylation
    'DNA methylation': '#8c564b',
    
    # CTCF
    'CTCF signal': '#e377c2',
    'CTCF peaks': '#000000',
    
    # EP300
    'EP300 signal': '#7f7f7f',
    'EP300 peaks': '#000000',
    
    # H3K27ac
    'H3K27ac signal': '#bcbd22',
    'H3K27ac peaks': '#000000',
    
    # H3K27me3
    'H3K27me3 signal': '#17becf',
    'H3K27me3 peaks': '#000000',
    
    # H3K36me3
    'H3K36me3 signal': '#aec7e8',
    'H3K36me3 peaks': '#000000',
    
    # H3K4me1
    'H3K4me1 signal': '#ffbb78',
    'H3K4me1 peaks': '#000000',
    
    # H3K4me3
    'H3K4me3 signal': '#98df8a',
    'H3K4me3 peaks': '#000000',
    
    # H3K9ac
    'H3K9ac signal': '#ff9896',
    'H3K9ac peaks': '#000000',
    
    # H3K9me3
    'H3K9me3 signal': '#c5b0d5',
    'H3K9me3 peaks': '#000000',
    
    # PolR2A
    'PolR2A signal': '#9edae5',
    'PolR2A peaks': '#000000',
    
    # POLR2AphosphoS5
    'POLR2AphosphoS5 signal': '#cedb9c',
    'POLR2AphosphoS5 peaks': '#000000',
    
    # Expression
    'Cancer expression': '#8b0000',
    'Normal tissue expression': '#ff6347',
    'GTEx tissue expression': '#c5b0d5',
    
    # TF Binding
    'Transcription Factor Binding Site Density': '#c49c94',
    'Transcription Factor Binding Site': '#8b6d63',
    'ETS Transcription Factor Binding Site Density': '#bd9e39',
    'ETS Transcription Factor Binding Site': '#654321',
    
    # eQTL
    'GTEx eQTL': '#f7b6d2',
    'GTEx eQTLs -log10(p-value)': '#ff69b4',
    'GTEx eQTLs slope': '#db7093',
    
    # Conservation
    'Conservation phyloP100way': '#c7c7c7',
    
    # Non-B DNA (each type gets a different color)
    'NonB-DNA a-phased repeats': '#e41a1c',
    'NonB-DNA direct repeats': '#377eb8',
    'NonB-DNA g-quadruplex-forming repeats': '#4daf4a',
    'NonB-DNA inverted repeats': '#984ea3',
    'NonB-DNA mirror repeats': '#ff7f00',
    'NonB-DNA short tandem repeats': '#ffff33',
    'NonB-DNA zDNA motifs': '#a65628',
    
    # Coding Effect
    'Coding effect score': '#808080',
}

# Feature group colors (for group-level coloring fallback)
GROUP_COLORS = {
    'GC Content': '#1f77b4',
    'Distance': '#ff7f0e',
    'Replication': '#2ca02c',
    'HiC': '#d62728',
    'DNA Accessibility': '#9467bd',
    'DNA Methylation': '#8c564b',
    'CTCF': '#e377c2',
    'EP300': '#7f7f7f',
    'H3K27ac': '#bcbd22',
    'H3K27me3': '#17becf',
    'H3K36me3': '#aec7e8',
    'H3K4me1': '#ffbb78',
    'H3K4me3': '#98df8a',
    'H3K9ac': '#ff9896',
    'H3K9me3': '#c5b0d5',
    'PolR2A': '#9edae5',
    'POLR2AphosphoS5': '#cedb9c',
    'Expression': '#e377c2',
    'TF Binding': '#c49c94',
    'eQTL': '#f7b6d2',
    'Conservation': '#c7c7c7',
    'Non-B DNA': '#dbdb8d',
    'Coding Effect': '#808080',
    'Other': '#808080'
}

# Feature group patterns
FEATURE_GROUPS = {
    'GC Content': [r'^GCcontent'],
    'Distance': [r'^Distance'],
    'Replication': [r'^Replication'],
    'HiC': [r'^HiC'],
    'DNA Accessibility': [r'^DNAse', r'^ATAC'],
    'DNA Methylation': [r'^DNA methylation'],
    'CTCF': [r'^CTCF'],
    'EP300': [r'^EP300'],
    'H3K27ac': [r'^H3K27ac'],
    'H3K27me3': [r'^H3K27me3'],
    'H3K36me3': [r'^H3K36me3'],
    'H3K4me1': [r'^H3K4me1'],
    'H3K4me3': [r'^H3K4me3'],
    'H3K9ac': [r'^H3K9ac'],
    'H3K9me3': [r'^H3K9me3'],
    'PolR2A': [r'^PolR2A\b'],
    'POLR2AphosphoS5': [r'^POLR2AphosphoS5'],
    'Expression': [r'expression'],
    'TF Binding': [r'Transcription Factor', r'^ETS'],
    'eQTL': [r'eQTL'],
    'Conservation': [r'^Conservation'],
    'Non-B DNA': [r'^NonB-DNA'],
    'Coding Effect': [r'^Coding effect']
}


def extract_scale(feature_name):
    """Extract genomic scale from feature name."""
    for scale in reversed(SCALE_ORDER):  # Check longer patterns first
        if scale in feature_name:
            return scale
    return None


def get_feature_group(feature_name):
    """Determine which group a feature belongs to."""
    for group_name, patterns in FEATURE_GROUPS.items():
        for pattern in patterns:
            if re.search(pattern, feature_name, re.IGNORECASE):
                return group_name
    return 'Other'


def get_sub_feature_color(feature_name):
    """Get color based on sub-feature type."""
    # Try exact match first
    if feature_name in SUB_FEATURE_COLORS:
        return SUB_FEATURE_COLORS[feature_name]
    
    # Try matching without scale suffix
    for sub_feat, color in SUB_FEATURE_COLORS.items():
        if feature_name.startswith(sub_feat):
            return color
    
    # Fallback to group color
    group = get_feature_group(feature_name)
    return GROUP_COLORS.get(group, '#808080')


def parse_feature_annotations(features):
    """Parse feature names into annotations DataFrame."""
    annotations = []
    for feat in features:
        scale = extract_scale(feat)
        group = get_feature_group(feat)
        sub_color = get_sub_feature_color(feat)
        
        annotations.append({
            'feature': feat,
            'scale': scale,
            'group': group,
            'sub_color': sub_color,
        })
    
    return pd.DataFrame(annotations)


def load_rf_gini_data(data_dir, tissues=None):
    """Load rf_gini data from CSV files."""
    if tissues is None:
        tissues = TISSUES
    
    data_dir = Path(data_dir)
    rf_gini_data = {}
    
    for tissue in tissues:
        possible_files = [
            data_dir / f"rf_gini_{tissue}.csv",
            data_dir / f"{tissue}_rf_gini.csv",
            data_dir / f"{tissue}.csv"
        ]
        
        file_path = None
        for pf in possible_files:
            if pf.exists():
                file_path = pf
                break
        
        if file_path is None:
            print(f"Warning: No CSV found for {tissue}, skipping")
            continue
        
        try:
            df = pd.read_csv(file_path, index_col=0)
            rf_gini_data[tissue] = df
            print(f"Loaded: {tissue} ({len(df)} features, {len(df.columns)} chromosomes)")
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
    
    return rf_gini_data


def create_rf_gini_heatmap(gini_data, output_path=None, figsize=(16, 18), dpi=150):
    """
    Create RF Gini importance heatmap matching the original R plot.
    """
    tissues_to_plot = [t for t in TISSUES if t in gini_data]
    n_tissues = len(tissues_to_plot)
    
    if n_tissues == 0:
        raise ValueError("No tissue data to plot")
    
    # Get feature list from first tissue (assuming all same)
    first_tissue = tissues_to_plot[0]
    features = gini_data[first_tissue].index.tolist()
    n_features = len(features)
    
    # Get chromosomes
    chromosomes = gini_data[first_tissue].columns.tolist()
    n_chr = len(chromosomes)
    
    # Parse annotations (keep original order from CSV)
    annotations_df = parse_feature_annotations(features)
    
    # Reverse feature order to match R plot (Coding effect at top)
    features_reversed = features[::-1]
    annotations_df = annotations_df.iloc[::-1].reset_index(drop=True)
    
    # Get group boundaries and ranges
    boundaries = []
    group_ranges = {}
    current_group = None
    start = 0
    
    for i, (_, row) in enumerate(annotations_df.iterrows()):
        if row['group'] != current_group:
            if current_group is not None:
                boundaries.append(i)
                group_ranges[current_group] = (start, i - 1)
            current_group = row['group']
            start = i
    if current_group is not None:
        group_ranges[current_group] = (start, len(annotations_df) - 1)
    
    # Build combined matrix (features x (tissues * chromosomes))
    total_cols = n_chr * n_tissues
    combined_matrix = np.full((n_features, total_cols), np.nan)
    
    col_offset = 0
    tissue_boundaries = [0]
    
    for tissue in tissues_to_plot:
        df = gini_data[tissue]
        # Reorder rows to match reversed features
        df_reordered = df.reindex(features_reversed)
        
        # Scale values per tissue (z-score, then shift to positive)
       # vals = df_reordered.values.astype(float)
        #vals_scaled = (vals - np.nanmean(vals)) / np.nanstd(vals) if np.nanstd(vals) > 0 else vals
        
        # Scale values per tissue (divide by SD only, no centering - matches R's scale(x, center=F))
        vals = df_reordered.values.astype(float)
       # vals_scaled = vals / np.nanstd(vals) if np.nanstd(vals) > 0 else vals
        vals_scaled = vals / np.nanstd(vals) if np.nanstd(vals) > 0 else vals
        
        
        combined_matrix[:, col_offset:col_offset + n_chr] = vals_scaled
        col_offset += n_chr
        tissue_boundaries.append(col_offset)
    
    # Create figure
    fig = plt.figure(figsize=figsize)
    
    # Layout parameters
    hm_left = 0.14
    hm_bottom = 0.08
    hm_width = 0.75
    hm_height = 0.82
    bar_width = 0.012
    gap = 0.003
    
    # Colormap: grey90 to red (matching R scale_fill_gradient)
    #heatmap_cmap = LinearSegmentedColormap.from_list(
   #     'grey_red', ['#e5e5e5', '#fee5d9', '#fcbba1', '#fc9272', '#fb6a4a', '#de2d26', '#a50f15'],
      #  N=256
    # )
     
    heatmap_cmap = LinearSegmentedColormap.from_list(
    'grey_red', ['#e5e5e5', '#ff0000'],  # grey90 to red, matching R
    N=256
    ) 
    
    heatmap_cmap.set_bad(color='#d3d3d3')  # NA values as grey
    
    # Main heatmap
    ax_heatmap = fig.add_axes([hm_left, hm_bottom, hm_width, hm_height])
    
   # vmin, vmax = np.nanpercentile(combined_matrix, [2, 98])
    vmin, vmax = 0, 8
    
    im = ax_heatmap.imshow(
        combined_matrix, cmap=heatmap_cmap, aspect='auto',
        interpolation='nearest', extent=[0, total_cols, n_features, 0],
        vmin=vmin, vmax=vmax
    )
    
    # Tissue separators and labels
    for i, boundary in enumerate(tissue_boundaries[1:-1], 1):
        ax_heatmap.axvline(x=boundary, color='white', linewidth=2)
    
    for i, tissue in enumerate(tissues_to_plot):
        center = (tissue_boundaries[i] + tissue_boundaries[i + 1]) / 2
        ax_heatmap.text(center, -0.5, TISSUE_NAMES.get(tissue, tissue),
                        ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Group separators
    for boundary in boundaries:
        ax_heatmap.axhline(y=boundary, color='white', linewidth=0.5)
    
    ax_heatmap.set_xlim(0, total_cols)
    ax_heatmap.set_ylim(n_features, 0)
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel('Chromosome', fontsize=12, fontweight='bold', labelpad=10)
    ax_heatmap.set_ylabel('', fontsize=12, fontweight='bold', labelpad=60)
    
    # Scale annotation bar
    scale_colors = []
    for _, row in annotations_df.iterrows():
        if row['scale'] is None:
            scale_colors.append(NO_SCALE_COLOR)
        else:
            scale_idx = SCALE_VALUES.get(row['scale'], 0)
            scale_colors.append(SCALE_CMAP(scale_idx / (len(SCALE_ORDER) - 1)))
    
    ax_row_scale = fig.add_axes([hm_left - bar_width - gap, hm_bottom, bar_width, hm_height])
    for i, color in enumerate(scale_colors):
        ax_row_scale.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='none'))
    for boundary in boundaries:
        ax_row_scale.axhline(y=boundary, color='black', linewidth=0.5)
    ax_row_scale.set_xlim(0, 1)
    ax_row_scale.set_ylim(0, n_features)
    ax_row_scale.invert_yaxis()
    ax_row_scale.axis('off')
    
    # Group annotation bar (using sub-feature colors)
    group_colors_list = annotations_df['sub_color'].tolist()
    
    ax_row_group = fig.add_axes([hm_left - 2*bar_width - 2*gap, hm_bottom, bar_width, hm_height])
    for i, color in enumerate(group_colors_list):
        ax_row_group.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='none'))
    for boundary in boundaries:
        ax_row_group.axhline(y=boundary, color='black', linewidth=0.5)
    ax_row_group.set_xlim(0, 1)
    ax_row_group.set_ylim(0, n_features)
    ax_row_group.invert_yaxis()
    ax_row_group.axis('off')
    
    # Group labels (on the left)
    ax_labels = fig.add_axes([0.01, hm_bottom, 0.10, hm_height])
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(0, n_features)
    ax_labels.invert_yaxis()
    for group_name, (start, end) in group_ranges.items():
        center = (start + end) / 2 + 0.5
        ax_labels.text(0.95, center, group_name, ha='right', va='center', fontsize=11, fontweight='bold')
    ax_labels.axis('off')
    
    # Colorbar for heatmap
    ax_cbar = fig.add_axes([hm_left + hm_width + 0.02, hm_bottom + hm_height * 0.6, 0.015, hm_height * 0.3])
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label('Scaled Gini\nImportance', fontsize=10)
    
    # Scale legend (moved further right)
    scale_patches = []
    for i, scale in enumerate(SCALE_ORDER):
        color = SCALE_CMAP(i / (len(SCALE_ORDER) - 1))
        scale_patches.append(mpatches.Patch(color=color, label=scale))
    scale_patches.append(mpatches.Patch(color=NO_SCALE_COLOR, label='No scale'))
    
    fig.legend(handles=scale_patches, title='Genomic Scale', loc='upper right',
               bbox_to_anchor=(0.99, 0.40), fontsize=11, title_fontsize=12)
    
    # Feature sub-type legend (collect unique sub-features)
    # unique_sub_features = {}
    # for _, row in annotations_df.iterrows():
    #     # Get base feature name without scale
    #     feat = row['feature']
    #     for scale in SCALE_ORDER:
    #         if feat.endswith(scale):
    #             feat = feat.replace(f' {scale}', '').strip()
    #             break
    #     if feat not in unique_sub_features:
    #         unique_sub_features[feat] = row['sub_color']
    # 
    # sub_feature_patches = [mpatches.Patch(color=c, label=n) for n, c in unique_sub_features.items()]
    
    
    # Feature sub-type legend (collect unique sub-features, consolidate peaks)
    unique_sub_features = {}
    has_peaks = False
    for _, row in annotations_df.iterrows():
        # Get base feature name without scale
        feat = row['feature']
        for scale in SCALE_ORDER:
            if feat.endswith(scale):
                feat = feat.replace(f' {scale}', '').strip()
                break
        # Skip individual peak entries, track that we have peaks
        if 'peaks' in feat.lower():
            has_peaks = True
            continue
        if feat not in unique_sub_features:
            unique_sub_features[feat] = row['sub_color']
    
    sub_feature_patches = [mpatches.Patch(color=c, label=n) for n, c in unique_sub_features.items()]
    # Add single "Peaks" entry
    if has_peaks:
        sub_feature_patches.append(mpatches.Patch(color='#000000', label='Peaks'))
    
    fig.legend(handles=sub_feature_patches, title='Features', loc='lower center',
               bbox_to_anchor=(0.5, -0.12), fontsize=11, title_fontsize=12, ncol=4)
    
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"Saved: {output_path}")
    else:
        plt.show()
    
    return fig


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate RF Gini importance heatmap')
    parser.add_argument('--data-dir', '-d', required=True, help='Directory with rf_gini_{tissue}.csv files')
    parser.add_argument('--output', '-o', default='./rf_gini_heatmap.png', help='Output path')
    parser.add_argument('--tissues', '-t', nargs='*', help='Tissues to include')
    
    args = parser.parse_args()
    
    # Load data
    gini_data = load_rf_gini_data(args.data_dir, args.tissues)
    
    if not gini_data:
        print("No data loaded. Make sure CSV files exist with pattern: rf_gini_{tissue}.csv")
        return
    
    # Create plot
    create_rf_gini_heatmap(gini_data, args.output)


if __name__ == "__main__":
    main()
