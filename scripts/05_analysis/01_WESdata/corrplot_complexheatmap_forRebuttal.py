#!/usr/bin/env python3
"""

step 1: run script in R to extract correlation matrices per tissue: /cellfile/cellnet/MutationModel/scripts/05_analysis/01_WESdata/01b_modelEvaluationWEX_predictorCorrPlots_forRebuttal.R
step 2:
Usage: python3 /cellfile/cellnet/MutationModel/scripts/05_analysis/01_WESdata/corrplot_complexheatmap_forRebuttal.py fig/modelEvaluation/WES_predictor_cors_csv/cors_skin.csv --tissue 'Skin'
 --output fig/modelEvaluation/
ComplexHeatmap-style correlation plot for genomic features.
Handles both original feature names AND p2P-mapped pretty names.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import re

# =============================================================================
# CONSTANTS
# =============================================================================

SCALE_ORDER = ['1Mb', '100kb', '10kb', '1kb', '100bp', '10bp', '1bp']
SCALE_VALUES = {'1bp': 0, '10bp': 1, '100bp': 2, '1kb': 3, '10kb': 4, '100kb': 5, '1Mb': 6}

SCALE_CMAP = LinearSegmentedColormap.from_list(
    'genomic_scale', 
    ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725'],
    N=len(SCALE_ORDER)
)
NO_SCALE_COLOR = '#d3d3d3'
PEAKS_COLOR = '#000000'

# =============================================================================
# FEATURE DEFINITIONS - Supports both original and p2P-mapped names
# =============================================================================

FEATURE_GROUPS = {
    'GC Content': {
        'patterns': [r'GCcontent'],
        'color': '#1f77b4',
        'sub_features': {
            'GC Content': {'patterns': [r'GCcontent'], 'color': '#1f77b4'}
        }
    },
    'Distance': {
        'patterns': [r'Distance to', r'dist_'],
        'color': '#ff7f0e',
        'sub_features': {
            'Centromere Distance': {'patterns': [r'centromere'], 'color': '#ff7f0e'},
            'Telomere Distance': {'patterns': [r'telomere'], 'color': '#d62728'}
        }
    },
    'Replication': {
        'patterns': [r'^Replication'],
        'color': '#2ca02c',
        'sub_features': {
            'Replication Timing': {'patterns': [r'WaveSignal', r'timing'], 'color': '#2ca02c'},
            'Replication Peaks': {'patterns': [r'Replication Peaks', r'Replication_peaks'], 'color': '#98df8a'},
            'Replication Valleys': {'patterns': [r'Valleys'], 'color': '#006400'}
        }
    },
    'HiC': {
        'patterns': [r'^HiC'],
        'color': '#d62728',
        'sub_features': {
            'HiC Compartment': {'patterns': [r'compartment', r'PCA'], 'color': '#d62728'},
            'HiC Interactions': {'patterns': [r'interactions'], 'color': '#ff9896'}
        }
    },
    'DNA Accessibility': {
        'patterns': [r'DNAse', r'ATAC'],
        'color': '#9467bd',
        'sub_features': {
            'DNAse Signal': {'patterns': [r'DNAse.*signal'], 'color': '#9467bd'},
            'DNAse Peaks': {'patterns': [r'DNAse.*peaks'], 'color': PEAKS_COLOR},
            'ATAC-seq Signal': {'patterns': [r'ATAC.*signal'], 'color': '#c5b0d5'},
            'ATAC-seq Peaks': {'patterns': [r'ATAC.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'DNA Methylation': {
        'patterns': [r'methylation'],
        'color': '#8c564b',
        'sub_features': {
            'DNA Methylation': {'patterns': [r'methylation'], 'color': '#8c564b'}
        }
    },
    'CTCF': {
        'patterns': [r'^CTCF'],
        'color': '#e377c2',
        'sub_features': {
            'CTCF Signal': {'patterns': [r'CTCF.*signal'], 'color': '#e377c2'},
            'CTCF Peaks': {'patterns': [r'CTCF.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'EP300': {
        'patterns': [r'^EP300'],
        'color': '#7f7f7f',
        'sub_features': {
            'EP300 Signal': {'patterns': [r'EP300.*signal'], 'color': '#7f7f7f'},
            'EP300 Peaks': {'patterns': [r'EP300.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K27ac': {
        'patterns': [r'^H3K27ac'],
        'color': '#bcbd22',
        'sub_features': {
            'H3K27ac Signal': {'patterns': [r'H3K27ac.*signal'], 'color': '#bcbd22'},
            'H3K27ac Peaks': {'patterns': [r'H3K27ac.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K27me3': {
        'patterns': [r'^H3K27me3'],
        'color': '#17becf',
        'sub_features': {
            'H3K27me3 Signal': {'patterns': [r'H3K27me3.*signal'], 'color': '#17becf'},
            'H3K27me3 Peaks': {'patterns': [r'H3K27me3.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K36me3': {
        'patterns': [r'^H3K36me3'],
        'color': '#aec7e8',
        'sub_features': {
            'H3K36me3 Signal': {'patterns': [r'H3K36me3.*signal'], 'color': '#aec7e8'},
            'H3K36me3 Peaks': {'patterns': [r'H3K36me3.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K4me1': {
        'patterns': [r'^H3K4me1'],
        'color': '#ffbb78',
        'sub_features': {
            'H3K4me1 Signal': {'patterns': [r'H3K4me1.*signal'], 'color': '#ffbb78'},
            'H3K4me1 Peaks': {'patterns': [r'H3K4me1.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K4me3': {
        'patterns': [r'^H3K4me3'],
        'color': '#98df8a',
        'sub_features': {
            'H3K4me3 Signal': {'patterns': [r'H3K4me3.*signal'], 'color': '#98df8a'},
            'H3K4me3 Peaks': {'patterns': [r'H3K4me3.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K9ac': {
        'patterns': [r'^H3K9ac'],
        'color': '#ff9896',
        'sub_features': {
            'H3K9ac Signal': {'patterns': [r'H3K9ac.*signal'], 'color': '#ff9896'},
            'H3K9ac Peaks': {'patterns': [r'H3K9ac.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'H3K9me3': {
        'patterns': [r'^H3K9me3'],
        'color': '#c5b0d5',
        'sub_features': {
            'H3K9me3 Signal': {'patterns': [r'H3K9me3.*signal'], 'color': '#c5b0d5'},
            'H3K9me3 Peaks': {'patterns': [r'H3K9me3.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'PolR2A': {
        'patterns': [r'^PolR2A\b'],
        'color': '#9edae5',
        'sub_features': {
            'PolR2A Signal': {'patterns': [r'PolR2A.*signal'], 'color': '#9edae5'},
            'PolR2A Peaks': {'patterns': [r'PolR2A.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'POLR2AphosphoS5': {
        'patterns': [r'^POLR2AphosphoS5'],
        'color': '#cedb9c',
        'sub_features': {
            'POLR2AphosphoS5 Signal': {'patterns': [r'POLR2AphosphoS5.*signal'], 'color': '#cedb9c'},
            'POLR2AphosphoS5 Peaks': {'patterns': [r'POLR2AphosphoS5.*peaks'], 'color': PEAKS_COLOR}
        }
    },
    'Expression': {
        'patterns': [r'expression'],
        'color': '#c5b0d5',
        'sub_features': {
            'GTEx Expression': {'patterns': [r'GTEx.*expression', r'healthy'], 'color': '#c5b0d5'},
            'Cancer Expression': {'patterns': [r'Cancer.*expression', r'cancer'], 'color': '#8b0000'}
        }
    },
    'TF Binding': {
        'patterns': [r'Transcription Factor', r'^TF_binding', r'^ETS'],
        'color': '#c49c94',
        'sub_features': {
            'TF Binding Sites': {'patterns': [r'^Transcription Factor Binding Site$', r'^TF_binding$'], 'color': '#c49c94'},
            'TF Binding Density': {'patterns': [r'Transcription Factor.*Density'], 'color': '#8c6d31'},
            'ETS TF Binding Sites': {'patterns': [r'^ETS Transcription Factor Binding Site$', r'^ETS_TF'], 'color': '#654321'},
            'ETS TF Binding Density': {'patterns': [r'ETS.*Density'], 'color': '#bd9e39'}
        }
    },
    'eQTL': {
        'patterns': [r'eQTL'],
        'color': '#f7b6d2',
        'sub_features': {
            'eQTL': {'patterns': [r'^GTEx eQTL$', r'^eQTL$'], 'color': '#f7b6d2'},
            'eQTL p-value': {'patterns': [r'p-value', r'pval'], 'color': '#ff69b4'},
            'eQTL slope': {'patterns': [r'slope'], 'color': '#db7093'}
        }
    },
    'Conservation': {
        'patterns': [r'Conservation', r'conservation', r'phyloP'],
        'color': '#c7c7c7',
        'sub_features': {
            'Conservation': {'patterns': [r'Conservation', r'conservation', r'phyloP'], 'color': '#c7c7c7'}
        }
    },
    'Non-B DNA': {
        'patterns': [r'NonB-DNA', r'^aPhased', r'^direct_repeats', r'^g_quadruplex', 
                     r'^inverted_repeats', r'^mirror_repeats', r'^short_tandem', r'^zDNA'],
        'color': '#dbdb8d',
        'sub_features': {
            'A-phased Repeats': {'patterns': [r'a-phased', r'aPhased'], 'color': '#dbdb8d'},
            'Direct Repeats': {'patterns': [r'direct'], 'color': '#b5cf6b'},
            'G-quadruplex': {'patterns': [r'quadruplex'], 'color': '#e7cb94'},
            'Inverted Repeats': {'patterns': [r'inverted'], 'color': '#ce6dbd'},
            'Mirror Repeats': {'patterns': [r'mirror'], 'color': '#9c9ede'},
            'Short Tandem Repeats': {'patterns': [r'short tandem', r'short_tandem'], 'color': '#d6616b'},
            'Z-DNA': {'patterns': [r'zDNA'], 'color': '#e7969c'}
        }
    },
    'Coding Effect': {
        'patterns': [r'Coding effect', r'^effect$'],
        'color': '#808080',
        'sub_features': {
            'Coding Effect': {'patterns': [r'Coding effect', r'^effect$'], 'color': '#808080'}
        }
    }
}


def extract_scale(feature_name):
    """Extract genomic scale from feature name."""
    for scale in SCALE_ORDER:
        if scale in feature_name:
            return scale
    return None


def get_feature_group(feature_name):
    """Determine which group a feature belongs to."""
    for group_name, group_info in FEATURE_GROUPS.items():
        for pattern in group_info['patterns']:
            if re.search(pattern, feature_name, re.IGNORECASE):
                return group_name
    return 'Other'


def get_sub_feature(feature_name, group):
    """Determine sub-feature and color."""
    if group not in FEATURE_GROUPS:
        return feature_name, '#808080'
    
    group_info = FEATURE_GROUPS[group]
    for sub_name, sub_info in group_info.get('sub_features', {}).items():
        for pattern in sub_info['patterns']:
            if re.search(pattern, feature_name, re.IGNORECASE):
                return sub_name, sub_info['color']
    
    return group, group_info.get('color', '#808080')


def parse_feature_annotations(features):
    """Parse feature names into annotations DataFrame."""
    annotations = []
    for feat in features:
        scale = extract_scale(feat)
        group = get_feature_group(feat)
        sub_feature, sub_color = get_sub_feature(feat, group)
        
        annotations.append({
            'feature': feat,
            'scale': scale,
            'group': group,
            'sub_feature': sub_feature,
            'sub_color': sub_color
        })
    
    return pd.DataFrame(annotations)


def get_group_boundaries(annotations_df):
    """Get indices where feature groups change."""
    boundaries = []
    current_group = None
    for i, (_, row) in enumerate(annotations_df.iterrows()):
        if row['group'] != current_group:
            if current_group is not None:
                boundaries.append(i)
            current_group = row['group']
    return boundaries


def get_group_ranges(annotations_df):
    """Get start and end indices for each group."""
    ranges = {}
    current_group = None
    start = 0
    
    for i, (_, row) in enumerate(annotations_df.iterrows()):
        if row['group'] != current_group:
            if current_group is not None:
                ranges[current_group] = (start, i - 1)
            current_group = row['group']
            start = i
    
    if current_group is not None:
        ranges[current_group] = (start, len(annotations_df) - 1)
    
    return ranges


def create_corrplot(cors_matrix, annotations_df, tissue_name='', 
                    output_path=None, figsize=(14, 14), dpi=150):
    """Create ComplexHeatmap-style correlation plot."""
    n = len(cors_matrix)
    
    fig = plt.figure(figsize=figsize)
    
    # Get colors
    scale_colors = []
    for _, row in annotations_df.iterrows():
        if row['scale'] is None:
            scale_colors.append(NO_SCALE_COLOR)
        else:
            scale_idx = SCALE_VALUES.get(row['scale'], 0)
            scale_colors.append(SCALE_CMAP(scale_idx / (len(SCALE_ORDER) - 1)))
    
    group_colors = annotations_df['sub_color'].tolist()
    group_boundaries = get_group_boundaries(annotations_df)
    group_ranges = get_group_ranges(annotations_df)
    
    # PuOr colormap
    corr_cmap = LinearSegmentedColormap.from_list(
        'PuOr_custom',
        ['#2d004b', '#542788', '#8073ac', '#b2abd2', '#d8daeb', '#f7f7f7',
         '#fee0b6', '#fdb863', '#e08214', '#b35806', '#7f3b08'],
        N=256
    )
    
    # Layout
    hm_left = 0.15
    hm_bottom = 0.15
    hm_size = 0.65
    bar_width = 0.012
    gap = 0.003
    
    # Main heatmap
    ax_heatmap = fig.add_axes([hm_left, hm_bottom, hm_size, hm_size])
    
    mask = np.triu(np.ones_like(cors_matrix, dtype=bool), k=1)
    cors_masked = np.ma.masked_array(cors_matrix.values, mask=mask)
    
    im = ax_heatmap.imshow(
        cors_masked, cmap=corr_cmap, vmin=-1, vmax=1, 
        aspect='auto', interpolation='nearest', extent=[0, n, n, 0]
    )
    
    for boundary in group_boundaries:
        ax_heatmap.axhline(y=boundary, color='white', linewidth=2)
        ax_heatmap.axvline(x=boundary, color='white', linewidth=2)
    
    ax_heatmap.set_xlim(0, n)
    ax_heatmap.set_ylim(n, 0)
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    
    # Row annotations
    ax_row_scale = fig.add_axes([hm_left - bar_width - gap, hm_bottom, bar_width, hm_size])
    for i, color in enumerate(scale_colors):
        ax_row_scale.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='white', linewidth=0.3))
    for boundary in group_boundaries:
        ax_row_scale.axhline(y=boundary, color='black', linewidth=1.5)
    ax_row_scale.set_xlim(0, 1)
    ax_row_scale.set_ylim(0, n)
    ax_row_scale.invert_yaxis()
    ax_row_scale.axis('off')
    
    ax_row_group = fig.add_axes([hm_left - 2*bar_width - 2*gap, hm_bottom, bar_width, hm_size])
    for i, color in enumerate(group_colors):
        ax_row_group.add_patch(plt.Rectangle((0, i), 1, 1, facecolor=color, edgecolor='white', linewidth=0.3))
    for boundary in group_boundaries:
        ax_row_group.axhline(y=boundary, color='black', linewidth=1.5)
    ax_row_group.set_xlim(0, 1)
    ax_row_group.set_ylim(0, n)
    ax_row_group.invert_yaxis()
    ax_row_group.axis('off')
    
    ax_row_labels = fig.add_axes([0.01, hm_bottom, 0.11, hm_size])
    ax_row_labels.set_xlim(0, 1)
    ax_row_labels.set_ylim(0, n)
    ax_row_labels.invert_yaxis()
    for group_name, (start, end) in group_ranges.items():
        center = (start + end) / 2 + 0.5
        ax_row_labels.text(0.95, center, group_name, ha='right', va='center', fontsize=11, fontweight='bold')
    ax_row_labels.axis('off')
    
    # Column annotations
    ax_col_scale = fig.add_axes([hm_left, hm_bottom + hm_size + gap, hm_size, bar_width])
    for i, color in enumerate(scale_colors):
        ax_col_scale.add_patch(plt.Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='white', linewidth=0.3))
    for boundary in group_boundaries:
        ax_col_scale.axvline(x=boundary, color='black', linewidth=1.5)
    ax_col_scale.set_xlim(0, n)
    ax_col_scale.set_ylim(0, 1)
    ax_col_scale.axis('off')
    
    ax_col_group = fig.add_axes([hm_left, hm_bottom + hm_size + bar_width + 2*gap, hm_size, bar_width])
    for i, color in enumerate(group_colors):
        ax_col_group.add_patch(plt.Rectangle((i, 0), 1, 1, facecolor=color, edgecolor='white', linewidth=0.3))
    for boundary in group_boundaries:
        ax_col_group.axvline(x=boundary, color='black', linewidth=1.5)
    ax_col_group.set_xlim(0, n)
    ax_col_group.set_ylim(0, 1)
    ax_col_group.axis('off')
    
    ax_col_labels = fig.add_axes([hm_left, hm_bottom + hm_size + 2*bar_width + 3*gap, hm_size, 0.12])
    ax_col_labels.set_xlim(0, n)
    ax_col_labels.set_ylim(0, 1)
    for group_name, (start, end) in group_ranges.items():
        center = (start + end) / 2 + 0.5
        ax_col_labels.text(center, 0.05, group_name, ha='left', va='bottom', fontsize=11, fontweight='bold', rotation=45)
    ax_col_labels.axis('off')
    
    # Colorbar
    ax_cbar = fig.add_axes([hm_left + hm_size + 0.02, hm_bottom, 0.02, hm_size])
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label('Correlation', fontsize=11)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    
    # Legends
    sub_feature_colors = {}
    has_peaks = False
    for _, row in annotations_df.iterrows():
        if 'peaks' in row['sub_feature'].lower() or 'Peaks' in row['sub_feature']:
            has_peaks = True
        else:
            if row['sub_feature'] not in sub_feature_colors:
                sub_feature_colors[row['sub_feature']] = row['sub_color']
    
    sub_feature_patches = [mpatches.Patch(color=c, label=n) for n, c in sub_feature_colors.items()]
    if has_peaks:
        sub_feature_patches.append(mpatches.Patch(color=PEAKS_COLOR, label='Peaks'))
    
    scale_patches = []
    for i, scale in enumerate(reversed(SCALE_ORDER)):
        color = SCALE_CMAP((len(SCALE_ORDER) - 1 - i) / (len(SCALE_ORDER) - 1))
        scale_patches.append(mpatches.Patch(color=color, label=scale))
    scale_patches.append(mpatches.Patch(color=NO_SCALE_COLOR, label='No scale'))
    
    fig.legend(handles=sub_feature_patches, title='Features', loc='upper left',
               bbox_to_anchor=(0.01, 0.13), fontsize=9, title_fontsize=10, ncol=4)
    fig.legend(handles=scale_patches, title='Genomic Scale', loc='upper right',
               bbox_to_anchor=(0.75, 0.13), fontsize=9, title_fontsize=10)
    
    #if tissue_name:
     #   fig.suptitle(f'Feature Correlation Matrix - {tissue_name}', fontsize=14, fontweight='bold', y=0.98)
    
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"Saved: {output_path}")
    else:
        plt.show()
    
    return fig


def process_tissue(csv_path, tissue_name='', output_dir=None, exclude_features=None):
    """Main function to process a correlation matrix CSV."""
    cors = pd.read_csv(csv_path, index_col=0)
    
    if exclude_features:
        cols = [c for c in cors.columns if c not in exclude_features]
        cors = cors.loc[cols, cols]
    
    features = cors.columns.tolist()
    annotations_df = parse_feature_annotations(features)
    
    # Sort by group
    group_order = list(FEATURE_GROUPS.keys())
    annotations_df['group_order'] = annotations_df['group'].apply(
        lambda x: group_order.index(x) if x in group_order else len(group_order)
    )
    annotations_df['scale_order'] = annotations_df['scale'].apply(
        lambda x: -SCALE_VALUES.get(x, -1) if x else 999
    )
    annotations_df = annotations_df.sort_values(['group_order', 'sub_feature', 'scale_order'])
    
    sorted_features = annotations_df['feature'].tolist()
    cors = cors.loc[sorted_features, sorted_features]
    annotations_df = annotations_df.reset_index(drop=True)
    
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / f"predictor_corrplot_{tissue_name.lower().replace(' ', '_')}_forRebuttal.png"
    else:
        output_path = None
    
    fig = create_corrplot(cors, annotations_df, tissue_name, output_path)
    return fig, cors, annotations_df


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate correlation plots')
    parser.add_argument('csv_path', help='Path to correlation matrix CSV')
    parser.add_argument('--tissue', '-t', default='', help='Tissue name')
    parser.add_argument('--output', '-o', default='./output', help='Output directory')
    parser.add_argument('--exclude', '-e', nargs='*', default=None, help='Features to exclude')
    
    args = parser.parse_args()
    
    process_tissue(
        csv_path=args.csv_path,
        tissue_name=args.tissue,
        output_dir=args.output,
        exclude_features=args.exclude
    )
