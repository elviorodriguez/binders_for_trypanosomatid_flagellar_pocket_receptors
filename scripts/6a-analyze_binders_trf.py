#!/usr/bin/env python3

import os
import json
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.PDBIO import PDBIO
import re
import glob
import argparse
from tqdm import tqdm
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import plotly
from scipy.spatial.distance import euclidean


def get_chain_length(pdb_file, chain_id='A'):
    """Get the length of a specific chain from a PDB file."""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("temp", pdb_file)
        chain = structure[0][chain_id]
        return len(list(chain.get_residues()))
    except Exception as e:
        print(f"Error parsing {pdb_file}: {e}")
        return None

def calculate_rmsd(ref_pdb, target_pdb, ref_chain_id='C', target_chain_id='A'):
    """Calculate RMSD between two protein structures for specific chains."""
    parser = PDBParser(QUIET=True)
    
    try:
        ref_structure = parser.get_structure("reference", ref_pdb)
        target_structure = parser.get_structure("target", target_pdb)
        
        # Use separate chain IDs for each structure
        ref_chain = ref_structure[0][ref_chain_id]
        target_chain = target_structure[0][target_chain_id]
        
        ref_atoms = []
        target_atoms = []
        
        ref_residues = list(ref_chain.get_residues())
        target_residues = list(target_chain.get_residues())
        
        min_length = min(len(ref_residues), len(target_residues))
        
        for i in range(min_length):
            try:
                ref_atoms.append(ref_residues[i]["CA"])
                target_atoms.append(target_residues[i]["CA"])
            except KeyError:
                continue
        
        super_imposer = Superimposer()
        super_imposer.set_atoms(ref_atoms, target_atoms)
        
        return super_imposer.rms
    except Exception as e:
        print(f"Error calculating RMSD between {ref_pdb} and {target_pdb}: {e}")
        return None

#def extract_design_info(design_path):
#    """Extract design ID and target from path."""
#    # Pattern to match design_X.Y__vs__TARGET
#    match = re.search(r'design_(\d+)\.(\d+)__vs__(.+?)(?:/|$)', design_path)
#    if match:
#        backbone_id = match.group(1)
#        sequence_id = match.group(2)
#        target = match.group(3)
#        return backbone_id, sequence_id, target
#    return None, None, None
    
def extract_design_info(design_path):
    """Extract design ID and target from path."""
    # Pattern now captures any prefix before the numeric backbone_id.sequence_id
    match = re.search(r'(.+?)(\d+)\.(\d+)__vs__(.+?)(?:/|$)', design_path)
    if match:
        prefix = match.group(1)        # e.g. "HRGx_" or "design_"
        backbone_id = match.group(2)   # e.g. "42"
        sequence_id = match.group(3)   # e.g. "3"
        target = match.group(4)
        return prefix, backbone_id, sequence_id, target
    return None, None, None

def add_pae_rank(df):
    """
    Add pae_rank column to the dataframe based on min_interaction_pae values.
    Each design gets ranks 1-5 based on min_interaction_pae (lower is better).
    """
    # Create a new column for pae_rank
    # Group by design_id and assign ranks within each group based on min_interaction_pae
    df['pae_rank'] = df.groupby('design_id')['min_interaction_pae'].rank(method='first')
    return df


def score_designs(df):
    """
    Score miniprotein binder designs based on min_interaction_pae and rmsd_to_rfdiff.
    
    The function calculates the Euclidean distance from the ideal point (0,0) for each
    model rank, applies weights based on rank importance, and computes both raw and
    normalized scores.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the design data with columns:
        - design_id: unique identifier for each design
        - pae_rank: rank of the model (1-5)
        - min_interaction_pae: minimum interaction PAE value
        - rmsd_to_rfdiff: RMSD against the initial RFdiffusion backbone
    
    Returns:
    --------
    pandas.DataFrame
        The input DataFrame with additional columns:
        - raw_score: unscaled score (higher is better)
        - norm_score: normalized score between 0-1 (higher is better)
        - combined_score: a combined metric considering additional factors
    """
    # Define weights for each rank
    weights = {1: 0.1, 2: 0.08, 3: 0.06, 4: 0.04, 5: 0.02}
    
    # Create a copy of the input DataFrame to avoid modifying the original
    result_df = df.copy()
    
    # Calculate scores by design
    design_scores = {}
    
    # Get unique design IDs
    design_ids = df['design_id'].unique()
    
    for design_id in design_ids:
        # Get data for this design
        design_data = df[df['design_id'] == design_id]
        
        # Check if we have data for all ranks 1-5
        weighted_sum = 0
        total_weight = 0
        
        # Process each rank
        for rank in range(1, 6):
            rank_data = design_data[design_data['pae_rank'] == rank]
            
            if not rank_data.empty:
                # Get the first (and should be only) row for this rank
                row = rank_data.iloc[0]
                
                # Calculate Euclidean distance from ideal point (0,0)
                distance = euclidean(
                    [0, 0],
                    [row['min_interaction_pae'], row['rmsd_to_rfdiff']]
                )
                
                # Add weighted distance
                weighted_sum += distance * weights[rank]
                total_weight += weights[rank]
        
        # Calculate raw score (higher is better) - inverse of weighted distance
        if total_weight > 0 and weighted_sum > 0:
            raw_score = 1.0 / weighted_sum
        else:
            # Handle edge case where we might have division by zero
            raw_score = 0.0
            
        # Store the score for this design
        design_scores[design_id] = raw_score
    
    # Add raw scores to the DataFrame
    result_df['raw_score'] = result_df['design_id'].map(design_scores)
    
    # # Normalize scores to 0-1 range
    # min_score = result_df['raw_score'].min()
    # max_score = result_df['raw_score'].max()
    
    # if max_score > min_score:  # Avoid division by zero
    #     result_df['norm_score'] = (result_df['raw_score'] - min_score) / (max_score - min_score)
    # else:
    #     result_df['norm_score'] = 0.5  # Default value if all scores are the same
    
    # # Calculate a combined score that also considers pLDDT and intrachain PAE
    # # This is an alternative scoring approach that might better predict success
    # result_df['combined_score'] = (
    #     result_df['norm_score'] * 0.6 +  # Main score component
    #     result_df['mean_binder_plddt'] / 100.0 * 0.3 +  # Higher pLDDT is better
    #     (1.0 - result_df['mean_binder_intrachain_pae'] / result_df['mean_binder_intrachain_pae'].max()) * 0.1  # Lower intrachain PAE is better
    # )
    
    return result_df

def analyze_design(af_folder, rf_folder, output_csv):
    """Analyze all designs in the folder structure."""
    results = []
    design_folders = [f for f in os.listdir(af_folder) if os.path.isdir(os.path.join(af_folder, f))]
    
    for design_dir in tqdm(design_folders, desc="Analyzing designs"):
        # backbone_id, sequence_id, target = extract_design_info(design_dir)
        prefix, backbone_id, sequence_id, target = extract_design_info(design_dir)

        if backbone_id is None:
            continue
            
        design_path = os.path.join(af_folder, design_dir)
        model_files = glob.glob(os.path.join(design_path, f"*_unrelaxed_rank_*_alphafold2_multimer_*.pdb"))
        score_files = glob.glob(os.path.join(design_path, f"*_scores_rank_*_alphafold2_multimer_*.json"))
        
        # Match the RF diffusion backbone
        rf_backbone = os.path.join(rf_folder, f"{prefix}{backbone_id}.pdb")
        if not os.path.exists(rf_backbone):
            print(f"Warning: RF diffusion backbone not found: {rf_backbone}")
            continue
            
        # Loop through each model rank
        for i, (model_file, score_file) in enumerate(zip(sorted(model_files), sorted(score_files))):
            # Extract rank from filename
            rank_match = re.search(r'rank_(\d+)', os.path.basename(model_file))
            if not rank_match:
                continue
            rank = rank_match.group(1)
                
            # Get chain lengths
            binder_length = get_chain_length(model_file, chain_id='A')
            target_length = get_chain_length(model_file, chain_id='C')
            
            if binder_length is None or target_length is None:
                continue
            
            # Load scores from JSON
            try:
                with open(score_file, 'r') as f:
                    scores = json.load(f)
                    
                # Calculate metrics
                plddt = np.array(scores.get('plddt', []))
                pae = np.array(scores.get('pae', []))
                
                if len(plddt) == 0 or pae.shape[0] == 0:
                    print(f"Warning: Empty pLDDT or PAE data in {score_file}")
                    continue
                    
                # Get binder pLDDT (binder is chain A)
                binder_plddt = plddt[:binder_length]
                mean_binder_plddt = np.mean(binder_plddt)
                
                # Calculate intrachain PAE for binder (chain A)
                binder_pae = pae[:binder_length, :binder_length]
                mean_binder_intrachain_pae = np.mean(binder_pae)
                
                # Calculate interchain PAE (minimum interaction PAE)
                interchain_pae = pae[:binder_length, binder_length:binder_length+target_length]
                min_interaction_pae = np.min(interchain_pae)
                
                # Calculate RMSD between AF2 binder and RFdiffusion backbone
                rmsd = calculate_rmsd(rf_backbone, model_file, ref_chain_id='C', target_chain_id='A')
                
                # Collect results
                result = {
                    #'design_id': f"{backbone_id}.{sequence_id}",
                    'design_id': f"{prefix}{backbone_id}.{sequence_id}",
                    'target': target,
                    'rank': int(rank),
                    'mean_binder_plddt': mean_binder_plddt,
                    'mean_binder_intrachain_pae': mean_binder_intrachain_pae,
                    'min_interaction_pae': min_interaction_pae,
                    'rmsd_to_rfdiff': rmsd,
                    'binder_length': binder_length,
                    'target_length': target_length,
                    'model_file': os.path.basename(model_file),
                    'ptm': scores.get('ptm', None),
                    'iptm': scores.get('iptm', None),
                }
                results.append(result)
                
            except Exception as e:
                print(f"Error processing {score_file}: {e}")
    
    # Create DataFrame and save to CSV
    if results:
        df = pd.DataFrame(results)
        df.sort_values(['design_id', 'rank'], inplace=True)
        df = add_pae_rank(df)
        df = score_designs(df)
        df.to_csv(output_csv, index=False)
        print(f"Results saved to {output_csv}")
        return df
    else:
        print("No results found.")
        return None

def plot_metric_histograms(df, metric, cutoff, output_file, title_prefix="", x_label="", color_scheme=None):
    if color_scheme is None:
        color_scheme = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Create subplots with shared x-axis (removed subplot titles)
    fig = make_subplots(rows=5, cols=1, shared_xaxes=True, 
                        subplot_titles=[""] * 5,  # Empty strings to remove subtitles
                        vertical_spacing=0.05)
    
    # Calculate bin sizes based on the data range
    min_val = df[metric].min()
    max_val = df[metric].max()
    nbins = 30
    
    # Calculate bin edges once
    bin_edges = np.linspace(min_val, max_val, nbins + 1)
    
    # Storage for max frequencies and histogram data
    max_freqs = []
    all_hist_data = []
    
    # First pass to calculate histograms and find max frequency
    for rank in range(1, 6):
        rank_data = df[df['pae_rank'] == rank][metric]
        
        if len(rank_data) == 0:
            all_hist_data.append(None)
            continue
        
        # Calculate histogram without normalization
        hist, _ = np.histogram(rank_data, bins=bin_edges)
        
        # Normalize by the count to get probability
        hist_norm = hist / len(rank_data)
        
        # Store the normalized histogram
        all_hist_data.append(hist_norm)
        
        # Track max frequency
        if len(hist_norm) > 0:
            max_freqs.append(np.max(hist_norm))
    
    # Get overall max frequency with a small margin
    if max_freqs:
        max_freq = max(max_freqs) * 1.1
    else:
        max_freq = 1.0
    
    # Second pass to create plots with consistent y-axis
    for i, rank in enumerate(range(1, 6)):
        rank_data = df[df['pae_rank'] == rank][metric]
        
        if len(rank_data) == 0:
            continue
        
        # Create histogram using pre-calculated edges
        fig.add_trace(
            go.Histogram(
                x=rank_data,
                marker_color=color_scheme[i],
                name=f"PAE Rank {rank}",
                nbinsx=nbins,
                histnorm='probability',
                autobinx=False,
                xbins=dict(
                    start=min_val,
                    end=max_val,
                    size=(max_val - min_val) / nbins
                )
            ),
            row=i+1, col=1
        )
        
        # Add vertical line at cutoff
        fig.add_shape(
            type="line",
            x0=cutoff, y0=0,
            x1=cutoff, y1=max_freq,
            line=dict(color="red", width=2, dash="dash"),
            row=i+1, col=1
        )
        
        # Add annotation for count and percentage
        count_below = len(rank_data[rank_data <= cutoff])
        count_above = len(rank_data[rank_data > cutoff])
        pct_below = 100 * count_below / len(rank_data)
        
        # Modified annotations with increased offset
        fig.add_annotation(
            x=cutoff, y=0.9 * max_freq,
            text=f"≤ {cutoff}: {count_below} ({pct_below:.1f}%)",
            showarrow=True,
            arrowhead=1,
            ax=-60,  # Increased from -40 to -60
            ay=-20,
            row=i+1, col=1
        )
        
        fig.add_annotation(
            x=cutoff, y=0.8 * max_freq,
            text=f"> {cutoff}: {count_above} ({100-pct_below:.1f}%)",
            showarrow=True,
            arrowhead=1,
            ax=60,  # Increased from 40 to 60
            ay=-20,
            row=i+1, col=1
        )
        
        # Set consistent x and y ranges
        fig.update_xaxes(
            range=[min_val, max_val], 
            row=i+1, col=1,
            constrain='domain'
        )
        
        fig.update_yaxes(
            range=[0, max_freq], 
            row=i+1, col=1,
            constrain='domain'
        )
    
    # Update layout with centered title
    fig.update_layout(
        title=f"{title_prefix} by PAE Rank",
        title_x=0.5,  # Center the main title
        autosize=True,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=False
    )
    
    # Set x-axis title only for the bottom plot
    fig.update_xaxes(title_text=x_label, row=5, col=1)
    
    # Set y-axis titles with PAE rank information
    for i in range(1, 6):
        fig.update_yaxes(title_text=f"PAE rank {i}<br>Frequency", row=i, col=1)  # Added rank to y-axis
    
    # Save figure
    with open(output_file, 'w') as f:
        html = fig.to_html(config={
            'responsive': True,
            'scrollZoom': True,
            'displayModeBar': True,
            'modeBarButtonsToRemove': ['lasso2d', 'select2d'],
            'doubleClick': 'reset'
        }, include_plotlyjs=True)
        
        # Add responsive JS
        responsive_js = """
        <script>
        window.addEventListener('resize', function() {
            var gd = document.querySelector('.plotly-graph-div');
            if (gd) Plotly.Plots.resize(gd);
        });
        </script>
        """
        html = html.replace('</body>', f'{responsive_js}</body>')
        f.write(html)
    
    return fig

def plot_3d_connections(df, output_file):
    """
    Creates a 3D scatter plot with min_interaction_pae as x-axis, 
    rmsd_to_rfdiff as y-axis, and pae_rank as z-axis.
    Points from the same design are connected with lines.
    """
    # Create a figure
    fig = go.Figure()
    
    # Color options setup
    color_options = {
        'mean_binder_plddt': {'scale': px.colors.sequential.Viridis, 'min': df['mean_binder_plddt'].min(), 'max': df['mean_binder_plddt'].max()},
        'mean_binder_intrachain_pae': {'scale': px.colors.sequential.Plasma, 'min': df['mean_binder_intrachain_pae'].min(), 'max': df['mean_binder_intrachain_pae'].max()},
        'binder_length': {'scale': px.colors.sequential.Cividis, 'min': df['binder_length'].min(), 'max': df['binder_length'].max()},
        'ptm': {'scale': px.colors.sequential.Turbo, 'min': df['ptm'].min(), 'max': df['ptm'].max()},
        'iptm': {'scale': px.colors.sequential.Magma, 'min': df['iptm'].min(), 'max': df['iptm'].max()}
    }

    # Prepare data containers
    x, y, z = [], [], []
    hovertext = []
    color_data = {metric: [] for metric in color_options}
    design_indices = {}  # Dictionary to store indices by design_id for hover highlighting
    all_ids = []  # Store design_ids in order of appearance

    # Preprocess data
    for design_id in tqdm(df['design_id'].unique(), desc="Processing designs"):
        design_df = df[df['design_id'] == design_id].sort_values('pae_rank')
        points = len(design_df)
        
        # Store start index for this design
        start_idx = len(x)
        
        # Add coordinates
        x.extend(design_df['min_interaction_pae'].tolist())
        y.extend(design_df['rmsd_to_rfdiff'].tolist())
        z.extend(design_df['pae_rank'].tolist())
        
        # Store indices for hover highlighting
        design_indices[design_id] = list(range(start_idx, start_idx + points))
        all_ids.extend([design_id] * points)
        
        # Add NaN to break lines
        x.append(np.nan)
        y.append(np.nan)
        z.append(np.nan)
        all_ids.append(None)  # No design ID for the NaN separator

        # Create hover text
        for _, row in design_df.iterrows():
            hover_str = (
                f"Design: {row['design_id']}<br>"
                f"Target: {row['target']}<br>"
                f"PAE Rank: {row['pae_rank']}<br>"
                f"Rank: {row['rank']}<br>"
                f"Min Interaction PAE: {row['min_interaction_pae']:.2f}<br>"
                f"RMSD to RFdiff: {row['rmsd_to_rfdiff']:.2f}<br>"
                f"Mean Binder pLDDT: {row['mean_binder_plddt']:.2f}<br>"
                f"Mean Binder IntraChain PAE: {row['mean_binder_intrachain_pae']:.2f}<br>"
                f"Binder Length: {row['binder_length']}<br>"
                f"Target Length: {row['target_length']}<br>"
                f"pTM: {row['ptm']:.2f}<br>"
                f"ipTM: {row['iptm']:.2f}<br>"
                f"Score: {row['raw_score']:.2f}"
            )
            hovertext.append(hover_str)
            for metric in color_options:
                color_data[metric].append(row[metric])
        
        # Add NaN values for color arrays
        hovertext.append('')
        for metric in color_options:
            color_data[metric].append(np.nan)

    # Add main scatter trace
    initial_metric = 'mean_binder_plddt'
    scatter_trace = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='lines+markers',
        line=dict(width=2, color='rgba(150,150,150,0.5)'),
        marker=dict(
            size=5,  # Fixed marker size
            color=color_data[initial_metric],
            colorscale=color_options[initial_metric]['scale'],
            cmin=color_options[initial_metric]['min'],
            cmax=color_options[initial_metric]['max'],
            opacity=0.8,
            colorbar=dict(title=initial_metric)
        ),
        hovertext=hovertext,
        hoverinfo='text',
        customdata=all_ids,  # Store design IDs for highlighting
        name='Design Points',
        showlegend=False
    )
    fig.add_trace(scatter_trace)

    # Calculate max values for axes, ensuring they're floats
    max_vals = {
        'x': float(df['min_interaction_pae'].max()),
        'y': float(df['rmsd_to_rfdiff'].max()),
        'z': float(df['pae_rank'].max())
    }
    
    # Add explicit axis arrows using scatter3d instead of cone
    # X-axis arrow
    fig.add_trace(go.Scatter3d(
        x=[0, max_vals['x']],
        y=[0, 0],
        z=[1, 1],
        mode='lines',
        line=dict(color='black', width=6),
        hoverinfo='none',
        showlegend=False
    ))
    
    # Y-axis arrow
    fig.add_trace(go.Scatter3d(
        x=[0, 0],
        y=[0, max_vals['y']],
        z=[1, 1],
        mode='lines',
        line=dict(color='black', width=6),
        hoverinfo='none',
        showlegend=False
    ))
    
    # Z-axis arrow
    fig.add_trace(go.Scatter3d(
        x=[0, 0],
        y=[0, 0],
        z=[1, max_vals['z']],
        mode='lines',
        line=dict(color='black', width=6),
        hoverinfo='none',
        showlegend=False
    ))

    # Create dropdown menu options
    dropdown_options = []
    for metric in color_options:
        dropdown_option = dict(
            method='update',
            label=metric,
            args=[
                {'marker.color': [color_data[metric]],
                 'marker.colorscale': [color_options[metric]['scale']],
                 'marker.cmin': [color_options[metric]['min']],
                 'marker.cmax': [color_options[metric]['max']],
                 'marker.colorbar.title': metric},
                {'title': f'3D Binder Design Space - Color by {metric}'}
            ]
        )
        dropdown_options.append(dropdown_option)

    # Update layout
    fig.update_layout(
        title='3D Binder Design Space - Color by mean_binder_plddt',
        title_x=0.5,
        scene=dict(
            xaxis_title='Min Interaction PAE',
            yaxis_title='RMSD to RFdiff',
            zaxis_title='PAE Rank',
            xaxis=dict(range=[-0.5, max_vals['x']], showspikes=False),
            yaxis=dict(range=[-0.5, max_vals['y']], showspikes=False),
            zaxis=dict(range=[0.9, 5.5], dtick=1, showspikes=False)
        ),
        updatemenus=[dict(
            type='dropdown',
            showactive=True,
            buttons=dropdown_options,
            x=0.0,  # Positioned at left
            xanchor='left',
            y=1.04,  # Positioned above the plot
            yanchor='top',
            pad={"r": 10, "t": 10},
            bgcolor="white",
            bordercolor="lightgrey",
            font=dict(size=12)
        )],
        autosize=True,
        margin=dict(l=0, r=0, b=0, t=40),
        scene_camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
    )
    
    # Save figure with optimized JavaScript for hover highlighting
    with open(output_file, 'w') as f:
        plot_html = plotly.offline.plot(fig, include_plotlyjs='cdn', output_type='div', config={'responsive': True})
        
        # Prepare a more optimized version of design indices for JS
        js_design_indices = {}
        # Only include necessary information - map design ID to list of point indices
        for design_id, indices in design_indices.items():
            # Use string keys for JS
            js_design_indices[str(design_id)] = indices
        
        # Insert custom JavaScript for hover highlighting - optimized version
        highlight_script = f"""
        <script>
        document.addEventListener('DOMContentLoaded', function() {{
            var graphDiv = document.querySelector('.js-plotly-plot');
            if (!graphDiv) return;
            
            // Design indices mapping - optimized to only include necessary data
            const designIndices = {json.dumps(js_design_indices)};
            
            // Store original styles
            const originalStyles = {{
                lineColor: 'rgba(150,150,150,0.5)',
                lineWidth: 2,
                markerSize: 5,
                markerOpacity: 0.8,
                markerLineColor: 'rgba(0,0,0,0)',
                markerLineWidth: 0
            }};
            
            // Track state to prevent unnecessary updates
            let highlightedDesign = null;
            let isUpdating = false;
            
            // Throttle function to limit the rate at which a function can fire
            function throttle(func, limit) {{
                let inThrottle;
                return function() {{
                    const args = arguments;
                    const context = this;
                    if (!inThrottle) {{
                        func.apply(context, args);
                        inThrottle = true;
                        setTimeout(() => inThrottle = false, limit);
                    }}
                }};
            }}
            
            // Optimized highlight function
            function highlightDesign(design_id) {{
                if (isUpdating || highlightedDesign === design_id) return;
                isUpdating = true;
                
                try {{
                    // Get the indices for this design
                    const indices = designIndices[design_id];
                    if (!indices || indices.length === 0) {{
                        isUpdating = false;
                        return;
                    }}
                    
                    // Create arrays for the updates, initialized with defaults
                    const pointCount = graphDiv.data[0].x.length;
                    const sizes = Array(pointCount).fill(5);
                    const opacities = Array(pointCount).fill(0.8);
                    const borderColors = Array(pointCount).fill('rgba(0,0,0,0)');
                    const borderWidths = Array(pointCount).fill(0);
                    const lineColors = Array(pointCount).fill(originalStyles.lineColor);
                    const lineWidths = Array(pointCount).fill(originalStyles.lineWidth);
                    
                    // Only update the specific points and connecting lines we care about
                    indices.forEach((idx, i) => {{
                        if (idx < pointCount) {{
                            // Highlight markers
                            sizes[idx] = 14;
                            opacities[idx] = 1.0;
                            borderColors[idx] = 'black';
                            borderWidths[idx] = 3;  // Thicker borders as requested
                            
                            // Highlight connecting lines
                            if (i < indices.length - 1) {{
                                const nextIdx = indices[i + 1];
                                // Set line color for this segment
                                lineColors[idx] = 'black';
                                lineWidths[idx] = 4;  // Thicker line for visibility
                            }}
                        }}
                    }});
                    
                    // Apply all updates at once
                    const update = {{
                        'marker.size': [sizes],
                        'marker.opacity': [opacities],
                        'marker.line.color': [borderColors],
                        'marker.line.width': [borderWidths],
                        'line.color': [lineColors],
                        'line.width': [lineWidths]
                    }};
                    
                    // Apply the updates
                    Plotly.restyle(graphDiv, update, [0]).then(() => {{
                        highlightedDesign = design_id;
                        isUpdating = false;
                    }}).catch(() => {{
                        isUpdating = false;
                    }});
                }} catch (e) {{
                    console.error("Error in highlighting:", e);
                    isUpdating = false;
                }}
            }}
            
            // Reset function
            function resetHighlight() {{
                if (!highlightedDesign || isUpdating) return;
                isUpdating = true;
            
                const pointCount = graphDiv.data[0].x.length;
                
                // Reset everything to original values
                Plotly.restyle(graphDiv, {{
                    'marker.size': originalStyles.markerSize,
                    'marker.opacity': originalStyles.markerOpacity,
                    'marker.line.color': Array(pointCount).fill('rgba(0,0,0,0)'),
                    'marker.line.width': Array(pointCount).fill(0),
                    'line.color': Array(pointCount).fill(originalStyles.lineColor),
                    'line.width': Array(pointCount).fill(originalStyles.lineWidth)
                }}, [0]).then(() => {{
                    highlightedDesign = null;
                    isUpdating = false;
                }}).catch(() => {{
                    isUpdating = false;
                }});
            }}
        
            // Use throttled versions of our functions
            const throttledHighlight = throttle(highlightDesign, 100);
            const throttledReset = throttle(resetHighlight, 100);
            
            // Attach optimized event handlers
            graphDiv.on('plotly_hover', function(data) {{
                if (!data.points || data.points.length === 0) return;
                const point = data.points[0];
                if (!point.customdata) return;
                
                throttledHighlight(point.customdata);
            }});
            
            graphDiv.on('plotly_unhover', function() {{
                throttledReset();
            }});
        }});
        </script>
        """
        
        # Insert the script right before the closing body tag
        html_with_script = plot_html + highlight_script
        f.write(html_with_script)
    
    return fig

def create_main_report(output_files, output_html):
    """
    Creates a main HTML report with a side menu to navigate between plots.
    """
    # Generate the HTML content
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Binder Design Analysis Report</title>
        <style>
            html, body {
                margin: 0;
                padding: 0;
                width: 100%;
                height: 100%;
                overflow: hidden;
                font-family: Arial, sans-serif;
            }
            .container {
                display: flex;
                width: 100%;
                height: 100vh;
                overflow: hidden;
            }
            #sidebar {
                width: 250px;
                background-color: #f5f5f5;
                height: 100%;
                overflow-y: auto;
                padding: 20px;
                box-sizing: border-box;
                flex-shrink: 0;
                box-shadow: 2px 0 5px rgba(0,0,0,0.1);
                z-index: 100;
            }
            #content {
                flex-grow: 1;
                height: 100%;
                overflow: hidden;
                position: relative;
            }
            iframe {
                width: 100%;
                height: 100%;
                border: none;
                position: absolute;
                top: 0;
                left: 0;
                right: 0;
                bottom: 0;
            }
            #sidebar h2 {
                margin-top: 0;
                border-bottom: 1px solid #ddd;
                padding-bottom: 10px;
            }
            #sidebar ul {
                list-style-type: none;
                padding: 0;
                margin-top: 20px;
            }
            #sidebar li {
                margin-bottom: 10px;
            }
            #sidebar a {
                text-decoration: none;
                color: #0066cc;
                padding: 8px 5px;
                display: block;
                border-radius: 5px;
                transition: background-color 0.2s;
            }
            #sidebar a:hover {
                background-color: #e0e0e0;
            }
            .active {
                background-color: #e0e0e0;
                font-weight: bold;
            }
            .toggle-sidebar {
                position: absolute;
                top: 10px;
                left: 10px;
                z-index: 200;
                background: #f5f5f5;
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 5px 10px;
                display: none;
                cursor: pointer;
            }
            @media (max-width: 768px) {
                #sidebar {
                    position: fixed;
                    left: -250px;
                    transition: left 0.3s ease;
                }
                #sidebar.show {
                    left: 0;
                }
                #content {
                    margin-left: 0;
                    width: 100%;
                }
                .toggle-sidebar {
                    display: block;
                }
            }
        </style>
    </head>
    <body>
        <div class="container">
            <div id="sidebar">
                <h2>Analysis Plots</h2>
                <ul>
    """
    
    # Add menu items
    first_item = True
    for name, file_path in output_files.items():
        file_name = os.path.basename(file_path)
        active_class = "active" if first_item else ""
        html_content += f'<li><a href="javascript:void(0)" class="{active_class}" onclick="loadIframe(\'{file_name}\', \'{name}\')">{name}</a></li>\n'
        first_item = False
    
    # Add iframe and JavaScript
    first_file = os.path.basename(list(output_files.values())[0]) if output_files else ""
    first_name = list(output_files.keys())[0] if output_files else ""
    
    html_content += f"""
                </ul>
            </div>
            <button class="toggle-sidebar" onclick="toggleSidebar()">☰</button>
            <div id="content">
                <iframe id="plot-frame" src="{first_file}" allowfullscreen></iframe>
            </div>
        </div>
        <script>
            function loadIframe(file, name) {{
                const iframe = document.getElementById('plot-frame');
                iframe.src = file;
                document.title = 'Binder Analysis: ' + name;
                
                // Update active class
                const links = document.querySelectorAll('#sidebar a');
                links.forEach(link => link.classList.remove('active'));
                
                // Find the clicked link
                const clickedLink = Array.from(links).find(link => 
                    link.textContent === name);
                if (clickedLink) clickedLink.classList.add('active');
                
                // Close sidebar on mobile after selection
                if (window.innerWidth <= 768) {{
                    document.getElementById('sidebar').classList.remove('show');
                }}
            }}
            
            function toggleSidebar() {{
                document.getElementById('sidebar').classList.toggle('show');
            }}
            
            // Force iframe content to resize when iframe loads
            document.getElementById('plot-frame').addEventListener('load', function() {{
                try {{
                    const innerWindow = this.contentWindow;
                    if (innerWindow && innerWindow.Plotly) {{
                        const gd = innerWindow.document.querySelector('.plotly-graph-div');
                        if (gd) {{
                            innerWindow.Plotly.Plots.resize(gd);
                        }}
                    }}
                }} catch (e) {{
                    console.error('Error resizing plot:', e);
                }}
            }});
            
            // Handle window resize
            window.addEventListener('resize', function() {{
                try {{
                    const iframe = document.getElementById('plot-frame');
                    const innerWindow = iframe.contentWindow;
                    if (innerWindow && innerWindow.Plotly) {{
                        const gd = innerWindow.document.querySelector('.plotly-graph-div');
                        if (gd) {{
                            innerWindow.Plotly.Plots.resize(gd);
                        }}
                    }}
                }} catch (e) {{
                    console.error('Error resizing plot:', e);
                }}
            }});
            
            // Initialize
            document.addEventListener('DOMContentLoaded', function() {{
                document.title = 'Binder Analysis: {first_name}';
            }});
        </script>
    </body>
    </html>
    """
    
    # Write the HTML file
    with open(output_html, 'w') as f:
        f.write(html_content)
    
    return output_html

def generate_all_plots(df, output_dir="plots"):
    """
    Generate all plots and save them to the specified directory.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        The dataframe containing the data
    output_dir : str, optional
        Directory to save the plots to
    
    Returns:
    --------
    dict
        Dictionary with plot names as keys and file paths as values
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # # Add pae_rank to the dataframe
    # df = add_pae_rank(df)
    
    # Dictionary to store output files
    output_files = {}
    
    # Generate miPAE histograms
    output_file = os.path.join(output_dir, "mipae_histograms.html")
    plot_metric_histograms(
        df, "min_interaction_pae", 10, output_file,
        title_prefix="Minimum Interaction PAE", x_label="Min Interaction PAE"
    )
    output_files["Minimum Interaction PAE Histograms"] = output_file
    
    # Generate RMSD histograms
    output_file = os.path.join(output_dir, "rmsd_histograms.html")
    plot_metric_histograms(
        df, "rmsd_to_rfdiff", 3, output_file,
        title_prefix="RMSD to RFdiffusion", x_label="RMSD (Å)"
    )
    output_files["RMSD to RFdiffusion Histograms"] = output_file
    
    # Generate pLDDT histograms
    output_file = os.path.join(output_dir, "plddt_histograms.html")
    plot_metric_histograms(
        df, "mean_binder_plddt", 60, output_file,
        title_prefix="Mean Binder pLDDT", x_label="Mean pLDDT"
    )
    output_files["Mean Binder pLDDT Histograms"] = output_file
    
    # Generate ipTM histograms
    output_file = os.path.join(output_dir, "iptm_histograms.html")
    plot_metric_histograms(
        df, "iptm", 0.5, output_file,
        title_prefix="Interface pTM Score", x_label="ipTM"
    )
    output_files["Interface pTM Histograms"] = output_file
    
    # Generate intrachain PAE histograms
    output_file = os.path.join(output_dir, "intrachain_pae_histograms.html")
    plot_metric_histograms(
        df, "mean_binder_intrachain_pae", 10, output_file,
        title_prefix="Mean Binder Intrachain PAE", x_label="Mean Intrachain PAE"
    )
    output_files["Mean Binder Intrachain PAE Histograms"] = output_file
    
    # Generate 3D scatter plot
    output_file = os.path.join(output_dir, "3d_scatter.html")
    plot_3d_connections(df, output_file)
    output_files["3D Design Space Visualization"] = output_file
    
    # Create main report
    main_html = os.path.join(output_dir, "binder_analysis_report.html")
    create_main_report(output_files, main_html)
    
    print(f"All plots generated. Main report saved to {main_html}")
    return main_html


def main():
    parser = argparse.ArgumentParser(description='Analyze protein binder designs')
    parser.add_argument('--af_folder', required=True, help='Path to AlphaFold2 results folder')
    parser.add_argument('--rf_folder', required=True, help='Path to RFdiffusion designs folder')
    parser.add_argument('--output_dir', default = "analysis_plots", help='Output CSV file')
    
    args = parser.parse_args()
    
    # Progress
    print(f"Analyzing designs in {args.af_folder}")
    print(f"Using RF diffusion designs from {args.rf_folder}")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate and save the metrics
    default_csv=args.output_dir + '/design_analysis_results.csv'
    df = analyze_design(args.af_folder, args.rf_folder, default_csv)

    # Generate and save plots
    generate_all_plots(df, output_dir= args.output_dir)
        
    # Save ranking file
    default_ranking_file = args.output_dir + '/designs_ranking.csv'
    scoring_columns = ['design_id', 'raw_score']				# Select columns to keep
    filtered_df = df[scoring_columns]
    unique_df = filtered_df.drop_duplicates()					# Drop duplicate rows to keep only unique entries
    sorted_by_raw = unique_df.sort_values(by='raw_score', ascending=False)	# Sort by 'raw_score'
    sorted_by_raw.to_csv(default_ranking_file, index=False)

    
    if df is not None:
        # Generate summary statistics
        print("\nSummary Statistics:")
        print(f"Total designs analyzed: {df['design_id'].nunique()}")
        print("\nMean pLDDT statistics:")
        print(df.groupby('design_id')['mean_binder_plddt'].agg(['mean', 'min', 'max']).describe())
        
        print("\nMinimum Interaction PAE statistics:")
        print(df.groupby('design_id')['min_interaction_pae'].agg(['mean', 'min', 'max']).describe())
        
        print("\nRMSD to RF diffusion statistics:")
        print(df.groupby('design_id')['rmsd_to_rfdiff'].agg(['mean', 'min', 'max']).describe())
        
        # Find best designs by different metrics
        print("\nTop 5 designs by mean pLDDT (rank 1 models only):")
        top_plddt = df[df['rank'] == 1].sort_values('mean_binder_plddt', ascending=False).head(5)
        print(top_plddt[['design_id', 'mean_binder_plddt', 'min_interaction_pae', 'rmsd_to_rfdiff']])
        
        print("\nTop 5 designs by min interaction PAE (rank 1 models only):")
        top_pae = df[df['rank'] == 1].sort_values('min_interaction_pae').head(5)
        print(top_pae[['design_id', 'mean_binder_plddt', 'min_interaction_pae', 'rmsd_to_rfdiff']])
        
        print("\nTop 5 designs by low RMSD to RF diffusion (rank 1 models only):")
        top_rmsd = df[df['rank'] == 1].sort_values('rmsd_to_rfdiff').head(5)
        print(top_rmsd[['design_id', 'mean_binder_plddt', 'min_interaction_pae', 'rmsd_to_rfdiff']])

if __name__ == "__main__":
    main()
