
# import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean



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

# Example usage
# scored_df = score_designs(df)

df = pd.read_csv("./design_analysis_results.csv", sep = ",")

# ['design_id', 'target', 'rank', 'mean_binder_plddt',
# 'mean_binder_intrachain_pae', 'min_interaction_pae', 'rmsd_to_rfdiff',
# 'binder_length', 'target_length', 'model_file', 'ptm', 'iptm']
df.columns



# First, select the columns you want to keep
scoring_columns = ['design_id', 'raw_score']
filtered_df = df[scoring_columns]

# Then, drop duplicate rows to keep only unique entries
unique_df = filtered_df.drop_duplicates()

# Sort by 'raw_score'
sorted_by_raw = unique_df.sort_values(by='raw_score', ascending=False)
sorted_by_raw.to_csv("designs_ranking.csv", index = False)

