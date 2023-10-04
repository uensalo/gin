import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde


def plot_side_by_side_scatterplots(filepath1, filepath2, threshold, output_path):
    """
    Plot two scatterplots side by side based on vertex length (x-axis) and vertex outdegree (y-axis) for given graphs.
    
    Args:
    - filepath1 (str): Path to the first graph file.
    - filepath2 (str): Path to the second graph file.
    - threshold (int): Threshold for label length.
    - output_path (str): Path to save the scatterplots.
    """
    # Function to prepare data for plotting
    def prepare_data(filepath):
        data = pd.read_csv(filepath, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
        vertex_data = data[data['Type'] == 'V']
        vertex_data = vertex_data[vertex_data['Target'].str.len() <= threshold]
        outdegree = data[data['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
        merged_data = vertex_data.merge(outdegree, left_on='Source', right_index=True, how='left').fillna(0)
        merged_data['Length'] = merged_data['Target'].str.len()
        return merged_data
    
    # Prepare data for both graphs
    merged_data1 = prepare_data(filepath1)
    merged_data2 = prepare_data(filepath2)
    
    # Compute common x and y limits for both scatterplots
    x_lim = (min(merged_data1['Length'].min(), merged_data2['Length'].min()),
             max(merged_data1['Length'].max(), merged_data2['Length'].max()))
    y_lim = (min(merged_data1['Outdegree'].min(), merged_data2['Outdegree'].min()),
             max(merged_data1['Outdegree'].max(), merged_data2['Outdegree'].max()))
    
    # Plot side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))
    
    ax1.scatter(merged_data1['Length'], merged_data1['Outdegree'], alpha=0.5)
    ax1.set_xscale('log')
    ax1.set_xlim(x_lim)
    ax1.set_ylim(y_lim)
    ax1.set_xlabel('Vertex Length')
    ax1.set_ylabel('Outdegree')
    ax1.set_title('Scatterplot for GRCh38-20-0.10b')
    ax1.grid(True, which="both", ls="--", c='0.7')
    
    ax2.scatter(merged_data2['Length'], merged_data2['Outdegree'], alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_xlim(x_lim)
    ax2.set_ylim(y_lim)
    ax2.set_xlabel('Vertex Length')
    ax2.set_ylabel('Outdegree')
    ax2.set_title('Scatterplot for gencode.v40')
    ax2.grid(True, which="both", ls="--", c='0.7')
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def compute_kde_values(dataframe, bandwidth = 5):
    """
    Compute KDE values for a given dataframe of vertex lengths and outdegrees.
    
    Args:
    - dataframe (pd.DataFrame): DataFrame with 'Length' and 'Outdegree' columns.
    
    Returns:
    - np.array: Grid of KDE values.
    - (min_length, max_length): Minimum and maximum vertex lengths.
    - (min_outdegree, max_outdegree): Minimum and maximum outdegrees.
    """
    # KDE computation
    kde = gaussian_kde(np.vstack([dataframe['Length'], dataframe['Outdegree']]), bw_method=bandwidth)
    
    # Define grid
    x_grid = np.linspace(dataframe['Length'].min(), dataframe['Length'].max(), 100)
    y_grid = np.linspace(0, dataframe['Outdegree'].max(), 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = kde.evaluate(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    
    return Z, (x_grid.min(), x_grid.max()), (y_grid.min(), y_grid.max())

def plot_side_by_side_kde_heatmaps(filepath1, filepath2, threshold, output_path, bandwidth=2):
    """
    Plot two KDE-smoothed heatmaps side by side based on vertex length (x-axis) and vertex outdegree (y-axis) for given graphs.
    
    Args:
    - filepath1 (str): Path to the first graph file.
    - filepath2 (str): Path to the second graph file.
    - threshold (int): Threshold for label length.
    - output_path (str): Path to save the heatmap.
    - bandwidth (float): Bandwidth for KDE. Smaller values produce more localized smoothing.
    """
    data1 = pd.read_csv(filepath1, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
    data2 = pd.read_csv(filepath2, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
    
    # Corrected filtering
    vertex_data1 = data1[(data1['Type'] == 'V') & (data1['Target'].str.len() <= threshold)]
    vertex_data2 = data2[(data2['Type'] == 'V') & (data2['Target'].str.len() <= threshold)]
    
    # Calculate outdegree for each vertex and provide a name
    outdegree1 = data1[data1['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
    outdegree2 = data2[data2['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
    
    # Merge vertex data and outdegree
    merged_data1 = vertex_data1.merge(outdegree1, left_on='Source', right_index=True, how='left').fillna(0)
    merged_data2 = vertex_data2.merge(outdegree2, left_on='Source', right_index=True, how='left').fillna(0)
    
    merged_data1['Length'] = merged_data1['Target'].str.len()
    merged_data2['Length'] = merged_data2['Target'].str.len()
    
    # Compute KDE values
    Z1, (min_length1, max_length1), (min_outdegree1, max_outdegree1) = compute_kde_values(merged_data1)
    Z2, (min_length2, max_length2), (min_outdegree2, max_outdegree2) = compute_kde_values(merged_data2)
    
    # Set common scale for color based on the maximum KDE value across both datasets
    vmax = max(Z1.max(), Z2.max())
    
    # Plot side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))
    
    sns.heatmap(Z1, cmap="YlGnBu", cbar_kws={'label': 'Density'}, ax=ax1, vmax=vmax)
    ax1.set_title('GRCh38-20-0.10b')
    ax1.set_xlabel('Vertex Length')
    ax1.set_ylabel('Outdegree')
    
    sns.heatmap(Z2, cmap="YlGnBu", cbar_kws={'label': 'Density'}, ax=ax2, vmax=vmax)
    ax2.set_title('gencode.v40')
    ax2.set_xlabel('Vertex Length')
    ax2.set_ylabel('Outdegree')
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

plot_side_by_side_scatterplots('../input/GRCh38-20-0.10b.fmdg', '../input/gencode.v40.fmdg', 5000000, '../plot/scatter_degree_len.png')
plot_side_by_side_kde_heatmaps('../input/GRCh38-20-0.10b.fmdg', '../input/gencode.v40.fmdg', 5000000, '../plot/heatmap_degree_len.png')
