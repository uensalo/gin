import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde



def parse_graph_data(filepath, threshold):
    """
    Parse graph data from a file and return filtered vertex lengths and outdegrees.
    
    Args:
    - filepath (str): Path to the graph file.
    - threshold (int): Threshold for label length.
    
    Returns:
    - pd.DataFrame: DataFrame with 'Length' and 'Outdegree' columns.
    """
    data = pd.read_csv(filepath, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
    
    vertex_data = data[(data['Type'] == 'V') & (data['Target'].str.len() <= threshold)]
    outdegree = data[data['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
    merged_data = vertex_data.merge(outdegree, left_on='Source', right_index=True, how='left').fillna(0)
    merged_data['Length'] = merged_data['Target'].str.len()
    
    return merged_data[['Length', 'Outdegree']]


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
    ax1.set_yscale('log')
    ax1.set_xlim(x_lim)
    ax1.set_ylim(y_lim)
    ax1.set_xlabel('Vertex Length')
    ax1.set_ylabel('Outdegree')
    ax1.set_title('Scatterplot for GRCh38-20-0.10b')
    ax1.grid(True, which="both", ls="--", c='0.7')
    
    ax2.scatter(merged_data2['Length'], merged_data2['Outdegree'], alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(x_lim)
    ax2.set_ylim(y_lim)
    ax2.set_xlabel('Vertex Length')
    ax2.set_ylabel('Outdegree')
    ax2.set_title('Scatterplot for gencode.v40')
    ax2.grid(True, which="both", ls="--", c='0.7')
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def compute_kde_values(dataframe, bandwidth=5):
    """
    Compute KDE values for a given dataframe of vertex lengths and outdegrees.
    
    Args:
    - dataframe (pd.DataFrame): DataFrame with 'Length' and 'Outdegree' columns.
    - bandwidth (float): Bandwidth for KDE.
    
    Returns:
    - np.array: Grid of KDE values.
    - (min_length, max_length): Minimum and maximum vertex lengths.
    - (min_outdegree, max_outdegree): Minimum and maximum outdegrees.
    """
    kde = gaussian_kde(np.vstack([dataframe['Length'], dataframe['Outdegree']]), bw_method=bandwidth)
    x_grid = np.linspace(dataframe['Length'].min(), dataframe['Length'].max(), 100)
    y_grid = np.linspace(0, dataframe['Outdegree'].max(), 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = kde.evaluate(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    return Z, (x_grid.min(), x_grid.max()), (y_grid.min(), y_grid.max())

# def plot_side_by_side_kde_heatmaps(filepath1, filepath2, threshold, output_path, bandwidth=2):
#     data1 = pd.read_csv(filepath1, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
#     data2 = pd.read_csv(filepath2, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
    
#     vertex_data1 = data1[(data1['Type'] == 'V') & (data1['Target'].str.len() <= threshold)]
#     vertex_data2 = data2[(data2['Type'] == 'V') & (data2['Target'].str.len() <= threshold)]
    
#     outdegree1 = data1[data1['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
#     outdegree2 = data2[data2['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
    
#     merged_data1 = vertex_data1.merge(outdegree1, left_on='Source', right_index=True, how='left').fillna(0)
#     merged_data2 = vertex_data2.merge(outdegree2, left_on='Source', right_index=True, how='left').fillna(0)
    
#     merged_data1['Length'] = merged_data1['Target'].str.len()
#     merged_data2['Length'] = merged_data2['Target'].str.len()
    
#     Z1, (min_length1, max_length1), (min_outdegree1, max_outdegree1) = compute_kde_values(merged_data1, bandwidth)
#     Z2, (min_length2, max_length2), (min_outdegree2, max_outdegree2) = compute_kde_values(merged_data2, bandwidth)
    
#     Z1 /= Z1.sum()
#     Z2 /= Z2.sum()
    
#     vmax = max(Z1.max(), Z2.max())

#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8))
    
#     sns.heatmap(np.log(Z1 + 1e-10), cmap="YlGnBu", cbar_kws={'label': 'Log Density'}, ax=ax1, vmax=np.log(vmax + 1e-10))
#     ax1.set_title('GRCh38-20-0.10b')
#     ax1.set_xlabel('Log Vertex Length')
#     ax1.set_ylabel('Log Outdegree')
#     ax1.set_xlim(np.log(min_length1), np.log(max_length1))
#     ax1.set_ylim(1, max_outdegree1)
    
#     sns.heatmap(np.log(Z2 + 1e-10), cmap="YlGnBu", cbar_kws={'label': 'Log Density'}, ax=ax2, vmax=np.log(vmax + 1e-10))
#     ax2.set_title('gencode.v40')
#     ax2.set_xlabel('Log Vertex Length')
#     ax2.set_ylabel('Log Outdegree')
#     ax2.set_xlim(np.log(min_length2), np.log(max_length2))
#     ax2.set_ylim(1, max_outdegree2)
    
#     plt.tight_layout()
#     plt.savefig(output_path)
#     plt.close()

def plot_side_by_side_heatmaps(filepath1, filepath2, threshold, output_path):
    """
    Plot two heatmaps side by side based on vertex length (x-axis) and vertex outdegree (y-axis) for given graphs.
    
    Args:
    - filepath1 (str): Path to the first graph file.
    - filepath2 (str): Path to the second graph file.
    - threshold (int): Threshold for label length.
    - output_path (str): Path to save the heatmap.
    """
    data1 = parse_graph_data(filepath1, threshold)
    data2 = parse_graph_data(filepath2, threshold)
    
    # Determine global min and max for vertex lengths and outdegrees
    global_min_length = min(data1['Length'].min(), data2['Length'].min())
    global_max_length = max(data1['Length'].max(), data2['Length'].max())
    global_max_outdegree = max(data1['Outdegree'].max(), data2['Outdegree'].max())
    
    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 8), sharex=True, sharey=True)
    
    # Using 2D histogram to generate heatmap data
    c1 = ax1.hist2d(np.log(data1['Length']), np.log1p(data1['Outdegree']), bins=[100, 100], 
                   range=[[np.log(global_min_length), np.log(global_max_length)], [0, np.log1p(global_max_outdegree)]], 
                   cmap="YlGnBu", density=True)
    ax1.set_title('GRCh38-20-0.10b')
    ax1.set_xlabel('Log Vertex Length')
    ax1.set_ylabel('Log Outdegree')
    
    c2 = ax2.hist2d(np.log(data2['Length']), np.log1p(data2['Outdegree']), bins=[100, 100], 
                   range=[[np.log(global_min_length), np.log(global_max_length)], [0, np.log1p(global_max_outdegree)]], 
                   cmap="YlGnBu", density=True)
    ax2.set_title('gencode.v40')
    ax2.set_xlabel('Log Vertex Length')
    
    # Add a colorbar to the plot
    cbar = fig.colorbar(c1[3], ax=[ax1, ax2], orientation='vertical', fraction=0.05)
    cbar.ax.set_ylabel('Density')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

#plot_side_by_side_scatterplots('../input/GRCh38-20-0.10b.fmdg', '../input/gencode.v40.fmdg', 5000000, '../plot/scatter_degree_len.png')
#plot_side_by_side_kde_heatmaps('../input/GRCh38-20-0.10b.fmdg', '../input/gencode.v40.fmdg', 5000000, '../plot/heatmap_degree_len.png')
plot_side_by_side_heatmaps('../input/GRCh38-20-0.10b.fmdg', '../input/gencode.v40.fmdg', 5000000, '../plot/heatmap_degree_len.png')
