import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde



def parse_graph_data(filepath, threshold):
    """
    Parse graph data from a file and return filtered vertex lengths and degrees (sum of indegrees and outdegrees).
    
    Args:
    - filepath (str): Path to the graph file.
    - threshold (int): Threshold for label length.
    
    Returns:
    - pd.DataFrame: DataFrame with 'Length' and 'Degree' columns.
    """
    data = pd.read_csv(filepath, sep="\t", header=None, names=["Type", "Source", "Target"], low_memory=False)
    
    vertex_data = data[(data['Type'] == 'V') & (data['Target'].str.len() <= threshold)]
    indegree = data[data['Type'] == 'E'].groupby('Target').size().rename('Indegree')
    outdegree = data[data['Type'] == 'E'].groupby('Source').size().rename('Outdegree')
    
    # Convert the 'Source' column to string type to match the index type of indegree and outdegree
    vertex_data.loc[:, 'Source'] = vertex_data['Source'].astype(str)
    # Convert the indices of indegree and outdegree to strings
    indegree.index = indegree.index.astype(str)
    outdegree.index = outdegree.index.astype(str)
    
    merged_data = vertex_data.merge(outdegree, left_on='Source', right_index=True, how='left').fillna(0)
    merged_data = merged_data.merge(indegree, left_on='Source', right_index=True, how='left').fillna(0)
    
    # Compute total degree as sum of indegree and outdegree
    merged_data['Degree'] = merged_data['Indegree'] + merged_data['Outdegree']
    merged_data['Length'] = merged_data['Target'].str.len()
    
    return merged_data[['Length', 'Degree']]


def plot_side_by_side_scatterplots(filepath1, filepath2, threshold, output_path):
    # Function to prepare data for plotting
    merged_data1 = parse_graph_data(filepath1, threshold)
    merged_data2 = parse_graph_data(filepath2, threshold)
    
    # Compute common x and y limits for both scatterplots
    x_min = min(merged_data1['Length'][merged_data1['Length'] > 0].min(), 
                merged_data2['Length'][merged_data2['Length'] > 0].min())
    y_min = min(merged_data1['Degree'][merged_data1['Degree'] > 0].min(), 
                merged_data2['Degree'][merged_data2['Degree'] > 0].min())

    x_lim = (x_min, max(merged_data1['Length'].max(), merged_data2['Length'].max()))
    y_lim = (y_min, max(merged_data1['Degree'].max(), merged_data2['Degree'].max()))

    # Bin sizes
    x_bins = np.logspace(np.log10(x_lim[0]), np.log10(x_lim[1]), 100)
    y_bins = np.arange(y_lim[0], y_lim[1] + 1, 1)

    # Create a 2x2 subplot grid
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(24, 16))

    # Scatterplots
    ax1.scatter(merged_data1['Length'], merged_data1['Degree'], alpha=0.5)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(x_lim)
    ax1.set_ylim(y_lim)
    #ax1.set_xlabel('Label Length', fontsize=14)
    ax1.set_ylabel('Vertex Degree', fontsize=14)
    ax1.set_title('GRCh38-20-0.10b')
    ax1.grid(True, which="both", ls="--", c='0.7')
    
    ax2.scatter(merged_data2['Length'], merged_data2['Degree'], alpha=0.5)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(x_lim)
    ax2.set_ylim(y_lim)
    #ax2.set_xlabel('Label Length', fontsize=14)
    ax2.set_title('gencode.v40.annotation')
    ax2.grid(True, which="both", ls="--", c='0.7')

    # Heatmaps with normalized colorbars
    weights1 = np.ones_like(merged_data1['Length']) / float(len(merged_data1['Length']))
    cax1 = ax3.hist2d(merged_data1['Length'], merged_data1['Degree'], bins=[x_bins, y_bins], cmap='Blues', 
                      cmin=1/float(len(merged_data1['Length'])), weights=weights1)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlim(x_lim)
    ax3.set_ylim(y_lim)
    ax3.set_xlabel('Label Length', fontsize=14)
    ax3.set_ylabel('Vertex Degree', fontsize=14)
    #ax3.set_title('GRCh38-20-0.10b Heatmap', fontsize=14)
    ax3.grid(True, which="both", ls="--", c='0.7')
    
    weights2 = np.ones_like(merged_data2['Length']) / float(len(merged_data2['Length']))
    cax2 = ax4.hist2d(merged_data2['Length'], merged_data2['Degree'], bins=[x_bins, y_bins], cmap='Blues', 
                      cmin=1/float(len(merged_data2['Length'])), weights=weights2)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlim(x_lim)
    ax4.set_ylim(y_lim)
    ax4.set_xlabel('Label Length', fontsize=14)
    #ax4.set_title('gencode.v40.annotation Heatmap', fontsize=14)
    ax4.grid(True, which="both", ls="--", c='0.7')
    
    # Colorbars for ax3 and ax4 below x-axis
    cbar_ax1 = fig.add_axes([ax3.get_position().x0, ax3.get_position().y0 -0.08, 
                             ax3.get_position().width, 0.02]) # position of the colorbar
    cbar1 = fig.colorbar(cax1[3], cax=cbar_ax1, orientation='horizontal')
    cbar1.set_label('Density', rotation=0)

    cbar_ax2 = fig.add_axes([ax4.get_position().x0, ax4.get_position().y0 -0.08, 
                             ax4.get_position().width, 0.02]) # position of the colorbar
    cbar2 = fig.colorbar(cax2[3], cax=cbar_ax2, orientation='horizontal')
    cbar2.set_label('Density', rotation=0)

    #plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


plot_side_by_side_scatterplots('../input/GRCh38-20-0.10b.fmdg', '../input/gencode.v40.fmdg', 10000000, '../plot/scatter_degree_len.png')
