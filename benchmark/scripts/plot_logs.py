import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def parse_file(filepath):
    basename = os.path.basename(filepath)

    # Common pattern extraction
    common_pattern = re.compile(r'(?P<type>\w+)_log_ptime_(?P<ptime>\d+)_pdepth_(?P<pdepth>\d+)')
    common_fields = common_pattern.match(basename)
    if common_fields is None:
        print(f"Common pattern did not match for filename: {basename}")
        return None
    data = common_fields.groupdict()
    data['ptime'] = int(data['ptime'])
    data['pdepth'] = int(data['pdepth'])

    # Additional type-specific fields
    if basename.startswith('find_'):
        if 'decode' not in basename:
            find_pattern = re.compile(
                r'sampling_rate_(?P<sampling_rate>\d+)_cache_(?P<cache>\d+)_fork_(?P<fork>-?\d+)_match_(?P<match>-?\d+)_threads_(?P<threads>\d+)_length_(?P<length>\d+)\.txt')
            find_fields = find_pattern.search(basename)
            data.update(find_fields.groupdict())
        else:
            find_pattern = re.compile(
                r'sampling_rate_(?P<sampling_rate>\d+)_cache_(?P<cache>\d+)_fork_(?P<fork>-?\d+)_match_(?P<match>-?\d+)_threads_(?P<threads>\d+)_length_(?P<length>\d+)_decode\.txt')
            find_fields = find_pattern.search(basename)
            data.update(find_fields.groupdict())


    elif basename.startswith('perm_'):
        perm_pattern = re.compile(r'threads_(?P<threads>\d+)\.txt')
        perm_fields = perm_pattern.search(basename)
        data.update(perm_fields.groupdict())

    elif basename.startswith('cache_'):
        cache_pattern = re.compile(r'depth_(?P<depth>\d+)\.txt')
        cache_fields = cache_pattern.search(basename)
        data.update(cache_fields.groupdict())

    elif basename.startswith('index_'):
        index_pattern = re.compile(r'sampling_rate_(?P<sampling_rate>\d+)\.txt')
        index_fields = index_pattern.search(basename)
        data.update(index_fields.groupdict())

    # Convert numerical fields to numbers
    for key, value in data.items():
        if key != 'type' and isinstance(value, str) and value.replace('.', '', 1).isdigit():
            data[key] = int(value) if '.' not in value else float(value)

    # Reading the file content
    with open(filepath, 'r') as file:
        for line in file:
            if "Parsing reference index" in line:
                data['parsing_reference'] = True
            if "Constructing a cache" in line:
                data['constructing_cache_depth'] = int(re.search(r'of depth (\d+)', line).group(1))
            if "Dumping cache" in line:
                data['dumping_cache_size'] = float(re.search(r'of size (\d+\.\d+)', line).group(1))
            if "Parsing input file" in line:
                data['input_file'] = re.search(r'Parsing input file (.+)', line).group(1)

            # Parsing type-specific fields
            if basename.startswith('find_'):
                if "Params: Index file name" in line:
                    data['index_file_name'] = re.search(r'Params: Index file name \(-r\): (.+)', line).group(1)
                if "Params: Read batch size" in line:
                    data['read_batch_size'] = int(re.search(r'Params: Read batch size \(-b\): (\d+)', line).group(1))
                if "Params: Threads" in line:
                    data['threads'] = int(re.search(r'Params: Threads \(-j\): (\d+)', line).group(1))
                if "Params: Maximum forks tracked" in line:
                    data['max_forks_tracked'] = int(re.search(r'Params: Maximum forks tracked \(-m\): (\d+|-1)', line).group(1))
                if "Params: Maximum matches decoded" in line:
                    data['max_matches_decoded'] = int(re.search(r'Params: Maximum matches decoded \(-M\): (\d+|-1)', line).group(1))
                if "Params: Cache path:" in line:
                    data['cache_path'] = re.search(r'Params: Cache path: \(-C\): (.+)', line).group(1)
                if "Index: Index parse time" in line:
                    data['index_parse_time'] = float(re.search(r'Index parse time \(s\): (\d+\.\d+)', line).group(1))
                if "Cache: Cache parse time in seconds:" in line:
                    data['cache_parse_time'] = float(re.search(r'Cache parse time in seconds: (\d+\.\d+)', line).group(1))
                if "Cache: Cache size in memory (MB)" in line:
                    data['cache_size_memory_MB'] = float(re.search(r'Cache size in memory \(MB\): (\d+\.\d+)', line).group(1))
                if "Decode: Total matches decoded:" in line:
                    data['total_matches_decoded'] = int(re.search(r'Total matches decoded: (\d+)', line).group(1))
                if "Decode: Total decoding time (s):" in line:
                    data['total_decoding_time'] = float(re.search(r'Total decoding time \(s\): (\d+\.\d+)', line).group(1))
                if "Decode: Matches decoded per second:" in line:
                    data['matches_decoded_per_second'] = float(re.search(r'Matches decoded per second: (\d+\.\d+)', line).group(1))
                if "Decode: Time per match decode (s):" in line:
                    data['time_per_match_decode'] = float(re.search(r'Time per match decode \(s\): (\d+\.\d+)', line).group(1))
                if "Find: Total queries processed" in line:
                    data['total_queries_processed'] = int(re.search(r'Total queries processed: (\d+)', line).group(1))
                if "Find: Total querying time" in line:
                    data['total_querying_time'] = float(re.search(r'Total querying time \(s\): (\d+\.\d+)', line).group(1))
                if "Find: Queries per second:" in line:
                    data['queries_per_second'] = float(re.search(r'Queries per second: (\d+\.\d+)', line).group(1))
                if "Find: Time per query (s)" in line:
                    data['time_per_query'] = float(re.search(r'Time per query \(s\): (\d+\.\d+)', line).group(1))
                if "Find: Number of matching forks" in line:
                    data['number_of_matching_forks'] = int(re.search(r'Number of matching forks: (\d+)', line).group(1))
                if "Find: Number of partial forks" in line:
                    data['number_of_partial_forks'] = int(re.search(r'Number of partial forks: (\d+)', line).group(1))
                if "Find: Number of matches" in line:
                    data['number_of_matches'] = int(re.search(r'Number of matches: (\d+)', line).group(1))
                if "Find: Number of fork advances per query" in line:
                    data['number_fork_advances_per_query'] = float(re.search(r'Number of fork advances per query: (\d+\.\d+)', line).group(1))
                elif "Find: Number of fork advances" in line:
                    data['number_fork_advances'] = int(re.search(r'Number of fork advances: (\d+)', line).group(1))

            elif basename.startswith('index_'):
                if "Number of vertices:" in line:
                    data['vertices'] = int(re.search(r'Number of vertices: (\d+)', line).group(1))
                if "Number of edges:" in line:
                    data['edges'] = int(re.search(r'Number of edges: (\d+)', line).group(1))

            elif basename.startswith('cache_'):
                if "Cache contains" in line:
                    data['cache_contains'] = int(re.search(r'Cache contains (\d+)', line).group(1))

            elif basename.startswith('perm_'):
                if "Optimization begins" in line:
                    data['initial_cost'] = float(re.search(r'initial cost (\d+\.\d+)', line).group(1))

    return data


def parse_directory_logs(dir_path):
    find_data = []
    index_data = []
    cache_data = []
    perm_data = []

    for filename in os.listdir(dir_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(dir_path, filename)
            data = parse_file(file_path)
            if data:
                if filename.startswith('find_'):
                    find_data.append(data)
                elif filename.startswith('index_'):
                    index_data.append(data)
                elif filename.startswith('cache_'):
                    cache_data.append(data)
                elif filename.startswith('perm_'):
                    perm_data.append(data)

    find_df = pd.DataFrame(find_data)
    index_df = pd.DataFrame(index_data)
    cache_df = pd.DataFrame(cache_data)
    perm_df = pd.DataFrame(perm_data)

    return find_df, index_df, cache_df, perm_df


def plot_principal(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Create a new column for total forks
    #find_df['total_forks'] = find_df['number_of_matching_forks'] + find_df['number_of_partial_forks']

    # Compute average forks per query
    #find_df['avg_forks_per_query'] = find_df['total_forks'] / find_df['total_queries_processed']

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Plot the data grouped by query length with regular lines
    for qlength in unique_lengths:
        subset = find_df[find_df['length'] == qlength]
        subset = subset.sort_values(by='cache')
        plt.loglog(subset['number_fork_advances_per_query'], subset['queries_per_second'], '-o', color=length_cmap(unique_lengths.index(qlength)), label=f'{qlength}')

    # Connect data points of the same cache size with dashed lines
    caches = sorted(find_df['cache'].unique())
    cache_cmap = plt.get_cmap('Dark2', len(caches))  # Using the Dark2 colormap for cache lines
    for cache in caches:
        subset = find_df[find_df['cache'] == cache]
        subset = subset.sort_values(by='length')
        plt.loglog(subset['number_fork_advances_per_query'], subset['queries_per_second'], '--', color=cache_cmap(caches.index(cache)))

    # Set the x and y scales to be logarithmic
    plt.xscale('log')
    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legends with explicit positioning
    cache_handles = [plt.Line2D([0], [0], color=cache_cmap(i), linestyle='--', label=f'{cache}') for i, cache in enumerate(caches)]
    legend1 = plt.legend(handles=cache_handles, loc='upper right', title="Cache Size")
    plt.gca().add_artist(legend1)  # To make sure first legend is not overwritten by the second
    legend2 = plt.legend(loc='upper right', bbox_to_anchor=(0.9025, 1), title="Query Length")

    plt.xlabel('Average effective query length')
    plt.ylabel('Queries matched per second')
    plt.title(f"Effective query length vs queries matched per second for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_partial_forks(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)
    find_df['avg_partial_forks_per_query'] = find_df['number_of_partial_forks'] / find_df['total_queries_processed']

    plt.figure(figsize=(12, 8))

    # Unique cache depths to iterate over and assign a color for each one
    cache_depths = sorted(find_df['cache'].unique())
    color_map = plt.get_cmap('tab10', len(cache_depths))

    for cache in cache_depths:
        subset = find_df[find_df['cache'] == cache]
        subset = subset.sort_values(by='length')
        plt.plot(subset['length'], subset['avg_partial_forks_per_query'], '-o', label=f"{cache}", color=color_map(cache_depths.index(cache)))

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    plt.xlabel('Query length')
    plt.ylabel('Average number of partially matching forks')
    plt.title(f"Average number of partially matching forks vs query length for {os.path.basename(directory_name).split('.fmdg')[0]}")
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_threads_rate(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    plt.figure(figsize=(12, 8))

    # Create a list of line styles
    line_styles = ['-', '-.', ':', '--']

    # Create a color map for the different thread values
    unique_threads = sorted(find_df['threads'].unique())
    thread_cmap = plt.get_cmap('tab10', len(unique_threads))

    # Create a marker list for the different sampling rates
    unique_sampling_rates = sorted(find_df['sampling_rate'].unique())
    markers = ['o', 's', 'D', 'X', '^']

    # Plot the data grouped by thread values and sampling rates
    for i, thread in enumerate(unique_threads):
        for j, sampling_rate in enumerate(unique_sampling_rates):
            subset = find_df[(find_df['threads'] == thread) & (find_df['sampling_rate'] == sampling_rate)]
            subset = subset.sort_values(by='length')
            plt.loglog(subset['length'], subset['queries_per_second'],
                       linestyle=line_styles[j % len(line_styles)],
                       marker=markers[i % len(markers)],
                       markersize=4,
                       color=thread_cmap(i / len(unique_threads)),
                       label=f'Threads: {thread}, Sampling Period: {sampling_rate}')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legends with explicit positioning
    thread_handles = [plt.Line2D([0], [0], color=thread_cmap(i / len(unique_threads)), marker=markers[i % len(markers)], linestyle='None', label=f'{thread}') for i, thread in enumerate(unique_threads)]
    legend1 = plt.legend(handles=thread_handles, loc='upper right', title="Threads")
    plt.gca().add_artist(legend1)
    sampling_rate_handles = [plt.Line2D([0], [0], color='black', linestyle=line_styles[i % len(line_styles)], label=f'{sampling_rate}') for i, sampling_rate in enumerate(unique_sampling_rates)]
    legend2 = plt.legend(handles=sampling_rate_handles, loc='upper right', bbox_to_anchor=(0.9025, 1), title="Sampling Period")

    plt.xlabel('Query length')
    plt.ylabel('Queries matched per second')
    plt.title(f"Query length vs queries matched per second for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_threads_rate2(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    plt.figure(figsize=(12, 8))

    # Create a list of line styles
    line_styles = ['-', '-.', ':', '--']

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Create a marker list for the different sampling rates
    unique_sampling_rates = sorted(find_df['sampling_rate'].unique())
    markers = ['o', 's', 'D', 'X', '^']

    # Plot the data grouped by query length and sampling rates
    for i, length in enumerate(unique_lengths):
        for j, sampling_rate in enumerate(unique_sampling_rates):
            subset = find_df[(find_df['length'] == length) & (find_df['sampling_rate'] == sampling_rate)]
            subset = subset.sort_values(by='threads')
            plt.loglog(subset['threads'], subset['queries_per_second'],
                       linestyle=line_styles[j % len(line_styles)],
                       marker=markers[i % len(markers)],
                       markersize=4,
                       color=length_cmap(i / len(unique_lengths)),
                       label=f'Query Length: {length}, Sampling Period: {sampling_rate}')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legends with explicit positioning
    length_handles = [plt.Line2D([0], [0], color=length_cmap(i / len(unique_lengths)), marker=markers[i % len(markers)], linestyle='None', label=f'{length}') for i, length in enumerate(unique_lengths)]
    legend1 = plt.legend(handles=length_handles, loc='upper right', title="Query Length")
    plt.gca().add_artist(legend1)
    sampling_rate_handles = [plt.Line2D([0], [0], color='black', linestyle=line_styles[i % len(line_styles)], label=f'{sampling_rate}') for i, sampling_rate in enumerate(unique_sampling_rates)]
    legend2 = plt.legend(handles=sampling_rate_handles, loc='upper right', bbox_to_anchor=(0.88, 1), title="Sampling Period")

    plt.xlabel('Number of threads')
    plt.ylabel('Queries matched per second')
    plt.title(f"Number of Threads vs Queries per Second for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_decode(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Plot the data grouped by query length
    for i, qlength in enumerate(unique_lengths):
        subset = find_df[find_df['length'] == qlength]
        subset = subset.sort_values(by='sampling_rate')
        plt.loglog(subset['sampling_rate'], subset['matches_decoded_per_second'],
                   linestyle='--',
                   color=length_cmap(unique_lengths.index(qlength)),
                   marker='o',
                   label=f'{qlength}')

    # Set the x and y scales to be logarithmic
    plt.xscale('log')
    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legend
    plt.legend(title="Query Length")

    plt.xlabel('Index sampling period')
    plt.ylabel('Decoding speed (matches per second)')
    plt.title(f"Index Sampling Period vs Decoding Speed for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_decode_diff(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Plot the data grouped by query length
    for i, qlength in enumerate(unique_lengths):
        subset = find_df[find_df['length'] == qlength]
        subset = subset.sort_values(by='sampling_rate')
        plt.loglog(subset['sampling_rate'], subset['matches_decoded_per_second'],
                   linestyle='--',
                   color=length_cmap(unique_lengths.index(qlength)),
                   marker='o',
                   label=f'{qlength}')

    # Set the x and y scales to be logarithmic
    plt.xscale('log')
    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legend
    plt.legend(title="Query Length")

    plt.xlabel('Index sampling period')
    plt.ylabel('Decoding speed (matches per second)')
    plt.title(f"Index Sampling Period vs Decoding Speed for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()

def plot_fm_gap(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Compute FM gap
    find_df['FM_gap'] = find_df['number_fork_advances_per_query'] - find_df['length']

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Plot the data grouped by query length
    for qlength in unique_lengths:
        subset = find_df[find_df['length'] == qlength]
        subset = subset.sort_values(by='cache')
        plt.plot(subset['cache'], subset['FM_gap'],
                 linestyle='-',
                 color=length_cmap(unique_lengths.index(qlength)),
                 marker='o',
                 label=f'{qlength}')

    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legend
    plt.legend(title="Query Length")

    plt.xlabel('Cache Depth')
    plt.ylabel('FM Gap')
    plt.title(f"Cache Depth vs FM Gap for {os.path.basename(directory_name).split('.fmdg')[0]}")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_fm_gap_ratio(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Compute FM gap
    find_df['FM_gap_ratio'] = (find_df['number_fork_advances_per_query'] - find_df['length']) / find_df['length']

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Plot the data grouped by query length
    for qlength in unique_lengths:
        subset = find_df[find_df['length'] == qlength]
        subset = subset.sort_values(by='cache')
        plt.plot(subset['cache'], subset['FM_gap_ratio'],
                 linestyle='-',
                 color=length_cmap(unique_lengths.index(qlength)),
                 marker='o',
                 label=f'{qlength}')

    # Annotate the points for the largest cache depth
    ax = plt.gca()
    fig = plt.gcf()
    fig_size_px = fig.get_size_inches() * fig.dpi
    fig_width_px, fig_height_px = fig_size_px

    # Find the largest cache depth
    max_cache = find_df['cache'].max()
    min_cache = find_df['cache'].min()

    # Get the subset of data for the largest cache depth
    max_cache_subset = find_df[find_df['cache'] == max_cache]
    min_cache_subset = find_df[find_df['cache'] == min_cache]

    # Sort the subset by FM gap ratio
    sorted_subset_first = min_cache_subset.sort_values(by='FM_gap_ratio', ascending=True)
    sorted_subset_last = max_cache_subset.sort_values(by='FM_gap_ratio', ascending=True)

    # Calculate the y_step dynamically
    y_step = 8 / len(sorted_subset_last)

    # Start annotating from the bottom right of the figure
    for i, (index, row) in enumerate(sorted_subset_last.iterrows()):
        x_text = plt.xlim()[1] + 1
        y_text = y_step*(np.exp(i-3.5))
        plt.annotate(f'{row["FM_gap_ratio"] * 100:.2f}%',
                     xy=(row['cache'], row['FM_gap_ratio']),
                     xytext=(x_text, y_text),
                     ha='center',
                     xycoords='data',
                     fontsize=8,
                     arrowprops=dict(arrowstyle="->"))


    # Start annotating from the bottom right of the figure
    for i, (index, row) in enumerate(sorted_subset_first.iterrows()):
        x_text = -2
        y_text = y_step*(np.exp(i))
        plt.annotate(f'{row["FM_gap_ratio"] * 100:.2f}%',
                     xy=(row['cache'], row['FM_gap_ratio']),
                     xytext=(x_text, y_text),
                     ha='center',
                     fontsize=8,
                     xycoords='data',
                     arrowprops=dict(arrowstyle="->"))

    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legend
    plt.legend(title="Query Length")

    plt.xlabel('Cache Depth')
    plt.ylabel('Relative Effective Query Length')
    plt.title(f"Cache Depth vs Relative Effective Query Length for {os.path.basename(directory_name).split('.fmdg')[0]}")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()



if __name__ == '__main__':
    p_c_scale = ((1e1,5e5),(1e1,5e5))
    t_r_scale = ((0,128),(1e2,1e6))
    dec_scale = ((16,64),(1e5,1e7))

    plot_principal('../log/gencode.v40.fmdg_principal', '../plot/gencode.v40.fmdg_principal.png', p_c_scale)
    plot_principal('../log/GRCh38-20-0.10b.fmdg_principal', '../plot/GRCh38-20-0.10b.fmdg_principal.png', p_c_scale)
    plot_fm_gap('../log/gencode.v40.fmdg_principal', '../plot/gencode.v40.fmdg_fm_gap.png')
    plot_fm_gap('../log/GRCh38-20-0.10b.fmdg_principal', '../plot/GRCh38-20-0.10b.fmdg_fm_gap.png')
    plot_fm_gap_ratio('../log/gencode.v40.fmdg_principal', '../plot/gencode.v40.fmdg_fm_gap_ratio.png')
    plot_fm_gap_ratio('../log/GRCh38-20-0.10b.fmdg_principal', '../plot/GRCh38-20-0.10b.fmdg_fm_gap_ratio.png')

    #plot_partial_forks('../log/gencode.v40.fmdg_c_l2_hard', '../plot/gencode.v40.fmdg_partial_forks_hard.png')
    #plot_partial_forks('../log/GRCh38-20-0.10b.fmdg_c_l2_hard', '../plot/GRCh38-20-0.10b.fmdg_partial_forks_hard.png')

    #plot_threads_rate2('../log/gencode.v40.fmdg_j_l4_hard', '../plot/gencode.v40.fmdg_threads_rates_hard.png', t_r_scale)
    #plot_threads_rate2('../log/GRCh38-20-0.10b.fmdg_j_l4_hard', '../plot/GRCh38-20-0.10b.fmdg_threads_rates_hard.png', t_r_scale)

    #plot_decode('../log/gencode.v40.fmdg_decode3_hard', '../plot/gencode.v40.fmdg_decode_hard3.png', dec_scale)
    #plot_decode('../log/GRCh38-20-0.10b.fmdg_decode3_hard', '../plot/GRCh38-20-0.10b.fmdg_decode_hard3.png', dec_scale)

    #plot_decode('../log/gencode.v40.fmdg_decode4_hard', '../plot/gencode.v40.fmdg_decode_hard.png', dec_scale)
    #plot_decode('../log/GRCh38-20-0.10b.fmdg_decode4_hard', '../plot/GRCh38-20-0.10b.fmdg_decode_hard.png', dec_scale)

    #plot_decode_diff('../log/gencode.v40.fmdg_decode4_hard', '../plot/gencode.v40.fmdg_decode_diff_hard.png')
    #plot_decode_diff('../log/GRCh38-20-0.10b.fmdg_decode4_hard', '../plot/GRCh38-20-0.10b.fmdg_decode_diff_hard.png')
