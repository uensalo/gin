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
        find_pattern = re.compile(
            r'sampling_rate_(?P<sampling_rate>\d+)_cache_(?P<cache>\d+)_fork_(?P<fork>-?\d+)_match_(?P<match>-?\d+)_threads_(?P<threads>\d+)_length_(?P<length>\d+)\.txt')
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
    find_df['total_forks'] = find_df['number_of_matching_forks'] + find_df['number_of_partial_forks']

    # Compute average forks per query
    find_df['avg_forks_per_query'] = find_df['total_forks'] / find_df['total_queries_processed']

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    # Plot the data grouped by query length with regular lines
    for qlength in unique_lengths:
        subset = find_df[find_df['length'] == qlength]
        subset = subset.sort_values(by='cache')
        plt.loglog(subset['avg_forks_per_query'], subset['queries_per_second'], '-o', color=length_cmap(unique_lengths.index(qlength)), label=f'{qlength}')

    # Connect data points of the same cache size with dashed lines
    caches = sorted(find_df['cache'].unique())
    cache_cmap = plt.get_cmap('Dark2', len(caches))  # Using the Dark2 colormap for cache lines
    for cache in caches:
        subset = find_df[find_df['cache'] == cache]
        subset = subset.sort_values(by='length')
        plt.loglog(subset['avg_forks_per_query'], subset['queries_per_second'], '--', color=cache_cmap(caches.index(cache)))

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

    plt.xlabel('Average Forks per Query (log scale)')
    plt.ylabel('Queries per Second (log scale)')
    plt.title(f"Average Forks per Query vs Queries per Second for {os.path.basename(directory_name)} (log-log scale)")
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

    plt.xlabel('Query Length')
    plt.ylabel('Average Number of Forks')
    plt.title(f"Average Number of Forks vs Query Length from {os.path.basename(directory_name)}")
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_threads(directory_name, output_filename, scale=None):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Create a new column for total forks
    find_df['total_forks'] = find_df['number_of_matching_forks'] + find_df['number_of_partial_forks']

    # Compute average forks per query
    find_df['avg_forks_per_query'] = find_df['total_forks'] / find_df['total_queries_processed']

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    #unique_lengths = sorted(find_df['length'].unique())
    #length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    unique_threads = sorted(find_df['threads'].unique())
    threads_cmap = plt.get_cmap('tab10', len(unique_threads))

    # Plot the data grouped by query length with regular lines
    for threads in unique_threads:
        subset = find_df[find_df['threads'] == threads]
        subset = subset.sort_values(by='length')
        plt.loglog(subset['length'], subset['queries_per_second'], '-o', color=threads_cmap(unique_threads.index(threads)), label=f'{threads}')

    # Set the x and y scales to be logarithmic
    plt.xscale('log')
    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    plt.legend(loc='upper right', title='No threads')
    plt.xlabel('Query Length (log scale)')
    plt.ylabel('Queries per Second (log scale)')
    plt.title(f"Query Length vs. Queries per Second {os.path.basename(directory_name)} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
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

    # Create a color map for the different thread values
    unique_threads = sorted(find_df['threads'].unique())
    thread_cmap = plt.get_cmap('tab10', len(unique_threads))

    # Create a marker list for different sampling rates
    unique_rates = sorted(find_df['sampling_rate'].unique())
    markers = ['o', 's', 'D', '^', 'v', '<', '>', 'P', 'X']

    # Plot the data grouped by thread values and sampling rates
    for thread in unique_threads:
        for rate in unique_rates:
            subset = find_df[(find_df['threads'] == thread) & (find_df['sampling_rate'] == rate)]
            subset = subset.sort_values(by='length')
            plt.loglog(subset['length'], subset['queries_per_second'], '-o', color=thread_cmap(unique_threads.index(thread)), marker=markers[unique_rates.index(rate)], label=f'Threads: {thread}, Sampling Rate: {rate}')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legend
    plt.legend(title="Threads and Sampling Rate")

    plt.xlabel('Query Length')
    plt.ylabel('Queries per Second')
    plt.title(f"Query Length vs Queries per Second for {os.path.basename(directory_name)}")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()

# Example usage:
# plot_threads_rate('path/to/directory', 'output.png')



if __name__ == '__main__':
    p_c_scale = ((5,2e5),(5,5e5))

    plot_principal('../log/gencode.v40.fmdg_c_l2_hard', '../plot/gencode.v40.fmdg_p_c2_hard.png', p_c_scale)
    plot_principal('../log/GRCh38-20-0.10b.fmdg_c_l2_hard', '../plot/GRCh38-20-0.10b.fmdg_p_c2_hard.png', p_c_scale)

    plot_partial_forks('../log/gencode.v40.fmdg_c_l2_hard', '../plot/gencode.v40.fmdg_partial_forks_hard.png')
    plot_partial_forks('../log/GRCh38-20-0.10b.fmdg_c_l2_hard', '../plot/GRCh38-20-0.10b.fmdg_partial_forks_hard.png')

    plot_threads('../log/gencode.v40.fmdg_j_l2_hard', '../plot/gencode.v40.fmdg_threads_hard.png')
    plot_threads('../log/GRCh38-20-0.10b.fmdg_j_l2_hard', '../plot/GRCh38-20-0.10b.fmdg_threads_hard.png')

    plot_threads_rate('../log/gencode.v40.fmdg_j_l3_hard', '../plot/gencode.v40.fmdg_threads_rates_hard.png')
    plot_threads_rate('../log/GRCh38-20-0.10b.fmdg_j_l3_hard', '../plot/GRCh38-20-0.10b.fmdg_threads_rates_hard.png')
