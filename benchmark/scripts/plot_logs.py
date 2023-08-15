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


def plot_c_l(directory_name, output_path):
    # Step 1: Parsing the logs to get the find_df
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # We're assuming 'qps' is already in the find_df, if not we'll calculate it.
    if 'qps' not in find_df.columns:
        find_df['qps'] = find_df['total_queries_processed'] / find_df['total_querying_time']

    # Step 2: Grouping by cache size
    grouped = find_df.groupby('cache')

    # Step 3: Plotting
    plt.figure(figsize=(10, 6))
    colors = plt.cm.jet(np.linspace(0, 1, len(grouped)))

    for (cache, group), color in zip(grouped, colors):
        group = group.sort_values(by='length')
        plt.loglog(group['length'], group['qps'], label=f'Cache Size: {cache}', color=color)

    plt.xlabel('Query Length (log scale)')
    plt.ylabel('Queries per Second (log scale)')
    directory_basename = os.path.basename(directory_name)
    plt.title(f'Effect of cache size on query speed for {directory_basename}')
    plt.legend()
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.savefig(output_path, format='png', dpi=300)
    plt.close()


def plot_m_l(directory_name, output_filename):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Parse logs in the directory
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Start plotting
    plt.figure(figsize=(10, 6))

    # Group by query length
    grouped = find_df.groupby('length')

    for name, group in grouped:
        # Find the total matches for m=-1
        total_matches = group[group['max_forks_tracked'] == -1]['number_of_matches'].values[0]
        # Sort the group by max_forks_tracked and filter out m=-1
        sorted_group = group[group['max_forks_tracked'] != -1].sort_values(by='max_forks_tracked')
        # Calculate the percentage of matches
        sorted_group['percentage_matches'] = 100 * sorted_group['number_of_matches'] / total_matches
        plt.plot(sorted_group['max_forks_tracked'], sorted_group['percentage_matches'], label=f'Length {name}', marker='o')

    plt.xscale('log', base=2)
    plt.xlabel('Max Forks Tracked (m)')
    plt.ylabel('Percentage of Matches Returned')
    plt.title(f"Percentage of Matches vs. Max Forks Tracked for {os.path.basename(directory_name)}")
    plt.legend()

    # Handle tick formatting
    ax = plt.gca()
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.get_yaxis().set_major_formatter(plt.ScalarFormatter())
    ax.xaxis.set_major_locator(plt.FixedLocator(sorted_group['max_forks_tracked'].unique()))
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_p_f(directory_name, output_filename):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Create a new column for total forks
    find_df['total_forks'] = find_df['number_of_matching_forks'] + find_df['number_of_partial_forks']

    plt.figure(figsize=(10, 6))

    # For each unique cache size
    for cache in sorted(find_df['cache'].unique()):
        subset = find_df[find_df['cache'] == cache]

        # Sort this subset by query length
        subset = subset.sort_values(by='length')

        plt.plot(subset['length'], subset['total_forks'], label=f'Cache Size: {cache}', marker='o')

    plt.xlabel('Query Length')
    plt.ylabel('Total Forks')
    plt.title(f"Query Length vs Total Forks for {os.path.basename(directory_name)}")
    plt.legend()
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()


def plot_p_c(directory_name, output_filename):
    # Assuming the first return value from your function is the DataFrame
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Create a new column for total forks
    find_df['total_forks'] = find_df['number_of_matching_forks'] + find_df['number_of_partial_forks']

    # Compute average forks per query
    find_df['avg_forks_per_query'] = find_df['total_forks'] / find_df['total_queries_processed']

    plt.figure(figsize=(12, 8))

    # First plot the data grouped by cache
    caches = sorted(find_df['cache'].unique())
    cache_handles = []
    for cache in caches:
        subset = find_df[find_df['cache'] == cache]
        subset = subset.sort_values(by='length')

        line, = plt.loglog(subset['avg_forks_per_query'], subset['queries_per_second'], marker='o')
        cache_handles.append(line)

    # Create a color map for the different query lengths, using a colormap that trends towards darker colors
    unique_lengths = sorted(find_df['length'].unique())
    cmap = plt.get_cmap('Dark2', len(unique_lengths))

    # Draw dashed lines connecting points of the same query length across cache groups
    lines_for_legend = []  # This will store handles for the legend
    for idx, qlength in enumerate(unique_lengths):
        line_added_to_legend = False  # Track if we added the line to the legend yet

        for i in range(len(caches) - 1):  # We subtract 1 because we'll be using i and i+1 as indices
            subset1 = find_df[(find_df['cache'] == caches[i]) & (find_df['length'] == qlength)]
            subset2 = find_df[(find_df['cache'] == caches[i+1]) & (find_df['length'] == qlength)]

            # Ensure there's a point in both caches to connect
            if not subset1.empty and not subset2.empty:
                line, = plt.plot([subset1['avg_forks_per_query'].values[0], subset2['avg_forks_per_query'].values[0]],
                                 [subset1['queries_per_second'].values[0], subset2['queries_per_second'].values[0]],
                                 '--',
                                 color=cmap(idx))

                # Add line to the legend only once per query length
                if not line_added_to_legend:
                    lines_for_legend.append(line)
                    line_added_to_legend = True

    # Create the legends with explicit positioning
    legend1 = plt.legend(cache_handles, [str(c) for c in caches], loc='upper right', title="Cache Size", bbox_to_anchor=(1, 1))
    legend2 = plt.legend(lines_for_legend, unique_lengths, loc='upper right', bbox_to_anchor=(1, 0.75), title="Query Length")

    plt.gca().add_artist(legend1)  # To make sure first legend is not overwritten by the second

    plt.xlabel('Average Forks per Query (log scale)')
    plt.ylabel('Queries per Second (log scale)')
    plt.title(f"Average Forks per Query vs Queries per Second for {os.path.basename(directory_name)} (log-log scale)")
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png')
    plt.close()



if __name__ == '__main__':
    plot_c_l('../log/gencode.v40.fmdg_c_l', '../plot/gencode.v40.fmdg_c_l.png')
    plot_c_l('../log/GRCh38-20-0.10b.fmdg_c_l', '../plot/GRCh38-20-0.10b.fmdg_c_l.png')

    plot_m_l('../log/gencode.v40.fmdg_m_l', '../plot/gencode.v40.fmdg_m_l.png')
    plot_m_l('../log/GRCh38-20-0.10b.fmdg_m_l', '../plot/GRCh38-20-0.10b.fmdg_m_l.png')

    plot_p_f('../log/gencode.v40.fmdg_c_l', '../plot/gencode.v40.fmdg_p_f.png')
    plot_p_f('../log/GRCh38-20-0.10b.fmdg_c_l', '../plot/GRCh38-20-0.10b.fmdg_p_f.png')

    plot_p_c('../log/gencode.v40.fmdg_c_l', '../plot/gencode.v40.fmdg_p_c.png')
    plot_p_c('../log/GRCh38-20-0.10b.fmdg_c_l', '../plot/GRCh38-20-0.10b.fmdg_p_c.png')
