import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from decimal import Decimal


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


def plot_principal(directory_name, output_filename, scale=None, title=False, xlabel=False, ylabel=False):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Parse all logs from directory
    find_df, _, _, _ = parse_directory_logs(directory_name)

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
    cache_cmap = plt.get_cmap('Dark2', len(caches))
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

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if xlabel:
        plt.xlabel('Average effective query length (bp)', fontsize=14)
    if ylabel:
        plt.ylabel('Queries matched per second', fontsize=14)
    if title:
        plt.title(f"Effective query length vs queries matched per second for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)", fontsize=14)
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png', bbox_inches='tight')
    plt.close()


def plot_fm_gap_cache(directory_name, output_filename, scale=None, title=False, xlabel=False, ylabel=False):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Parse all logs from directory
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

    plt.yscale('log')

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legend
    plt.legend(title="Query Length", loc='upper right', fontsize=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if xlabel:
        plt.xlabel('Cache Depth', fontsize=14)
    if ylabel:
        plt.ylabel('Relative Effective Query Length', fontsize=14)
    if title:
        plt.title(f"Cache Depth vs Relative Effective Query Length for {os.path.basename(directory_name).split('.fmdg')[0]}", fontsize=14)
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png', bbox_inches='tight')
    plt.close()


def plot_fm_gap_permutation(directory_name, output_filename, scale=None, title=False, xlabel=False, ylabel=False):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Parse all logs from directory
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
        subset = subset.sort_values(by='ptime')
        plt.plot(subset['ptime'], subset['FM_gap_ratio'],
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
    plt.legend(title="Query Length", loc='lower right', fontsize=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if xlabel:
        plt.xlabel('Permutation Time (s)', fontsize=14)
    if ylabel:
        plt.ylabel('Relative Effective Query Length', fontsize=14)
    if title:
        plt.title(f"Permutation Time vs Relative Effective Query Length for {os.path.basename(directory_name).split('.fmdg')[0]}", fontsize=14)
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png', bbox_inches='tight')
    plt.close()


def plot_baseline(directory_name, output_filename, scale=None, title=False, xlabel=False, ylabel=False):
    # Check if directory exists
    if not all([os.path.exists(directory_name + suffix) for suffix in ['_baseline_plain', '_baseline_permutation', '_baseline_cache', '_baseline_cache_permutation']]):
        print(f"Baselines for {directory_name} not found!")
        return

    # Parse the four types of benchmarks
    find_df_plain, _, _, _ = parse_directory_logs(directory_name + '_baseline_plain')
    find_df_perm, _, _, _ = parse_directory_logs(directory_name + '_baseline_permutation')
    find_df_cache, _, _, _ = parse_directory_logs(directory_name + '_baseline_cache')
    find_df_cache_permutation, _, _, _ = parse_directory_logs(directory_name + '_baseline_cache_permutation')

    dfs = {
        "plain": find_df_plain,
        "permutation": find_df_perm,
        "cache": find_df_cache,
        "cache+permutation": find_df_cache_permutation
    }

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df_plain['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    benchmark_colors = {
        "plain": "blue",
        "permutation": "green",
        "cache": "red",
        "cache+permutation": "purple"
    }

    # Connect data points of the same query length with regular lines
    for qlength in unique_lengths:
        x_data, y_data = [], []
        for benchmark in ["plain", "permutation", "cache", "cache+permutation"]:
            df = dfs[benchmark]
            subset = df[df['length'] == qlength]
            x_data.extend(subset['number_fork_advances_per_query'].values)
            y_data.extend(subset['queries_per_second'].values)
            plt.loglog(subset['number_fork_advances_per_query'], subset['queries_per_second'], 'o', color=benchmark_colors[benchmark], label=f'{benchmark} {qlength}' if qlength == unique_lengths[0] else "")
        plt.loglog(x_data, y_data, '-', color=length_cmap(unique_lengths.index(qlength)))

    # Connect data points of the same benchmark type with dashed lines
    for benchmark, df in dfs.items():
        sorted_df = df.sort_values(by='length')
        plt.loglog(sorted_df['number_fork_advances_per_query'], sorted_df['queries_per_second'], '--', color=benchmark_colors[benchmark], label=benchmark)

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    # Create the legends with explicit positioning
    benchmark_handles = [plt.Line2D([0], [0], marker='o', color=benchmark_colors[benchmark], linestyle='', label=f'{benchmark}') for benchmark in dfs]
    length_handles = [plt.Line2D([0], [0], color=length_cmap(i), linestyle='-', label=f'{qlength}') for i, qlength in enumerate(unique_lengths)]
    legend1 = plt.legend(handles=length_handles, loc='upper right', title="Query Length")
    plt.gca().add_artist(legend1)  # To ensure first legend is not overwritten by the second
    legend2 = plt.legend(handles=benchmark_handles, loc='upper right', bbox_to_anchor=(0.8750, 1), title="Index")

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if xlabel:
        plt.xlabel('Average effective query length (bp)', fontsize=14)
    if ylabel:
        plt.ylabel('Queries matched per second', fontsize=14)
    if title:
        plt.title(f"Effective query length vs queries matched per second for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)", fontsize=14)
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png', bbox_inches='tight')
    plt.close()



def plot_baseline_simple(directory_name, output_filename, scale=None, title=False, xlabel=False, ylabel=False):
    # Check if directory exists
    if not all([os.path.exists(directory_name + suffix) for suffix in ['_baseline_plain', '_baseline_permutation', '_baseline_cache', '_baseline_cache_permutation']]):
        print(f"Baselines for {directory_name} not found!")
        return

    # Parse the four types of benchmarks
    find_df_plain, _, _, _ = parse_directory_logs(directory_name + '_baseline_plain')
    find_df_perm, _, _, _ = parse_directory_logs(directory_name + '_baseline_permutation')
    find_df_cache, _, _, _ = parse_directory_logs(directory_name + '_baseline_cache')
    find_df_cache_permutation, _, _, _ = parse_directory_logs(directory_name + '_baseline_cache_permutation')

    dfs = {
        "plain": find_df_plain,
        "permutation": find_df_perm,
        "cache": find_df_cache,
        "cache+permutation": find_df_cache_permutation
    }

    plt.figure(figsize=(12, 8))

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df_plain['length'].unique())
    length_cmap = plt.get_cmap('tab10', len(unique_lengths))

    benchmark_colors = {
        "plain": "blue",
        "permutation": "green",
        "cache": "red",
        "cache+permutation": "purple"
    }

    # Plot each benchmark data
    for benchmark, df in dfs.items():
        sorted_df = df.sort_values(by='length')
        plt.loglog(sorted_df['length'], sorted_df['queries_per_second'], 'o-', color=benchmark_colors[benchmark], label=benchmark)

    # If a scale is provided, set the x and y limits
    if scale:
        plt.xlim(scale[0])
        plt.ylim(scale[1])

    plt.legend(loc='upper right', title="Index")

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if xlabel:
        plt.xlabel('Regular query length (bp)', fontsize=14)
    if ylabel:
        plt.ylabel('Queries matched per second', fontsize=14)
    if title:
        plt.title(f"Query length vs queries matched per second for {os.path.basename(directory_name).split('.fmdg')[0]} (log-log scale)", fontsize=14)
    plt.grid(True, which="both", ls="--", c='0.65')
    plt.savefig(output_filename, format='png', bbox_inches='tight')
    plt.close()


def tabulate_decode(directory_name, output_filename):
    # Check if directory exists
    if not os.path.exists(directory_name):
        print(f"Directory {directory_name} not found!")
        return

    # Parse all logs from directory
    find_df, _, _, _ = parse_directory_logs(directory_name)

    # Create a color map for the different query lengths
    unique_lengths = sorted(find_df['length'].unique())

    # get a handle to the output file
    f = open(output_filename, 'w')
    f.write('\\begin{tabular}{|c|c|c|c|}\n\\hline\n Query Length (bp) & Avg. Match Time & Avg. Decode Time & Avg. Count\\\\\\hline\n')

    # Plot the data grouped by query length
    for i, qlength in enumerate(unique_lengths):
        subset = find_df[find_df['length'] == qlength]
        ct = subset['total_matches_decoded'].iloc[0] / subset['total_queries_processed'].iloc[0]
        mdps = subset['total_decoding_time'].iloc[0] / subset['total_queries_processed'].iloc[0]
        mcps = subset['time_per_query'].iloc[0]
        f.write(f'{qlength} & {Decimal(mcps):.3E} & {Decimal(mdps):.3E} & {Decimal(ct):.3E} \\\\ \\hline \n')

    f.write('\\end{tabular}\n')
    f.close()


if __name__ == '__main__':
    principal_scale = ((1e1,5e5),(1,5e5))
    baseline_scale = ((1e1,5e5),(1,5e5))
    fm_gap_cache_scale = ((-1,13),(1e-2,5e5))
    fm_gap_permuatation_scale = ((0,3800),(1e-2,5e5))

    plot_baseline('../log/GRCh38-20-0.10b.fmdg', '../plot/GRCh38-20-0.10b.fmdg_baseline.png', baseline_scale, ylabel=True) #a
    plot_baseline('../log/gencode.v40.fmdg', '../plot/gencode.v40.fmdg_baseline.png', baseline_scale) #b

    plot_baseline_simple('../log/GRCh38-20-0.10b.fmdg', '../plot/GRCh38-20-0.10b.fmdg_baseline_simple.png', baseline_scale, ylabel=True) #?
    plot_baseline_simple('../log/gencode.v40.fmdg', '../plot/gencode.v40.fmdg_baseline_simple.png', baseline_scale) #?

    plot_principal('../log/GRCh38-20-0.10b.fmdg_principal', '../plot/GRCh38-20-0.10b.fmdg_principal.png', principal_scale, xlabel=True, ylabel=True) #c
    plot_principal('../log/gencode.v40.fmdg_principal', '../plot/gencode.v40.fmdg_principal.png', principal_scale, xlabel=True) #d

    plot_fm_gap_permutation('../log/GRCh38-20-0.10b.fmdg_permutation', '../plot/GRCh38-20-0.10b.fmdg_fm_gap_permutation.png', fm_gap_permuatation_scale, xlabel=True, ylabel=True) #e
    plot_fm_gap_permutation('../log/gencode.v40.fmdg_permutation', '../plot/gencode.v40.fmdg_fm_gap_permutation.png', fm_gap_permuatation_scale, xlabel=True,)  #f

    plot_fm_gap_cache('../log/GRCh38-20-0.10b.fmdg_principal', '../plot/GRCh38-20-0.10b.fmdg_fm_gap_cache.png', fm_gap_cache_scale, xlabel=True, ylabel=True) #g
    plot_fm_gap_cache('../log/gencode.v40.fmdg_principal', '../plot/gencode.v40.fmdg_fm_gap_cache.png', fm_gap_cache_scale, xlabel=True, ) #h

    tabulate_decode('../log/gencode.v40.fmdg_rate_decode', '../plot/gencode.v40.fmdg_rate_decode_table.tex')
    tabulate_decode('../log/GRCh38-20-0.10b.fmdg_rate_decode', '../plot/GRCh38-20-0.10b.fmdg_rate_decode_table.tex')
