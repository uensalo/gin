import os
import pandas as pd
import glob
import matplotlib.pyplot as plt

def parse_log_files(directory):
    log_files = glob.glob(os.path.join(directory, "*.log"))
    all_data = []

    for log_file in log_files:
        with open(log_file, 'r') as file:
            data = {"filename": os.path.splitext(os.path.basename(log_file))[0]}
            for line in file:
                if "Params:" in line:
                    key, value = line.split("Params:")[1].split(":")
                    data[key.strip()] = value.strip()
                elif "Index:" in line:
                    key, value = line.split("Index:")[1].split(":")
                    data[key.strip()] = value.strip()
                elif "Find:" in line:
                    key, value = line.split("Find:")[1].split(":")
                    data[key.strip()] = value.strip()

            data['Query length'] = int(data['filename'].split('_')[-1])
            all_data.append(data)

    return pd.DataFrame(all_data)

def plot_query_speed_vs_ml(df):
    # Extract the unique filenames, -m values and query lengths
    unique_filenames = df['filename'].unique()

    for filename in unique_filenames:
        subset = df[df['filename'] == filename]
        unique_m_values = subset['Maximum forks tracked (-m)'].unique()
        unique_query_lengths = subset['Query length'].unique()

        # Plot query length vs matches per second for each -m value
        plt.figure(figsize=(10, 8))
        for m in unique_m_values:
            subset_m = subset[subset['Maximum forks tracked (-m)'] == m]
            plt.plot(subset_m['Query length'], subset_m['Queries per second'], label=f'-m {m}')
        plt.legend()
        plt.xlabel('Query Length')
        plt.ylabel('Matches per Second')
        plt.title(f'Query Length vs Matches per Second for {filename}')
        plt.savefig(f'{filename}_query_length_vs_matches.png')

        # Plot percentage of queries returned vs -m for each query length
        plt.figure(figsize=(10, 8))
        for ql in unique_query_lengths:
            subset_ql = subset[subset['Query length'] == ql]
            original_matches = subset_ql[subset_ql['Maximum forks tracked (-m)'] == -1]['Number of matches']
            subset_ql['Percentage of queries returned'] = subset_ql['Number of matches'] / original_matches
            plt.plot(subset_ql['Maximum forks tracked (-m)'], subset_ql['Percentage of queries returned'], label=f'Query Length {ql}')
        plt.legend()
        plt.xlabel('Maximum Forks Tracked (-m)')
        plt.ylabel('Percentage of Queries Returned')
        plt.title(f'Percentage of Queries Returned vs -m for {filename}')
        plt.savefig(f'{filename}_percentage_of_queries_vs_m.png')


log_directory = "../log/query_speed_vs_ml"
df = parse_log_files(log_directory)
plot_query_speed_vs_ml(df)


