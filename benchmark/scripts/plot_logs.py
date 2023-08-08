import os
import pandas as pd
import glob
import matplotlib.pyplot as plt

def parse_log_files(directory):
    log_files = glob.glob(os.path.join(directory, "*.txt"))
    all_data = []

    for log_file in log_files:
        with open(log_file, 'r') as file:
            file_name = os.path.splitext(os.path.basename(log_file))[0]
            data = {"filename": (os.path.splitext(os.path.basename(log_file))[0]).split("_")[0]}
            for line in file:
                if "Params:" in line:
                    vals = line.split("Params:")[1].split(":")
                    if(len(vals) == 2):
                        key, value = vals
                        data[key.strip()] = value.strip()
                    else:
                        key, _, value = vals
                        data[key.strip()] = value.strip()
                elif "Index:" in line:
                    key, value = line.split("Index:")[1].split(":")
                    data[key.strip()] = value.strip()
                elif "Find:" in line:
                    key, value = line.split("Find:")[1].split(":")
                    data[key.strip()] = value.strip()

            data['Query length'] = int(file_name.split('_')[-1])
            all_data.append(data)

    retval = pd.DataFrame(all_data)
    retval['Query length'] = pd.to_numeric(retval['Query length'], errors='coerce')
    retval['Queries per second'] = pd.to_numeric(retval['Queries per second'], errors='coerce')
    retval['Maximum forks tracked (-m)'] = pd.to_numeric(retval['Maximum forks tracked (-m)'], errors='coerce')
    retval['Number of matches'] = pd.to_numeric(retval['Number of matches'], errors='coerce')
    return retval

def plot_query_speed_vs_ml(df):
    unique_filenames = df['filename'].unique()

    for filename in unique_filenames:
        subset = df[df['filename'] == filename]
        unique_m_values = sorted(subset['Maximum forks tracked (-m)'].unique())
        plt.figure(figsize=(10, 8))
        for m in unique_m_values:
            subset_m = subset[subset['Maximum forks tracked (-m)'] == m]
            subset_m = subset_m.sort_values(by='Query length')
            plt.plot(subset_m['Query length'], subset_m['Queries per second'], label=f'-m {m}')
        plt.legend()
        plt.xlabel('Query Length')
        plt.ylabel('Matches per Second')
        plt.title(f'Query Length vs Matches per Second for {filename}')
        plt.savefig(f'{filename}_query_length_vs_matches.png')

def plot_query_speed_vs_ml_log(df):
    df['Query length'] = pd.to_numeric(df['Query length'], errors='coerce')
    df['Queries per second'] = pd.to_numeric(df['Queries per second'], errors='coerce')

    unique_filenames = df['filename'].unique()
    for filename in unique_filenames:
        subset = df[df['filename'] == filename]
        unique_m_values = sorted(subset['Maximum forks tracked (-m)'].unique())

        plt.figure(figsize=(10, 8))
        for m in unique_m_values:
            subset_m = subset[subset['Maximum forks tracked (-m)'] == m]
            subset_m = subset_m.sort_values(by='Query length')
            plt.loglog(subset_m['Query length'], subset_m['Queries per second'], label=f'-m {m}')

        plt.legend()
        plt.xlabel('Query Length')
        plt.ylabel('Matches per Second')
        plt.title(f'Log-Log Plot: Query Length vs Matches per Second for {filename}')
        plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.tight_layout()
        plt.savefig(f'{filename}_query_length_vs_matches_loglog.png')

def plot_matches_returned_ml_log(df):
    unique_filenames = df['filename'].unique()

    for filename in unique_filenames:
        subset = df[df['filename'] == filename]
        unique_query_lengths = sorted(subset['Query length'].unique())

        plt.figure(figsize=(10, 8))
        for ql in unique_query_lengths:
            subset_ql = subset[subset['Query length'] == ql]
            subset_ql = subset_ql.sort_values(by='Maximum forks tracked (-m)')
            original_matches = subset_ql[subset_ql['Maximum forks tracked (-m)'] == -1]['Number of matches']
            subset_ql = subset_ql[subset_ql['Maximum forks tracked (-m)'] != -1]
            if original_matches.empty:
                continue
            subset_ql['Percentage of queries returned'] = (subset_ql['Number of matches'] / original_matches.iloc[0]) * 100
            plt.plot(subset_ql['Maximum forks tracked (-m)'], subset_ql['Percentage of queries returned'], label=f'Query Length {ql}')
        plt.legend()
        plt.xlabel('Maximum Forks Tracked (-m)')
        plt.ylabel('Percentage of Queries Returned')
        plt.title(f'Percentage of Queries Returned vs -m for {filename}')
        plt.savefig(f'{filename}_percentage_of_queries_vs_m.png')


log_directory = "../log/query_speed_vs_ml_64"
df = parse_log_files(log_directory)

plot_query_speed_vs_ml(df)
plot_query_speed_vs_ml_log(df)
plot_matches_returned_ml_log(df)