import os
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import lambertw

class LogData:
    def __init__(self):
        self.data = {}

    def __str__(self):
        return str(self.data)

def extract_value(pattern, content, default=None):
    match = re.search(pattern, content)
    return match.group(1) if match else default

def parse_file(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    data = LogData()

    fmd_index_flag = False
    parsing_timings_flag = False
    parsing_permutation = False

    for line in lines:

        if "[fmd:" in line:
            prefix, content = line.strip().split("]", 1)
            prefix = prefix.strip("[")

        if fmd_index_flag:
            if "Number of vertices:" in line:
                data.data["Vertices"] = int(extract_value(r'Number of vertices: (\d+)', line))
            elif "Number of edges:" in line:
                data.data["Edges"] = int(extract_value(r'Number of edges: (\d+)', line))
            elif "Average in-degree:" in line:
                data.data["AverageInDegree"] = float(extract_value(r'Average in-degree: (.+)', line))
            elif "Total length of vertex labels:" in line:
                data.data["VertexLabels"] = int(extract_value(r'Total length of vertex labels: (\d+)', line))
            elif "Average label length per vertex:" in line:
                data.data["AvgLabelLength"] = float(extract_value(r'Average label length per vertex: (.+)', line))
                fmd_index_flag = False  # reset the flag
            continue

        if parsing_timings_flag:
            if "Parsing:" in line:
                data.data["ParsingTime"] = float(extract_value(r'Parsing:  (.+)', line))
            elif "Indexing:" in line:
                data.data["IndexingTime"] = float(extract_value(r'Indexing: (.+)', line))
            elif "Writing:" in line:
                data.data["WritingTime"] = float(extract_value(r'Writing:  (.+)', line))
                parsing_timings_flag = False  # reset the flag
            continue

        if parsing_permutation:
            data.data["Iterations"] = int(re.search(r'Iteration (\d+)', line).group(1))
            data.data["BestCost"]  = float(re.search(r'best_cost = (\d+\.\d+)', line).group(1))
            data.data["CurCost"]  = float(re.search(r'cur_cost = (\d+\.\d+)', line).group(1))
            parsing_permutation = False
            continue

        if prefix == "fmd:permutation" and "Optimization begins" in content:
            data.data["InitialCost"] = float(extract_value(r'initial cost (.+?) for', content))
            data.data["PermutationDepth"] = int(extract_value(r'depth=(\d+)', content))
            data.data["PermutationTime"] = int(extract_value(r'for (\d+) seconds', content))
            parsing_permutation = True

        elif prefix == "fmd:index":
            if "Constructing FMD" in content:
                fmd_index_flag = True
                continue
            elif "Timings in seconds" in content:
                parsing_timings_flag = True
                continue
            elif "Resulting index size" in content:
                data.data["IndexSize"] = float(extract_value(r'Resulting index size: (\d+)', content))

            #if "Estimated peak memory" in content:
            #    data.data["MemoryRequirement"] = float(extract_value(r'Estimated peak memory requirement: (\d+)', content))
            #    data.data["IndexSize"] = float(extract_value(r'Resulting index size: (\d+) bytes', content))

        elif prefix == "fmd:query":
            if "Parsing and launching queries with" in content:
                data.data["QueryThreads"] = int(extract_value(r'with (\d+) threads', line))
                data.data["BatchSize"] = int(extract_value(r'and (\d+) batch size', line))
            elif "Total querying time in seconds:" in content:
                data.data["TotalQueryTime"] = float(extract_value(r'Total querying time in seconds: (.+)', line))
            elif "Queries processed:" in content:
                data.data["QueriesProcessed"] = int(extract_value(r'Queries processed: (.+)', line))
            elif "Average time per query:" in content:
                data.data["AvgQueryTime"] = float(extract_value(r'Average time per query: (.+)', line))
            elif "Total forks for matching queries:" in content:
                data.data["TotalForksMatching"] = int(extract_value(r'Total forks for matching queries: (.+)', line))
            elif "Total count for matching queries:" in content:
                data.data["TotalCountMatching"] = int(extract_value(r'Total count for matching queries: (.+)', line))
            elif "Total forks for partially matching queries:" in content:
                data.data["TotalForksPartialMatching"] = int(extract_value(r'Total forks for partially matching queries: (.+)', line))
            elif "Average forks per matching:" in content:
                data.data["AvgForksPerMatch"] = float(extract_value(r'Average forks per matching: (.+)', line))
            elif "Average forks per query:" in content:
                data.data["AvgForksPerQuery"] = float(extract_value(r'Average forks per query: (.+)', line))
            elif "Time per fork:" in content:
                data.data["TimePerFork"] = float(extract_value(r'Time per fork: (.+)', line))
            elif "Time per fork per thread:" in content:
                data.data["TimePerForkPerThread"] = float(extract_value(r'Time per fork per thread: (.+)', line))

    return data


def parse_directory(directory):
    data = {}
    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            filepath = os.path.join(directory, filename)
            match = re.search(r'query_length_(\d+)', filename)
            query_length = int(match.group(1)) if match else 8
            match = re.search(r'sample_rate_(\d+)', filename)
            sample_rate = int(match.group(1)) if match else 64
            data[filename] = parse_file(filepath)
            data[filename].data["QueryLength"] = query_length
            data[filename].data["SampleRate"] = sample_rate
    return data

def flatten_and_sort(data, sort_by='QueryThreads'):
    # Flatten the data
    flat_data = []
    for filename, filedata in data.items():
        filedata.data['filename'] = filename
        flat_data.append(filedata.data)

    # Sort the data
    sorted_data = sorted(flat_data, key=lambda x: x[sort_by])

    return sorted_data

def plot_data(x, y, xlabel='x', ylabel='y', title='x vs. y', filename=None, func=None):
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y)

    if func is not None:
        # perform the curve fit
        params, params_covariance = curve_fit(func, x, y)

        # create a range of x values from the minimum and maximum x values provided
        x_range = np.linspace(min(x), max(x), 500)

        # compute the y values for these x values given the fitted parameters
        y_range = func(x_range, *params)

        # plot the fitted function
        plt.plot(x_range, y_range, 'r-', label=f'fit: {params}')

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Scatter Plot')
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    #plt.show()
    if(filename is not None):
        plt.savefig(filename)


###############################################################################
#   Query Speed vs Query Length Plots
###############################################################################
query_speed_vs_query_length_data = parse_directory('../log/query_speed_vs_query_length')
query_speed_vs_query_length_data_flat = flatten_and_sort(query_speed_vs_query_length_data, 'QueryLength')
for item in query_speed_vs_query_length_data_flat:
    print(item)

#######################################
#   Total forks vs Query Length
#######################################
p1x = np.array([item['QueryLength'] for item in query_speed_vs_query_length_data_flat])
p1y = np.array([item['TotalForksMatching']+item['TotalForksPartialMatching'] for item in query_speed_vs_query_length_data_flat])
def curve1(x,a,b,c,d):
    return a*np.log(x)+c
plot_data(p1x,p1y, 'Query Length (bp)', 'Number of forks explored', '', None, curve1)

#######################################
#   Match return speed vs query length
#######################################
p2x = np.array([item['QueryLength'] for item in query_speed_vs_query_length_data_flat])
p2y = np.array([item['TotalForksMatching']/item['TotalQueryTime'] for item in query_speed_vs_query_length_data_flat])
plot_data(p2x,p2y,'Query length', 'Matches returned per second')

#######################################
#   Fork explore speed vs query length
#######################################
p2x = np.array([item['QueryLength'] for item in query_speed_vs_query_length_data_flat])
p2y = np.array([(item['TotalForksMatching']+item['TotalForksPartialMatching'])/item['TotalQueryTime'] for item in query_speed_vs_query_length_data_flat])
plot_data(p2x,p2y,'Query length', 'Forks explored per second')

#######################################
#   Query speed vs query length
#######################################
p3x = np.array([item['QueryLength'] for item in query_speed_vs_query_length_data_flat])
p3y = np.array([item['QueriesProcessed']/item['TotalQueryTime'] for item in query_speed_vs_query_length_data_flat])
plot_data(p3x,p3y,'Query length', 'Queries processed per second')

#######################################
#   Limited matches decoded
#######################################
#p4x = np.array([item['QueryLength'] for item in query_speed_vs_query_length_data_flat])
#p4y = np.array([item['QueriesProcessed']*50/item['TotalForksMatching']/item['TotalQueryTime'] for item in query_speed_vs_query_length_data_flat])
#plot_data(p4x,p4y,'Query length', 'Queries processed per second if queries limited to 50 matches')

###############################################################################
#   Query Speed vs Number of threads
###############################################################################
query_speed_vs_number_of_threads = parse_directory('../log/query_speed_vs_num_threads')
query_speed_vs_number_of_threads_flat = flatten_and_sort(query_speed_vs_number_of_threads, 'QueryThreads')
for item in query_speed_vs_query_length_data_flat:
    print(item)

#######################################
#   Match speed vs query length
#######################################
p2x = np.array([item['QueryThreads'] for item in query_speed_vs_number_of_threads_flat])
p2y = np.array([item['TotalForksMatching']/item['TotalQueryTime'] for item in query_speed_vs_number_of_threads_flat])
plot_data(p2x,p2y,'Number of threads', 'Matches returned per second')


#######################################
#   Query speed vs query length
#######################################
p3x = np.array([item['QueryThreads'] for item in query_speed_vs_number_of_threads_flat])
p3y = np.array([item['QueriesProcessed']/item['TotalQueryTime'] for item in query_speed_vs_number_of_threads_flat])
plot_data(p3x,p3y,'Number of threads', 'Queries processed per second')


###############################################################################
#   Query Speed vs Permutation time
###############################################################################
query_speed_vs_permutation_time = parse_directory('../log/query_speed_vs_permutation_time')
query_speed_vs_permutation_time_flat = flatten_and_sort(query_speed_vs_permutation_time, 'PermutationTime')
for item in query_speed_vs_permutation_time_flat:
    print(item)

#######################################
#   Matching forks vs Permutation time
#######################################
p2x = np.array([item['PermutationTime'] for item in query_speed_vs_permutation_time_flat])
p2y = np.array([item['TotalForksMatching']/item['QueriesProcessed'] for item in query_speed_vs_permutation_time_flat])
plot_data(p2x,p2y,'Permutation time', 'Average number of matching forks per query')


#######################################
#   Total forks vs Permutation time
#######################################
p2x = np.array([item['PermutationTime'] for item in query_speed_vs_permutation_time_flat])
p2y = np.array([(item['TotalForksMatching']+item['TotalForksPartialMatching'])/item['QueriesProcessed'] for item in query_speed_vs_permutation_time_flat])
plot_data(p2x,p2y,'Permutation time', 'Average number of forks per query')

#######################################
#   Final initial call vs speedup
#######################################
first_time = query_speed_vs_permutation_time_flat[0]['AvgQueryTime']
p2x = np.array([item['InitialCost'] / item['BestCost'] for item in query_speed_vs_permutation_time_flat])
p2y = np.array([first_time / item['AvgQueryTime'] for item in query_speed_vs_permutation_time_flat])
plot_data(p2x,p2y,'Initial Cost / Final Cost', 'Average query speedup')

#######################################
#   Query speed vs permutation time
#######################################
first_time = query_speed_vs_permutation_time_flat[0]['AvgQueryTime']
p2x = np.array([item['PermutationTime'] for item in query_speed_vs_permutation_time_flat])
p2y = np.array([item['QueriesProcessed']/item['TotalQueryTime'] for item in query_speed_vs_permutation_time_flat])
plot_data(p2x,p2y,'Permutation time', 'Queries Per Second')

#######################################
#   Match speed vs permutation time
#######################################
first_time = query_speed_vs_permutation_time_flat[0]['AvgQueryTime']
p2x = np.array([item['PermutationTime'] for item in query_speed_vs_permutation_time_flat])
p2y = np.array([item['TotalForksMatching']/item['TotalQueryTime'] for item in query_speed_vs_permutation_time_flat])
plot_data(p2x,p2y,'Permutation Time', 'Matches per second')


###############################################################################
#   Query Speed vs Permutation depth
###############################################################################
query_speed_vs_permutation_depth = parse_directory('../log/query_speed_vs_permutation_depth')
query_speed_vs_permutation_depth_flat = flatten_and_sort(query_speed_vs_permutation_depth, 'PermutationDepth')
for item in query_speed_vs_permutation_depth_flat:
    print(item)

#######################################
#   Matching forks vs Permutation time
#######################################
p2x = np.array([item['PermutationDepth'] for item in query_speed_vs_permutation_depth_flat])
p2y = np.array([item['TotalForksMatching']/item['QueriesProcessed'] for item in query_speed_vs_permutation_depth_flat])
plot_data(p2x,p2y,'Permutation Depth', 'Average number of matching forks per query')


#######################################
#   Total forks vs Permutation time
#######################################
p2x = np.array([item['PermutationDepth'] for item in query_speed_vs_permutation_depth_flat])
p2y = np.array([(item['TotalForksMatching']+item['TotalForksPartialMatching'])/item['QueriesProcessed'] for item in query_speed_vs_permutation_depth_flat])
plot_data(p2x,p2y,'Permutation Depth', 'Average number of forks per query')

#######################################
#   Query speed vs permutation depth
#######################################
first_time = query_speed_vs_permutation_depth_flat[0]['AvgQueryTime']
p2x = np.array([item['PermutationDepth'] for item in query_speed_vs_permutation_depth_flat])
p2y = np.array([item['QueriesProcessed']/item['TotalQueryTime'] for item in query_speed_vs_permutation_depth_flat])
plot_data(p2x,p2y,'Permutation Depth', 'Queries Per Second')

#######################################
#   Match speed vs permutation depth
#######################################
first_time = query_speed_vs_permutation_depth_flat[0]['AvgQueryTime']
p2x = np.array([item['PermutationDepth'] for item in query_speed_vs_permutation_depth_flat])
p2y = np.array([item['TotalForksMatching']/item['TotalQueryTime'] for item in query_speed_vs_permutation_depth_flat])
plot_data(p2x,p2y,'Permutation Depth', 'Matches per second')


###############################################################################
#   Query Speed vs Index sample size
###############################################################################
query_speed_vs_sample_rate = parse_directory('../log/query_speed_vs_sample_rate')
query_speed_vs_sample_rate_flat = flatten_and_sort(query_speed_vs_sample_rate, 'SampleRate')
for item in query_speed_vs_sample_rate_flat:
    print(item)

#######################################
#   Query speed vs log sample rate
#######################################
p2x = np.array([np.log(item['SampleRate'])/np.log(2) for item in query_speed_vs_sample_rate_flat])
p2y = np.array([item['QueriesProcessed']/item['TotalQueryTime'] for item in query_speed_vs_sample_rate_flat])
plot_data(p2x,p2y,'Log2 Sampling Rate', 'Queries per second')

#######################################
#   Match speed vs log sample rate
#######################################
p2x = np.array([np.log(item['SampleRate'])/np.log(2) for item in query_speed_vs_sample_rate_flat])
p2y = np.array([item['TotalForksMatching']/item['TotalQueryTime'] for item in query_speed_vs_sample_rate_flat])
plot_data(p2x,p2y,'Log2 Sampling Rate (log ranks per character)', 'Matches per second')

#######################################
#   Log index size vs log sample rate
#######################################
p2x = np.array([np.log(item['SampleRate'])/np.log(2) for item in query_speed_vs_sample_rate_flat])
p2y = np.array([np.log(item['IndexSize'])/np.log(2) for item in query_speed_vs_sample_rate_flat])
plot_data(p2x,p2y,'Log2 Sampling Rate (log ranks per character)', 'Log2 Index Size (log bytes)')