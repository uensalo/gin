import os
import re

def parse_experiment_logs(experiment_name):
    log_dir = f'../log/{experiment_name}'
    log_files = [f for f in os.listdir(log_dir) if f.endswith('.txt')]

    results = {}

    for log_file in log_files:
        file_path = os.path.join(log_dir, log_file)

        with open(file_path, 'r') as file:
            log_content = file.read()

            results[log_file] = {"permutation": {}, "index": {}, "query": {}}

            # Parse permutation data
            initial_cost = re.search(r"Optimization begins with initial cost ([\d.]+)\.", log_content)
            final_cost = re.search(r"best_cost = ([\d.]+), cur_cost", log_content)
            no_iterations = re.search(r"Iteration (\d+), best_cost", log_content)
            depth = re.search(r"for depth=(\d+)", log_content)
            time = re.search(r"for (\d+) seconds", log_content)

            if initial_cost:
                results[log_file]["permutation"]["initial_cost"] = float(initial_cost.group(1))
            if final_cost:
                results[log_file]["permutation"]["final_cost"] = float(final_cost.group(1))
            if no_iterations:
                results[log_file]["permutation"]["no_iterations"] = int(no_iterations.group(1))
            if depth:
                results[log_file]["permutation"]["depth"] = int(depth.group(1))
            if time:
                results[log_file]["permutation"]["time"] = int(time.group(1))

            # Parse index data
            num_vertices = re.search(r"Number of vertices: (\d+)", log_content)
            num_edges = re.search(r"Number of edges: (\d+)", log_content)
            total_label_len = re.search(r"Total length of vertex labels: (\d+)", log_content)
            index_size = re.search(r"Resulting index size: (\d+) bytes", log_content)
            parsing_time = re.search(r"Parsing:\s+([\d.]+)", log_content)
            indexing_time = re.search(r"Indexing:\s+([\d.]+)", log_content)
            writing_time = re.search(r"Writing:\s+([\d.]+)", log_content)

            if num_vertices:
                results[log_file]["index"]["num_vertices"] = int(num_vertices.group(1))
            if num_edges:
                results[log_file]["index"]["num_edges"] = int(num_edges.group(1))
            if total_label_len:
                results[log_file]["index"]["total_label_len"] = int(total_label_len.group(1))
            if index_size:
                results[log_file]["index"]["index_size"] = int(index_size.group(1))
            if parsing_time:
                results[log_file]["index"]["parsing_time"] = float(parsing_time.group(1))
            if indexing_time:
                results[log_file]["index"]["indexing_time"] = float(indexing_time.group(1))
            if writing_time:
                results[log_file]["index"]["writing_time"] = float(writing_time.group(1))

            # Parse query data
            total_query_time = re.search(r"Total querying time in seconds: ([\d.]+)", log_content)
            num_queries_processed = re.search(r"Queries processed: (\d+)", log_content)
            avg_query_time = re.search(r"Average time per query: ([\d.]+)", log_content)
            total_forks_matching = re.search(r"Total forks for matching queries: (\d+)", log_content)
            total_forks_partial_matching = re.search(r"Total forks for partially matching queries: (\d+)", log_content)
            avg_forks_matching = re.search(r"Average forks per matching query: ([\d.]+)", log_content)
            avg_forks_partial_matching = re.search(r"Average forks per partially matching query: ([\d.]+)", log_content)

            if total_query_time:
                results[log_file]["query"]["total_query_time"] = float(total_query_time.group(1))
            if num_queries_processed:
                results[log_file]["query"]["num_queries_processed"] = int(num_queries_processed.group(1))
            if avg_query_time:
                results[log_file]["query"]["avg_query_time"] = float(avg_query_time.group(1))
            if total_forks_matching:
                results[log_file]["query"]["total_forks_matching"] = int(total_forks_matching.group(1))
            if total_forks_partial_matching:
                results[log_file]["query"]["total_forks_partial_matching"] = int(total_forks_partial_matching.group(1))
            if avg_forks_matching:
                results[log_file]["query"]["avg_forks_matching"] = float(avg_forks_matching.group(1))
            if avg_forks_partial_matching:
                results[log_file]["query"]["avg_forks_partial_matching"] = float(avg_forks_partial_matching.group(1))

    return results

if __name__ == "__main__":
    import sys
    experiment_name = sys.argv[1]
    results = parse_experiment_logs(experiment_name)
    print(results)
