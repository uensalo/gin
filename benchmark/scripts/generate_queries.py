import sys
import random
import argparse
from collections import defaultdict

def sample_strings(input_file, length, num_samples, output_file, seed=None):
    if seed is not None:
        random.seed(seed)

    vertices = {}
    adjacency_list = defaultdict(list)
    total_label_length = 0
    reciprocals_sum = 0

    with open(input_file, 'r') as f:
        for line in f:
            type, id, data = line.strip().split('\t')
            id = int(id)
            if type == 'V':
                vertices[id] = data
                total_label_length += len(data)
                reciprocals_sum += 1 / len(data)
            elif type == 'E':
                _, from_vertex, to_vertex = line.strip().split('\t')
                adjacency_list[int(from_vertex)].append(int(to_vertex))

    probabilities = {k: len(v) / total_label_length for k, v in vertices.items()}
    inverse_probabilities = {k: (1 / len(v)) / reciprocals_sum for k, v in vertices.items()}

    with open(output_file, 'w') as f:
        for i in range(num_samples):
            if i < num_samples / 2:
                vertex = random.choices(list(vertices.keys()), weights=probabilities.values(), k=1)[0]
            else:
                vertex = random.choices(list(vertices.keys()), weights=inverse_probabilities.values(), k=1)[0]

            while True:
                path = [vertex]
                position = random.randint(0, len(vertices[vertex]) - 1)
                first_offset = position
                string = vertices[vertex][position:position + length]
                needed = length - len(string)

                while len(string) < length:
                    if not adjacency_list[vertex]:
                        break
                    vertex = random.choice(adjacency_list[vertex])
                    path.append(vertex)
                    needed = length - len(string)
                    string += vertices[vertex][:needed]

                last_offset = needed - 1 if needed > 0 else len(vertices[path[-1]]) - 1

                # Check if the string contains the letter 'N'
                if 'N' in string:
                    continue

                if len(string) == length:
                    f.write("{} ,{}\t{}\n".format('\t'.join(map(str, path)), first_offset, last_offset))
                    print(string)
                    break

        print("exit();")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str)
    parser.add_argument('length', type=int)
    parser.add_argument('num_samples', type=int)
    parser.add_argument('output_file', type=str)
    parser.add_argument('seed', type=int, nargs='?', default=None)
    args = parser.parse_args()
    sample_strings(args.input_file, args.length, args.num_samples, args.output_file, args.seed)

if __name__ == '__main__':
    main()
