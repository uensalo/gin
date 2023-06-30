import random
import sys

def generate_random_graph(num_vertices, num_edges, label_length, seed=None):
    # Seed the random number generator if a seed is provided
    if seed is not None:
        random.seed(seed)

    # Ensure it's possible to create the requested number of edges
    if num_edges > num_vertices * (num_vertices - 1) // 2:
        raise ValueError(f"Cannot create {num_edges} edges with {num_vertices} vertices")

    # Define Cantor pairing function
    def cantor_pairing(a, b):
        return (a+b)*(a+b+1)//2 + b if a != b else -1

    # Keep track of existing edges
    existing_edges = set()

    # Generate vertices with labels restricted to A,C,G,T
    vertices = [(i, ''.join(random.choices('ACGT', k=label_length))) for i in range(num_vertices)]

    # Generate edges
    edges = []
    while len(edges) < num_edges:
        a, b = random.sample(range(num_vertices), 2)
        pair = cantor_pairing(a, b)
        if pair not in existing_edges:
            edges.append((a, b))
            existing_edges.add(pair)

    # Print to stdout
    for vertex in vertices:
        print(f"V\t{vertex[0]}\t{vertex[1]}")
    for edge in edges:
        print(f"E\t{edge[0]}\t{edge[1]}")

if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python script.py num_vertices num_edges label_length [seed]")
        sys.exit(1)

    num_vertices = int(sys.argv[1])
    num_edges = int(sys.argv[2])
    label_length = int(sys.argv[3])
    seed = int(sys.argv[4]) if len(sys.argv) == 5 else None

    generate_random_graph(num_vertices, num_edges, label_length, seed)
