import argparse

def create_fmdg_file(mode, V, output_path):
    if mode not in ["connected", "bipartite"]:
        raise ValueError("Mode must be either 'connected' or 'bipartite'.")

    with open(output_path, "w") as file:
        if mode == "connected":
            for i in range(V):
                file.write(f"V\t{i}\t{'A' * (i+1)}\n")
            for i in range(V):
                for j in range(i+1, V):
                    file.write(f"E\t{i}\t{j}\n")
                    file.write(f"E\t{j}\t{i}\n")
        elif mode == "bipartite":
            if V % 2 != 0:
                raise ValueError("For bipartite mode, V must be an even number.")
            for i in range(V):
                file.write(f"V\t{i}\t{'A' * ((i % (V//2)) + 1)}\n")
            for i in range(V // 2):
                for j in range(V // 2, V):
                    file.write(f"E\t{i}\t{j}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create an FMDG file.")
    parser.add_argument('mode', type=str, help="Mode: 'connected' or 'bipartite'.")
    parser.add_argument('V', type=int, help="Number of vertices.")
    parser.add_argument('output_path', type=str, help="Output file path.")
    args = parser.parse_args()

    create_fmdg_file(args.mode, args.V, args.output_path)
