# GIN-TONIC: Non-hierarchical full-text indexing for graphs

**GIN-TONIC** (**G**raph **IN**dexing **T**hrough **O**ptimal **N**ear **I**nterval **C**ompaction), or `gin` for short, is a data structure inspired by the FM-Index for the indexing of directed, string labelled graphs. 
The data structure indexes all possible string walks on the graph and supports various implementations whose memory requirements scale linearly (or sub-linearly) 
scale with the size of the input graph.   

# Paper and Citing
If you'd like to include GIN-TONIC in your work, please cite the paper as follows.

```
@article{10.1093/nargab/lqae159,
    author = {{\"O}zt{\"u}rk, {\"U}nsal and Mattavelli, Marco and Ribeca, Paolo},
    title = {GIN-TONIC: non-hierarchical full-text indexing for graph genomes},
    journal = {NAR Genomics and Bioinformatics},
    volume = {6},
    number = {4},
    pages = {lqae159},
    year = {2024},
    month = {12},
    abstract = {This paper presents a new data structure, GIN-TONIC (Graph INdexing Through Optimal Near Interval Compaction), designed to index arbitrary string-labelled directed graphs representing, for instance, pangenomes or transcriptomes. GIN-TONIC provides several capabilities not offered by other graph-indexing methods based on the FM-Index. It is non-hierarchical, handling a graph as a monolithic object; it indexes at nucleotide resolution all possible walks in the graph without the need to explicitly store them; it supports exact substring queries in polynomial time and space for all possible walk roots in the graph, even if there are exponentially many walks corresponding to such roots. Specific ad-hoc optimizations, such as precomputed caches, allow GIN-TONIC to achieve excellent performance for input graphs of various topologies and sizes. Robust scalability capabilities and a querying performance close to that of a linear FM-Index are demonstrated for two real-world applications on the scale of human pangenomes and transcriptomes. Source code and associated benchmarks are available on GitHub.},
    issn = {2631-9268},
    doi = {10.1093/nargab/lqae159},
    url = {https://doi.org/10.1093/nargab/lqae159},
    eprint = {https://academic.oup.com/nargab/article-pdf/6/4/lqae159/61033863/lqae159.pdf},
}
```
The input data used in this paper is available via another GitHub repository, [gin-data](https://github.com/uensalo/gin-data). The repository contains two graph genomes in `.ging` format.

# Known Issues
- Please compile and use on Unix systems for the time being. This toolbox has not been properly tested with Windows and may provide incorrect results due different handling of line endings.

## Table of Contents

- [Compiling From Source](#compiling-from-source)
- [Reproducing Benchmarks and Data](#reproducing-benchmarks-and-data)
- [Usage](#usage)
  - [General](#general)
  - [Program Specific File Formats](#program-arguments-and-further-usage)
    - [.ging](#ging)
    - [.ginq](#ginq)
    - [.gini](#gini)
    - [.ginc](#ginc)
    - [.ginp](#ginp)
    - [.gine](#gine)
  - [Bare Minimum Indexing Pipeline](#bare-minimum-indexing-scheme)
  - [Full Indexing Scheme](#full-indexing-scheme)
- [Output Format Description of `gin query find` and `gin decode walks`](#output-format-description-of-gin-query-find-and-gin-decode-walks)    
  - [`gin query`](#ginquery)
  - [`gin decode`](#gindecode)
- [Program Arguments and Further Usage](#program-arguments-and-further-usage)
    - [gin:permutation](#ginpermutation)
    - [gin:index](#ginindex)
    - [gin:deindex](#gindeindex)
    - [gin:query](#ginquery)
    - [gin:decode](#gindecode)
    - [gin:utils](#ginutils)
    - [gin:help](#ginhelp)
- [License](#license)


## Compiling From Source
The following script can be used to build the binaries from source:
```bash
mkdir build && (cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make)
```
Submodules might not be initialised when building `gin` for the first time. To make sure, please navigate to the root directory of the project and run
```bash
git submodule update --init
```
If this also fails to clone the repositories, please clone `sdsl-lite` and `libdivsufsort-gin` under the directory `extern` manually and renate `sdsl-lite` to sdsl, i.e.
```bash
cd extern
git clone https://github.com/uensalo/libdivsufsort-gin
git clone https://github.com/simongog/sdsl-lite sdsl
```

The software package currently has six compile time CMake options.

- **`BUILD_OPENMP`** (`set(BUILD_OPENMP ON)`)
  Enables OpenMP support for parallelization; if disabled the -j parameter has no effect.

- **`BUILD_MARCH_NATIVE_FLAG`** (`set(BUILD_MARCH_NATIVE_FLAG ON)`)
  Compiles with the `-march=native` flag for optimizations specific to the host machine's architecture. Note that enabling this might produce binaries incompatible with different architectures.

- **`BUILD_BUILTIN_POPCOUNT`** (`set(BUILD_BUILTIN_POPCOUNT ON)`)
  Available when compiling with `gcc`. Uses the built-in popcount instead of the default, software implementation. May possibly improve indexing times.

- **`BUILD_DEBUG`** (`set(BUILD_DEBUG OFF)`)
  Build in debug mode, with `-O0` and libasan. Only interesting for development purposes.

- **`BUILD_DNA_FMI`** (`set(BUILD_DNA_FMI ON)`)
  Enables DNA-optimised cache friendly FM-Index. Recommended if the input graph contains only (A,C,G,T,N).

- **`BUILD_DNA_FMI_DR`** (`set(BUILD_DNA_FMI_DR ON)`)
  Enables double rank queries with the DNA-optimised FM-Index. `BUILD_DNA_FMI` Must be on.

- **`BUILD_SDSL`** (`set(BUILD_SDSL ON)`)
  Builds the project based on the `sdsl` implementation of an FM-Index. Renders sampling rate parameters useless. Indices generated with the option enabled or disabled have different bitstreams and are not compatible. Enabling the `sdsl` implementation can significantly decrease index size, but will scale poorly if the alphabet is large.

- **`BUILD_ORACLE`** (`set(BUILD_ORACLE ON)`)
  Builds the project with the Oracle Interval-Merge Tree, which further partitions the intervals stored in the tree per alphabet character to decrease the overhead due to having to filter partial matches. Experimental. Disabling this is recommended.

## Reproducing Benchmarks and Data
All relevant scripts to reproduce benchmarks and plots are given under the `benchmark` subdirectory. Firstly, navigate and bootstrap the directories for the benchmarks.
```bash
cd benchmark
cd scripts
./bootstrap.sh
```
This will create the necessary folders under which logs and plots will be placed. Then, under the same directory, run
```bash
./build_bin.sh
```
to build the binaries. Note that the binaries may fail to build for a variety of reasons, which range from not having the submodules initialised or `OpenMP` problems. If you encounter
such problems, please disable `OpenMP` in the `CMakeLists.txt` residing in the root directory.

Next, place the two input items with names `GRCh38-20-0.10b.ging` and `gencode.v40.ging` under `benchmark/input`. `GRCh38-20-0.10b.ging` can be produced via
`rgfa2ging` in the `gin` toolkit, and `gencode.v40.ging` can be produced via `gencode-2-ging` located under `tools` in the root directory of the project.
`gencode-2-ging` [requires BiOCamLib to be installed](https://github.com/PaoloRibeca/BiOCamLib). 

The input files to these programs can be obtained from the [minigraph repository](https://github.com/lh3/minigraph) and [GENCODE Release 40](https://www.gencodegenes.org/human/release_40.html).
For the pangenome, to simplify things, one might also convert nucleotide substitution codes to `N`s.

Then, run
```bash
./benchmark_pangenome_hard.sh
./benchmark_annotation_hard.sh
```
to benchmark `gin`. All related logs will be put under `benchmark/log`. Note that the permutation benchmarks are commented out, as they are memory intensive. Enable them at your own risk.

Finally, running `plot_logs.py` under `benchmark/scripts` will produce the relevant plots under `plots`.

## Usage

### General

The executable `gin` contains a suite of programs that allows string graphs to be indexed and queried. The program has 
six subprograms, namely `index`, `deindex`, `query`, `decode`, `permutation`, `utils`, `help`.
Subprograms may contain multiple modes of execution, which are discussed further into the section.

The indexing and querying algorithm has many moving parts, and setting up all the indices and files requires executing these subprograms in a particular order.
The bare minimum indexing scheme to index and query a graph can be achieved through the programs `index` and `query`, whereas the execution and the incorporation 
of the outputs of other programs might significantly boost the performance of indexing at the cost of extra resources. The bare minimum scheme is not recommended.

### Program Specific File Formats

`gin` works with program specific file formats. 
The file formats are as follows.

#### .ging

Human readable tab separated string graph file supplied as input to various programs under `gin`. The `.ging` format is defined as follows:

Each line defines either a string labelled vertex or a directed edge. Vertex lines start
with a `V` and edge lines start with an `E`. 

`V` lines contain three tab separated fields: first field is always `V`, the second field is a unique vertex ID
and the third field is the string label associated with the vertex. Example:
```
V	6	ATCAGTCATGCAGTACGT...
```

`E` lines contain three tab separated fields: first field is always `E`, the second field is the vertex ID of the source vertex, and the third field is the vertex ID
of the destination vertex. Example:
```
E	3	4
```

There are additional constraints on the format.
  - Vertex IDs appearing in the second column of `V` lines shall be unique.
  - For a file with `N` vertices, the vertex IDs shall be in the range `0,...,N-1`.
  - Source and destination vertex IDs in `E` lines shall be in the range `0,...,N-1`
  - Vertex labels may not contain the tab character or the newline character.
  - All `V` lines appear before `E` lines.

Example conformant `.ging` file:
```
V	0	CAGAC
V	1	TCTCGTA
V	2	AATCA
V	3	CGT
V	4	AAAGCATAT
E	0	3
E	1	0
E	1	2
E	1	3
E	1	4
E	2	4
E	3	0
E	3	4
E	4	0
E	4	3
```
which encodes the following string graph:

<img src="img/ging_example.png" alt="String graph encoded by the example .ging file." width="25%">

#### .ginq

Human readable query file format. One query per line, terminated with an exit prompt (`exit();` by default, can be changed arbitrarily in `gin.c`) to indicate the end of stream.

Example .ginq file:
```
CCGTAAAGC
GT
TC
TCAGA
AAAAGCAT
exit();
```

#### .gini
Binary graph index file.

#### .ginc
Binary precomputed suffix array range cache file.

#### .ginp
Human readable vertex permutation file. One integer per line. `V0` is permuted to the vertex index in the first line, `V1` to the second index, and so on.

Example .ginp file:
```
0
3
2
1
4
```

#### .gine
Binary encoding of an .ging file.

### Bare Minimum Indexing Scheme
To have a minimally working graph index, the following commands suffice:
```bash
./gin index -i <input_ging> -o <output_gini>
./gin query find -r <output_gini> -i <input_ginq> -o <output> --decode
```
The call to `gin index` produces a graph self-index over the input graph, and the call to `gin query` loads the index into memory and runs queries contained in the file `<input_ginq>`, and writes the results to `<output>`. 

### Full Indexing Scheme
To have a fully working graph index with all the components, consider the following commands:
```bash
./gin permutation -i <input_ging> -t <time> -e <init_temperature> -c <cooling> -d <constraint_set_depth> -o <output_perm>
./gin index -i <input_ging> -o <output_gini> -s <sa_sample_rate> -r <rank_sample_rate> -p <input_perm>
./gin query cache -r <input_gini> -c <cache_depth> -o <output_ginc>
./gin query find -r <output_gini> -C <input_ginc> -i <input_ginq> -o <output_roots> --decode
```
The call to `gin permutation` produces a permutation changing the order of vertices. The parameters `-t,-e,-c` dictate the parameters of optimisation of the permutation through simulated annealing. 
The extra call to `gin query cache` produces precomputed suffix array ranges to be used during querying, and will significantly speed up the querying if used. Depths up to `12` are feasible, and a depth between
10 and 12 is recommended.

It is also possible to decode complete walks by using the additional programs provided by the package. There is a perl decoder
under `decoder` with the name `gin_decode_paths.pl` that loads the corresponding `.ging` file in memory and reconstructs walks given the roots.
```bash
./gin_decode_paths.pl -r <input_ging> -i <input_roots> -o <output_walks>
```

Alternatively, `gin` itself also supports an efficient reconstruction scheme where the graph labels are bit-encoded. To decode walks, execute the following:
```bash
./gin decode encode -i <input_ging> -o <output_gine>
./gin decode walks -r <input_gine> -i <input_roots> -o <output_walks>
```
The first command produces a bit-encoded graph that is later used as a graph index in `decode walks`. __Note that the input roots produced by `./gin query find ... --decode`
must be produced without the verbose flag.__

## Output Format Description of `gin query find` and `gin decode walks`

### `gin query find`

The output of `gin query find` is structured in blocks, where each block represents information about a specific string query. The format of each block depends on the enabled flags, in particular `--verbose` and `--decode`.
Each query is responded to by a string query block, which is comprised of a string query line and tabbed lines. The syntax and semantics depend on the flags enabled.

#### 1. String Query Line

Every block starts with a string query (e.g., `CCGTAAAGC`).
Depending on the flags, this line might also contain additional information:
  - If `--verbose` is enabled, the string line contains some metadata, and is of the form `<string>:(c:<count>):(mf:<mfork>,pf:<pfork>):(a:<aval>>)` where
    - `<string>` denotes the query string.
    - `(c:<count>)` denotes the number of occurrences of walk roots.
    - `(mf:<mfork>,pf:<pfork>)` denotes the number of forks with a match (mfork) and forks with partial matches (pfork)
    - `(a:<adv>)` denotes the number of LF-mapping traversals to compute the suffix array ranges associated with the query.
  - Else if `--verbose` is disabled, then the line is of the form `<string>`
The `--decode` flag does not affect this line.

#### 2. Tabbed Lines

Following the string query line, there are zero or more tabbed lines.
Depending on whether `--decode` is enabled these lines represent one of the following.
  - `--decode` disabled: Each line denotes the suffix array match of the query denoted by a half open range in the form `(low,hi)` where `low` is inclusive, `hi` is exclusive.
  - `--decode` enabled: Each line is a vertex-offset pair of the walk roots, given in the format `(v:<vertex>,o:<offset>)`.
  
#### 3. Exit Prompt

The program outputs the exit prompt once it has received the exit prompt to signal other piped `gin` programs to terminate. 

#### 4. Examples

1. **Both `--decode` and `--verbose` enabled**:
```
CCGTAAAGC:(c:1):(mf:1,pf:5):(a:13)
	(v:0,o:4)
GT:(c:2):(mf:1,pf:0):(a:1)
	(v:1,o:4)
	(v:3,o:1)
TC:(c:5):(mf:2,pf:1):(a:3)
	(v:1,o:0)
	(v:1,o:2)
	(v:2,o:2)
	(v:3,o:2)
	(v:4,o:8)
TTTTTTTTT:(c:0):(mf:0,pf:1):(a:1)
	-
TCAGA:(c:2):(mf:1,pf:3):(a:7)
	(v:3,o:2)
	(v:4,o:8)
AAAAGCAT:(c:2):(mf:1,pf:1):(a:8)
	(v:1,o:6)
	(v:2,o:4)
exit();
```

2. **Only `--decode` enabled**:
```
CCGTAAAGC:
	(v:0,o:4)
GT:
	(v:1,o:4)
	(v:3,o:1)
TC:
	(v:1,o:0)
	(v:1,o:2)
	(v:2,o:2)
	(v:3,o:2)
	(v:4,o:8)
TTTTTTTTT:
	-
TCAGA:
	(v:3,o:2)
	(v:4,o:8)
AAAAGCAT:
	(v:1,o:6)
	(v:2,o:4)
exit();
```

3. **Only `--verbose` enabled**:
```
CCGTAAAGC:(c:1):(mf:1,pf:5):(a:13)
	(37,38)
GT:(c:2):(mf:1,pf:0):(a:1)
	(46,48)
TC:(c:5):(mf:2,pf:1):(a:3)
	(52,55)
	(48,50)
TTTTTTTTT:(c:0):(mf:0,pf:1):(a:1)
	-
TCAGA:(c:2):(mf:1,pf:3):(a:7)
	(48,50)
AAAAGCAT:(c:2):(mf:1,pf:1):(a:8)
	(26,28)
exit();
```

4. **Both flags disabled**:
```
CCGTAAAGC:
	(37,38)
GT:
	(46,48)
TC:
	(52,55)
	(48,50)
TTTTTTTTT:
	-
TCAGA:
	(48,50)
AAAAGCAT:
	(26,28)
exit();
```

### `gin decode walks`

Similarly to the output of `gin query find`, the output of `gin decode walks` is organised into blocks of string lines and tabbed lines, but flags don't change the output.

Every block starts with a string query (e.g., `CCGTAAAGC`), followed by zero or more tabbed lines.
Each tabbed line belonging to a string query block has the following syntax.
```
	(o1,oN);v1:...:vN
```
where `o1` denotes the offset of the walk into the label of `v1`, `oN` denotes the offset into the label of `vN`, and `v1,...,vN` denote the vertex IDs of the walk.

**Example output:**
```
CCGTAAAGC:
        (4,5);0:3:4
GT:
        (4,6);1
        (1,1);3
TC:
        (0,2);1
        (2,4);1
        (2,4);2
        (2,1);3:0
        (8,1);4:0
        (8,1);4:3
TTTTTTTTT:
	-
TCAGA:
        (2,4);3:0
        (8,4);4:0
AAAAGCAT:
        (6,7);1:4
        (4,7);2:4
```


## Program Arguments and Further Usage

### gin:permutation

`gin permutation` is designed to approximate a permutation of vertex indices so that they appear consecutively in the Burrows-Wheeler order. This task is computationally challenging as computing such a permutation is NP-Complete. The permutation program addresses this by employing a simulated annealing approach to generate an approximate permutation.

**Parameters:**
- **`--input` or `-i` (Optional):** Path to the input file in rGFA or ging format. If not specified, the program reads from the standard input (stdin).
- **`--rgfa` or `-g` (Optional Flag):** Indicates that the input file is an rGFA file. By default, this is set to false.
- **`--output` or `-o` (Optional):** Path to the output file, which contains one integer per line representing the vertex ID (vid). If not specified, the output is written to the standard output (stdout).
- **`--permutation` or `-p` (Optional):** Path to an initial permutation file to start the optimization process. If not provided, the identity permutation is used as the starting point.
- **`--path-span` or `-s` (Optional Flag):** When set, it disallows constraint sets to span multiple vertices. By default, this is set to false.
- **`--depth` or `-d` (Optional):** Specifies the maximum string length considered in constructing constraint sets. The default value is 4.
- **`--temperature` or `-e` (Optional):** Sets the initial temperature for the simulated annealing process. The default is 1e6.
- **`--cooling` or `-c` (Optional):** Determines the cooling factor for the annealing process. The default value is 0.95.
- **`--time` or `-t` (Optional):** Sets the time limit for the optimization process, in seconds. The default is 15 seconds.
- **`--update` or `-u` (Optional):** Specifies the time interval for informative prints during the optimization process, in seconds. The default interval is 3 seconds.
- **`--threads` or `-j` (Optional):** Number of threads to use for parallel cost computation. The default is 1 thread.
- **`--verbose` or `-v` (Optional):** Provides additional information about the indexing process, including time, progress, and memory requirements.

**Example invocation:**
```bash
./gin permutation -i mygraph.ging -o mygraph-perm.txt -t 300 -u 15 -j 8 -v
```

### gin:index

`gin index` indexes a string-labelled graph and produces a program-specific output for future querying. It accepts input files in rGFA or the program-specific .ging formats. Additionally, a vertex permutation can be provided to increase querying times, with an identity permutation used by default.

**Parameters:**
- **`--input` or `-i` (Optional):** Path to the input file in rGFA or ging format. Default: stdin.
- **`--rgfa` or `-g` (Optional Flag):** Indicates that the input file is an rGFA file. Default: false.
- **`--output` or `-o` (Optional):** Path to the output file, produced in binary gini format. Default: stdout.
- **`--permutation` or `-p` (Optional):** Path to the permutation file. Default: Identity permutation.
- **`--isa-sample-rate` or `-s` (Optional):** Sampling rate of the suffix array. A lower value increases query speeds but results in larger index files. Default: 32.
- **`--rank-sample-rate` or `-r` (Optional):** Frequency of rank caches. A lower value increases query speeds but results in larger index files. Default: 32.
- **`--verbose` or `-v` (Optional Flag):** Provides additional information about the indexing process.

**Example invocation:**
```bash
./gin index -i mygraph.ging -g -o mygraph.gini -p myperm -s 64 -r 64 -v
```

### gin:deindex

`gin deindex` decodes an index file (.gini) back into the input graph and a permutation file. This process can take some time.

**Parameters:**
- **`--input` or `-i` (Optional):** Path to the input file in gini format. Default: stdin.
- **`--output` or `-o` (Optional):** Path to the output file, produced in ging format. Default: stdout.
- **`--permutation` or `-p` (Optional):** Path to the output permutation file. For more help, see `gin permutation`. Default: stdout.
- **`--verbose` or `-v` (Optional Flag):** Provides more information (time) about the decoding process.

**Example invocation:**  
`gin deindex -i myindex.gini -o mygraph.ging -p myperm.ginp -v`


### gin:query

`gin query` loads a graph index in gini format into memory and runs provided queries on the graph. Queries are inputted as either a string per line or in FASTQ format. The program offers two modes: `cache` and `find`.

**Parameters:**
- **`--reference` or `-r` (Required for find, cache):** Path to the index file.
- **`--input` or `-i` (Optional for find):** Path to the input file containing string queries. Default: stdin.
- **`--fastq` or `-f` (Optional Flag for find):** Indicates if queries are in FASTQ format. Default: False.
- **`--output` or `-o` (Optional for find, cache):** Path to the output file. For `find`, outputs (vertex_id, index) or suffix array entries. For `cache`, outputs the cache binary. Default: stdout.
- **`--cache-depth` or `-c` (Optional for cache):** Specifies the depth of the cache. Default: 10.
- **`--cache` or `-C` (Optional for find):** Path to the index cache. Default: None.
- **`--max-forks` or `-m` (Optional for find):** Maximum forks to be tracked per query. `-1` for all forks. Default: -1.
- **`--max-matches` or `-M` (Optional for find):** Maximum number of matches decoded per query. `-1` for all matches. Default: -1.
- **`--decode` or `-d` (Optional Flag for find):** Decodes matches into text space. Default: False.
- **`--batch-size` or `-b` (Optional for find):** Number of queries to process at once. Default: 8.
- **`--threads` or `-j` (Optional for find, cache):** Number of threads for parallel querying. Default: 1.
- **`--verbose` or `-v` (Optional for find, cache):** Provides additional information about the indexing process.

**Example invocation (cache):**
```bash
./gin query cache -r myindex.gini -o myindex_cache.ginc -j 8 -c 10 -v
```

**Example invocation (find):**
```bash
./gin query find -r myindex.gini -i queries.fastq -f -C myindex_cache.ginc -o results.txt -j 8 -m -1 -M 10 -v
```

### gin:decode

`gin decode` is designed to handle bit encoded graphs, with two primary functionalities: encoding and walk enumeration.

**Modes:**

- `encode`: Encodes the input ging graph into a specific format using ceil(log2(s)) bits per character.
- `walks`: Enumerates all walks for a given string from (vertex, offset) pairs. Output format is `<string>:` followed by `(o1,oN);v1:...:vN` entries.

**Parameters:**
- **`--reference` or `-r` (Required for walks):** Path to the index file.
- **`--input` or `-i` (Optional for walks, encode):** Path to the input file containing walk roots or the input ging graph to be encoded. Walk roots produced from query find must be produced without the -v flag. Default: stdin.
- **`--output` or `-o` (Optional for walks, encode):** Path to the output file. Writes the bit encoded graph for `encode`, and the resulting walks for `walks`. Default: stdout.
- **`--batch-size` or `-b` (Optional for walks):** Number of queries to read and process at once. Default: 8.
- **`--threads` or `-j` (Optional for walks, encode):** Number of threads for parallel processing. Default: 1.
- **`--verbose` or `-v` (Optional for walks, encode):** Provides more information about the process.


**Example invocation (encode):**
```bash
./gin decode encode -i mygraph.ging -o mygraph.gine -v
```

**Example invocation (walks):**
```bash
./gin decode walks -r myindex.gine -i queries.query -o results.txt -j 8 -b 16 -v
```

### gin:utils

`gin utils` is a versatile tool for file conversion and analysis in the gin system. It supports four modes:

- **rgfa2ging:** Converts an rGFA file into a ging file. Note that naming conventions and stable sequences are discarded, and all sequences are assumed to be on the forward strand.
- **spectrum:** Extracts the k-mer spectrum from an input ging file.
- **find:** Finds (vertex ID, offset) pairs for query matches, similar to `gin:query`, but uses a brute-force approach. Not recommended.

**Parameters:**

- **`--input` or `-i`:** (Optional for rgfa2ging, fastq2query, find) Path to the input file. For rgfa2ging, input should be in rGFA format; for fastq2query, in FASTQ format; and for find, in query format. Default: stdin.
- **`--output` or `-o`:** (Optional for all modes) Path to the output file. Default: stdout.
- **`--reference` or `-r`:** (Required for find) Path to the input graph file in ging format. Default: stdin.
- **`--threads` or `-j`:** (Optional for find) Number of threads for parallel processing.
- **`--verbose` or `-v`:** (Optional) Provides additional details about the conversion process. Default: false.

**Example invocation (rgfa2ging):**
```bash
./gin utils rgfa2ging -i mygraph.rgfa -o mygraph.ging -v
```
**Example invocation (spectrum):**
```bash
./gin utils spectrum -i mygraph.rgfa -o spectrum.txt -v
```
**Example invocation (find):**
```bash
./gin utils find -i mygraph.gini -o matches.txt -v -j 4
```


### gin:help
`gin:help` provides general information about the functions and usage of other programs within the gin suite. It's designed to assist users in understanding the capabilities and options of various gin tools.
To obtain detailed help for a specific program in the gin toolkit, use the following command format:
```bash
./gin help <program_name>
```
Replace `<program_name>` with the name of the program you want to learn more about. For example, to get help for the `gin:index` program, you would use `gin help index`.

- When `gin help` is called, it simply prints general information about gin programs and then exits.
- Any additional parameters supplied in the command line are ignored.

The gin suite includes the following programs:

- `index`
- `query`
- `decode`
- `permutation`
- `utils`
- `help`

## License

Distributed under GPLv3. See `LICENSE` for more information.
