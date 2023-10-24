# GIN-TONIC: One-shot FM-Indexing for Arbitrary String Labelled Graphs

**GIN-TONIC** (**G**raph **I**ndexing **T**hrough **O**ptimal **N**ear **I**nterval **C**ompaction), or `gin` for short, is an implementation of FM-GIN, a data structure inspired by the FM-Index for the indexing of directed, string labelled graphs of (almost) arbitrary topology. 
The data structure indexes all possible string walks on the graph and supports various implementations whose memory requirements scale linearly (or sub-linearly) 
scale with the size of the input graph.   

## Table of Contents

- [Compiling From Source](#compiling-from-source)
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
  - [Full Indexing Pipeline](#full-indexing-pipeline)
- [Output Format Description of `gin query find` and `gin decode walks`](#output-format-description-of-gin-query-find-and-gin-decode-walks)    
  - [`gin query`](#ginquery)
  - [`gin decode`](#gindecode)
- [Program Arguments and Further Usage](#program-arguments-and-further-usage)
    - [gin:permutation](#ginpermutation)
    - [gin:index](#ginindex)
    - [gin:query](#ginquery)
    - [gin:decode](#gindecode)
    - [gin:utils](#ginutils)
    - [gin:help](#ginhelp)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Compiling From Source
The following script can be used to build the binaries from source:
```
mkdir build && (cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make)
```
The software package currently has four compile time CMake options.

- **`BUILD_OPENMP`** (`set(BUILD_OPENMP ON)`)
  Enables OpenMP support for parallelization; if disabled the -j parameter has no effect.

- **`BUILD_MARCH_NATIVE_FLAG`** (`set(BUILD_MARCH_NATIVE_FLAG ON)`)
  Compiles with the `-march=native` flag for optimizations specific to the host machine's architecture. Note that enabling this might produce binaries incompatible with different architectures.

- **`BUILD_TESTS`** (`set(BUILD_TESTS OFF)`)
  Disables building of the project's fuzz tests. Only interesting for development purposes.

- **`BUILD_DEBUG`** (`set(BUILD_DEBUG OFF)`)
  Build in debug mode, with `-O0` and libasan. Only interesting for development purposes.


## Usage

### General

The executable `gin` contains a suite of programs that allows string graphs to be indexed and queried. The program has 
six subprograms, namely `index`, `query`, `decode`, `permutation`, `utils`, `help`.
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
```
./gin index -i <input_ging> -o <output_gini>
./gin query -r <output_gini> -i <input_ginq> -o <output> --decode
```
The call to `gin index` produces a graph self-index over the input graph, and the call to `gin query` loads the index into memory and runs queries contained in the file `<input_ginq>`, and writes the results to `<output>`. 

### Full Indexing Pipeline

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

gin permutation approximates a permutation of vertex indices.

### gin:index

This program indexes a string labelled graph and produces a program-specific output for future querying.

**Parameters:**

- `--input` or `-i`: _Description here._
- `--rgfa` or `-g`: _Description here._
  ... _(continue for other parameters)_

### gin:query

gin query loads a graph index in gini format into memory and runs the provided queries on the graph.

**Parameters:**

- `--reference` or `-r`: _Description here._
  ... _(continue for other parameters)_

### gin:decode

gin decode loads a bit encoded graph into memory and enumerates full walks from the output of gin query find.

**Parameters:**

- `--reference` or `-r`: _Description here._
  ... _(continue for other parameters)_

**Parameters:**

- `--input` or `-i`: _Description here._
  ... _(continue for other parameters)_

### gin:utils

gin utils converts rgfa and FASTQ files to ging files and gin query files.

**Parameters:**

- `--input` or `-i`: _Description here._
  ... _(continue for other parameters)_

### gin:help

gin help prints general information about how other programs under gin work.

## Examples

_TODO: Provide some example commands and their expected outputs._

## Contributing

_TODO: Guidelines for those who want to contribute._

## License

_TODO: License details._
