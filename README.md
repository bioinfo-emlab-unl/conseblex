# ConSemblEX
An extension of ConSemble, a consensus-based ensemble approach to improve transcriptome assembly.
ConSemblex consists of three modules, Assembly, Read-processing, and Analysis;

- **Assembly**: This module is used to run individual assembly pipelines from given single-end or paired-end read files for specified assembly methods.  

- **Analysis**: This module is used to perform the final consensus assembly between given assembly files from chosen assemblies methods.

- **Read-processing**: This module compliments the **Assembly** module. It provides means to filter and normalize reads before assembly.   

More details about each module and examples runs are provided in the later sections of the document.

**README update**: last updated on September 18, 2022.

**Tested on**: Linux

### Assembly
### Dependancies
Ensure that all the individual assembly methods are available on your environment before running the Assembly pipeline.
1. [Python3 ^v3.9.x](https://www.python.org/downloads/)
2. [Biopython v1.78](https://biopython.org/)
3. [khmer v2.0](https://github.com/dib-lab/khmer)
4. [erne-filter v2.1](http://erne.sourceforge.net/index.php)

Only the exact versions of the above libraries have been tested. However, later versions are expected to work.

### Usage

#### *De novo* assembly input and options
```bash
usage: consemblex.py denovo [-h] [-a ASSEMBLYNAME] [-o OUTPUT] [-n | -c NORMCOV] [-f] [-m ASSEMBLERS] [-t] [-q Q_VALUE] [-s SINGLE | -p reads1 reads2]

optional arguments:
  -h, --help            show this help message and exit

output:
  -a ASSEMBLYNAME, --assemblyName ASSEMBLYNAME
                        specifies the assembly name (e.g testAssembly) to organize output files
  -o OUTPUT, --output OUTPUT
                        changes the base directory for output files (default is current working directory)

normalization:
  -n, --noNorm          skip digital normalization of filtered reads
  -c NORMCOV, --normCov NORMCOV
                        change the max kmer coverage for digital normalization - lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)

configurations:
  -f, --force           ignores the log file and runs entire pipeline
  -m ASSEMBLERS, --assemblers ASSEMBLERS
                        file containing assemblers and their respective configurations
  -t, --time            records the amount of time each stage of the pipeline takes

quality filter:
  -q Q_VALUE, --q_value Q_VALUE
                        changes q-value threshold for reads filtering by erne-filter (default 20)

reads:
  -s SINGLE, --single SINGLE
                        specifies the file including single-end reads
  -p reads1 reads2, --paired reads1 reads2
                        specifies the paired-end reads files
```

Help:
```bash
python consemblex.py d -h
```

OR

```bash
python consemblex.py denovo --help
```

*De novo* assembly:
```bash
python consemblex.py d -a col0 -o output -p data/col0/merged_norm_1.fq data/col0/merged_norm_2.fq
```
#### Guided assembly input and options
```bash
usage: consemblex.py guided [-h] [-a ASSEMBLYNAME] [-o OUTPUT] [-n | -c NORMCOV] [-f FORCE] [-m ASSEMBLERS] [-t] [-q Q_VALUE] [-r REFERENCE] [-s SINGLE | -p reads1 reads2]

optional arguments:
  -h, --help            show this help message and exit

output:
  -a ASSEMBLYNAME, --assemblyName ASSEMBLYNAME
                        specifies the assembly name (e.g testAssembly) to organize output files
  -o OUTPUT, --output OUTPUT
                        changes the base directory for output files (default is current working directory)

normalization:
  -n, --noNorm          skip digital normalization of filtered reads
  -c NORMCOV, --normCov NORMCOV
                        change the max kmer coverage for digital normalization - lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)

configurations:
  -f FORCE, --force FORCE
                        ignores the log file and runs entire pipeline
  -m ASSEMBLERS, --assemblers ASSEMBLERS
                        file containing assemblers and their respective configurations
  -t, --time            records the amount of time each stage of the pipeline takes

quality filter:
  -q Q_VALUE, --q_value Q_VALUE
                        changes q-value threshold for reads filtering by erne-filter (default 20)

transcriptome reference:
  -r REFERENCE, --reference REFERENCE
                        location for the reference Fasta file for read mapping

reads:
  -s SINGLE, --single SINGLE
                        specifies the file including single-end reads
  -p reads1 reads2, --paired reads1 reads2
                        specifies the paired-end reads files
```

Help:
```bash
python consemblex.py g -h
```

OR

```bash
python consemblex.py guided --help
```

Guided assembly:
```bash
python consemblex.py d -a col0 -o output -r Col0ref.faa -p data/col0/merged_norm_1.fq data/col0/merged_norm_2.fq
```

### Output files
/*TODO*/

## Analysis
### Dependancies
Ensure that the following libraries/programs exist on your environment before running the Analysis pipeline.
1. [Python3 ^v3.9.x](https://www.python.org/downloads/)
2. [Biopython v1.78](https://biopython.org/)
3. [ORFfinder linux-i64](https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/)

Only the exact versions of the above libraries have been tested. However, later versions are expected to work.

### Usage
#### Input and options
```bash
usage: consemblex.py analyze [-h] [-a ASSEMBLIES [ASSEMBLIES ...]] [-n NAME] [-w] [-b BENCHMARK] [-p]

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLIES [ASSEMBLIES ...], --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        specifies a list of the assemblies and their respective files i.e SPades spades.fa Trinity trinity.fa ...
  -n NAME, --name NAME  specifies the analysis name (e.g testAnalysis) to organize output files
  -w, --write           write intermediate output to files
  -b BENCHMARK, --benchmark BENCHMARK
                        bechmark transcriptome
  -p, --aminoacids      specifies whether the files are amino acid sequences files (default is nucleotide sequences fies)
```

Example usage:

Help:
```bash
python consemblex.py x -h
```

OR

```bash
python consemblex.py analyze --help
```

Analysis:
```bash
python consemblex.py analyze -n col0 -wp -b Col0ref.faa -a soap SOAP.fasta spades rnaSPAdes.fasta idba idba.fasta Trinity.fasta
```

### Output files
ConSemblEX creates the *output* folder in the current working directory. All the output files from the analysis are saved to the *analysis* folder within the *output* folder. Inside the *analysis* folder, a folder with the name provided at execution time is created in which all the final files are saved. Following the run example provided above, the output files would appear as shown below;
```bash
col0
├── 1
│   ├── consemblex1+.faa
│   ├── consemblex1+_hits.txt
│   ├── consemblex1+_longest.fasta
│   ├── consemblex1+_ref.txt
│   ├── consemblex1+_shortest.fasta
│   ├── consemblex1.faa
│   ├── consemblex1.fasta
│   ├── consemblex1_hits.txt
│   └── consemblex1_ref.txt
├── 2
│   ├── consemblex2+.faa
│   ├── consemblex2+_hits.txt
│   ├── consemblex2+_longest.fasta
│   ├── consemblex2+_ref.txt
│   ├── consemblex2+_shortest.fasta
│   ├── consemblex2.faa
│   ├── consemblex2.fasta
│   ├── consemblex2_hits.txt
│   └── consemblex2_ref.txt
├── 3
│   ├── consemblex3+.faa
│   ├── consemblex3+_hits.txt
│   ├── consemblex3+_longest.fasta
│   ├── consemblex3+_ref.txt
│   ├── consemblex3+_shortest.fasta
│   ├── consemblex3.faa
│   ├── consemblex3.fasta
│   ├── consemblex3_hits.txt
│   └── consemblex3_ref.txt
├── 4
│   ├── consemblex4+.faa
│   ├── consemblex4+_hits.txt
│   ├── consemblex4+_longest.fasta
│   ├── consemblex4+_ref.txt
│   ├── consemblex4+_shortest.fasta
│   ├── consemblex4.faa
│   ├── consemblex4.fasta
│   ├── consemblex4_hits.txt
│   └── consemblex4_ref.txt
├── 5
│   ├── consemblex5.faa
│   ├── consemblex5.fasta
│   ├── consemblex5_hits.txt
│   └── consemblex5_ref.txt
├── longestMergedORFs.faa
├── mergedORFs.faa
├── mergedTranscripts.faa
└── statistics.csv

```

- *1* contains all files for the one-way or greater consensus, *2* contains all files for the two-way or greater consensus and so on...
- *longestMergedORFs.fa*: longest ORFs among each set of ORFs for a specific nucleotide sequence.
- *mergedTranscripts.faa*: all nucleotide transcripts from the provided assembly methods.
- *mergedORFs.faa*: all ORFs resulting from running ORFfinder on *mergedTranscripts.faa*.
- In each numbered folder, the contained files are;
  - *consemblex{num}.faa*: consensus amino-acids represented by *contig-Ids* and actual sequences.
  - *consemblex{num}_hits.txt*: table of assembled contigs and their conrresponding benchmark transcript Ids. *Not hit* indicates that no matching benchmark transcript was found.
  - *consemblex{num}_longest.fasta*: representative longest nucleotides sequences among the assembled contigs.
  - *consemblex{num}_ref.txt*: table reference of assembled contigs, their corresponding assembly methods and the original nucleotide sequence Ids.  
  - *consemblex{num}_shortest.fasta*: representative lonshortestgest nucleotides sequences among the assembled contigs. 

### Read Processing
### Dependancies
Ensure that all the individual assembly methods are available on your environment before running the Read processing pipeline.
1. [Python3 ^v3.9.x](https://www.python.org/downloads/)
2. [Biopython v1.78](https://biopython.org/)
3. [khmer v2.0](https://github.com/dib-lab/khmer)
4. [erne-filter v2.1](http://erne.sourceforge.net/index.php)

Only the exact versions of the above libraries have been tested. However, later versions are expected to work.

### Usage
#### Input and options
```bash
usage: consemblex.py process_reads [-h] [-s SINGLE | -p reads1 reads2] [-o OUTPUT] [-n NAME] [-q Q_VALUE] [-c NORMCOV]

optional arguments:
  -h, --help            show this help message and exit

reads:
  -s SINGLE, --single SINGLE
                        specifies the file containing single-end reads
  -p reads1 reads2, --paired reads1 reads2
                        specifies the paired-end reads files

output:
  -o OUTPUT, --output OUTPUT
                        changes the base directory for output files (default is current working directory)
  -n NAME, --name NAME  specifies the process reads name (e.g col0) to organize output files

configurations:
  -q Q_VALUE, --q_value Q_VALUE
                        changes q-value threshold for reads filtering by erne-filter (default 20)
  -c NORMCOV, --normCov NORMCOV
                        change the max kmer coverage for digital normalization; lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)
```

Help:
```bash
python consemblex.py p -h
```

OR

```bash
python consemblex.py process_reads --help
```

Process reads:
```bash
python consemblex.py pr -n col0 -o process_reads -p data/col0-Alt-reads_1.fq data/col0-Alt-reads_2.fq
```

### Output files
/*TODO*/