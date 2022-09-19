__version__ = 2.0

import argparse
import os, re, time
import warnings
import logging
import importlib.util as implib

from Bio import SeqIO
from sys import maxsize
from lib.analysis import Analysis
from lib.util import Util as utillib
from lib.assembly import ProcessReads

# Package names
biopython_package = "Bio"
# biopython_package = "biopython"

# Global variables
basedir = os.getcwd()
logFile = "consemblex.log"
output = os.path.join(basedir, 'output/')
utillib.create_directory(directorypath=output)

## Process reads function
def run_process_reads(args):
    # Create process reads object
    process_reads = ProcessReads(args=args, basedir=output)

    # Create output folder
    utillib.create_directory(directorypath=process_reads.output)

    # Configure logging
    utillib.configure_logging(logfile=process_reads.logfile, level=logging.DEBUG)

    logging.info("Started processing of reads")
    tic = time.perf_counter()
    print("================================== Processing reads ==============================================")

    ## Normalize reads
    process_reads.normalize_reads()

    toc = time.perf_counter()
    logging.info(f"Finished processing reads in {toc-tic:0.4f} seconds")
    print("================================== Finished processing reads ==================================")

##  Denovo assembly function
def run_denovo_assembly(args):
    assemblyname = args.assemblyName
    outputdir = os.path.join(args.output, f"assembly/{assemblyname}")
    reads1 = args.paired[0]
    reads2 = args.paired[1]
    # Location of BayesDemovo script/executable if not in path
    bayesdenovo = "assemblers/BayesDenovo-v1/BayesDenovo_v1.pl"

    # Create output folder
    utillib.create_directory(directorypath=outputdir)
    
    # Configure logging
    assemblylog = os.path.join(outputdir, "assembly.log")
    utillib.configure_logging(logfile=assemblylog, level=logging.DEBUG)

    logging.info("Started running assembly")
    print("================================== Running assembly ==============================================")
    
    tic = time.perf_counter()

    # BayesDenovo assembly
    kmers = list(range(25, 45, 4)) # List of kmers

    for kmer in kmers:
        subdir = os.path.join(outputdir, f"kmer_{kmer}")
        # utillib.create_directory(directorypath=subdir)
        # Command to run
        msg = f"An error occured while running BayesDenovo assembly for kmer of size {kmer}"
        command = f"{bayesdenovo} -k {kmer} --seqType fq --left {reads1} --right {reads2} --CPU 8 --output {subdir} \
                    --bamfile {subdir}/bayesdenovo.sam --fafile {subdir}/bayesdenovo.fa"

        proc = utillib.execute_process(command=command, errormessage=msg)
    
    toc = time.perf_counter()

    logging.info(f"Finished running assembly in {toc-tic:0.4f} seconds")
    print("================================== Finished running assembly ==================================")

##  Guided assembly function
def run_guided_assembly(args):
    pass

## Analysis function
def run_analysis(args):
    ## Create analysis object
    analysis = Analysis(args=args, basedir=output)

    ## Create output folder
    utillib.create_directory(directorypath=analysis.output)

    ## Configure logging
    utillib.configure_logging(logfile=analysis.logfile, level=logging.DEBUG)

    logging.info("Started running analysis of assemblies")
    print("================================== Running analysis ==============================================")

    ## Merge assembly files 
    logging.info("Started merging assembly files")
    tic = time.perf_counter()
    analysis.merge_files()
    toc = time.perf_counter()
    print(f"Merged files in {toc - tic:0.4f} seconds")
    logging.info(f"Finished merging assembly files in {toc - tic:0.4f} seconds")
    logging.info("Assembly files merged to {}".format(analysis.mergedtranscripts_file))
                    
    ## Run ORFfinder
    logging.info("Started running of ORFfinder to find ORFs on merged file")
    tic = time.perf_counter()
    mergedORFs = os.path.join(analysis.output, 'mergedORFs.faa')
    orfloc = "singularity exec docker://unlhcc/orffinder ORFfinder"
    command = "{} -in {} -out {}".format(orfloc, analysis.mergedtranscripts_file, mergedORFs)
    analysis.orf_finder(command=command, outfile=mergedORFs)
    toc = time.perf_counter()
    print(f"Ran ORF finder in {toc - tic:0.4f} seconds")
    logging.info(f"Finished running of ORFfinder to find ORFs on merged file in {toc - tic:0.4f} seconds")
    logging.info("ORFs written to {}".format(mergedORFs))

    # Get longest ORFs from merged ORF file
    logging.info("Started finding longest ORFs in merged ORF file")
    tic = time.perf_counter()
    analysis.get_longest_orfs()
    toc = time.perf_counter()
    print(f"Ran longest ORF finder in {toc - tic:0.4f} seconds")
    logging.info(f"Finished finding longest ORFs in merged ORF file in {toc - tic:0.4f} seconds")
    logging.info("Longest ORFs written to {}".format(analysis.longestmergedORFs_file))

    ## cluster contigs that have identical proteins
    logging.info("Started clustering of contigs relative to longest ORFs found")
    tic = time.perf_counter()
    analysis.get_sequence_clusters()
    toc = time.perf_counter()
    print(f"Ran clustering of contigs in {toc - tic:0.4f} seconds")
    logging.info(f"Finished clustering of contigs relative to longest ORFs found in {toc - tic:0.4f} seconds")

    ## Get consensus contigs
    logging.info("Started finding consensus contigs")
    
    tic = time.perf_counter()

    nucleotidelongest, nucleotideshortest, aminocontigs, contigsref = analysis.get_consensus_sequences()
    if args.benchmark:
        logging.info("Started benchmarking of contigs")
        tic = time.perf_counter()
        asemblyhits = analysis.run_benchmark(contigs=aminocontigs)

        # Calculate statistics
        analysis.calculate_statistics()

        toc = time.perf_counter()
        print(f"Ran benchmark of contigs in {toc - tic:0.4f} seconds")
        logging.info(f"Finished benchmark of contigs relative to longest ORFs found in {toc - tic:0.4f} seconds")

    # Write to output files
    overlaps = len(analysis.assemblynames)
    for i in range(1,overlaps+1):
        dest = os.path.join(analysis.output, "{}".format(i))
        utillib.create_directory(directorypath=dest)

        # Longest sequences
        consemble_longest = os.path.join(analysis.output, "{}/consemblex{}+_longest.fasta".format(i, i))
        utillib.write_id_and_sequence(data=nucleotidelongest, file=consemble_longest, final=True, overlap=i, plusonly=False,
                                        overlapmax=overlaps)

        # Shortest sequences
        consemble_shortest = os.path.join(analysis.output, "{}/consemblex{}+_shortest.fasta".format(i, i))
        utillib.write_id_and_sequence(data=nucleotideshortest, file=consemble_shortest, final=True, overlap=i, plusonly=False,
                                        overlapmax=overlaps)

        # Amino sequences 
        consemble_fa = os.path.join(analysis.output, "{}/consemblex{}+.faa".format(i, i))
        utillib.write_id_and_sequence(data=aminocontigs, file=consemble_fa, final=True, overlap=i, plusonly=False, 
                                        overlapmax=overlaps)

        # Benchmark results
        if args.benchmark:
            consemble_fa_hits = os.path.join(analysis.output, "{}/consemblex{}+_hits.txt".format(i, i))
            utillib.write_benchmark_ref(data=asemblyhits, file=consemble_fa_hits, overlap=i, plusonly=False, overlapmax=overlaps)

        # Reference
        consemble_ref = os.path.join(analysis.output, "{}/consemblex{}+_ref.txt".format(i, i))
        utillib.write_ref(data=contigsref, names=analysis.contigsassemblers, file=consemble_ref, overlap=i, plusonly=False, 
                            overlapmax=overlaps)
    
    toc = time.perf_counter()
    print(f"Ran get consensus contigs in {toc - tic:0.4f} seconds")

    logging.info(f"Finished finding consensus contigs in {toc - tic:0.4f} seconds")
    
        # Close all open output files
        # list(map(lambda f: f.close(), openFiles))

    print("================================== Finished running analysis ==================================")
    logging.info("Finished running analysis of assemblies")

## Check if package is installed
def package_installed(packageName):
    spec = implib.find_spec(packageName)
    return spec

## Process command line arguments
parser = argparse.ArgumentParser(description="ConSemblEX: A consensus-based ensemble approach to improve transcriptome assembly")
parser.add_argument("-v", "--version", action="version", version="ConSemblEX pipeline version {}".format(__version__), 
                    help="print pipeline version")
subparsers = parser.add_subparsers(title="available commands", help="specify whether to run assembly or analysis")

## Add denovo command
denovo = subparsers.add_parser("denovo", aliases=['dn', 'd'], help="run denovo assembly methods")
denovo_output = denovo.add_argument_group("output")
denovo_output.add_argument("-a", "--assemblyName", help="specifies the assembly name (e.g testAssembly) to organize output files")
denovo_output.add_argument("-o", "--output", help="changes the base directory for output files (default is current working directory)")
# Mutually exclusive normalization
denovo_norm = denovo.add_argument_group("normalization")
denovo_norm_mutual = denovo_norm.add_mutually_exclusive_group()
denovo_norm_mutual.add_argument("-n","--noNorm", action="store_true", help="skip digital normalization of filtered reads")
denovo_norm_mutual.add_argument("-c", "--normCov", type=int, default=50, help="change the max kmer coverage for digital normalization - lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)")
denovo_config = denovo.add_argument_group("configurations")
denovo_config.add_argument("-f", "--force", action="store_true", help="ignores the log file and runs entire pipeline")
denovo_config.add_argument("-m", "--assemblers", help="file containing assemblers and their respective configurations")
denovo_config.add_argument("-t", "--time", action="store_true", help="records the amount of time each stage of the pipeline takes")
denovo_filter = denovo.add_argument_group("quality filter")
denovo_filter.add_argument("-q", "--q_value", type=int, default=20, help="changes q-value threshold for reads filtering by erne-filter (default 20)")
# Mutually exclusive single reads and paired-end reds
denovo_reads = denovo.add_argument_group("reads")
denovo_reads_mutual = denovo_reads.add_mutually_exclusive_group()
denovo_reads_mutual.add_argument("-s", "--single", help="specifies the file including single-end reads")
denovo_reads_mutual.add_argument("-p", "--paired", nargs=2, metavar=("reads1", "reads2"), help="specifies the paired-end reads files")
# Specify function to run
denovo.set_defaults(func=run_denovo_assembly)

## Add guided command
guided = subparsers.add_parser("guided", aliases=['gd', 'g'], help="run guided assembly methods")
guided_output = guided.add_argument_group("output")
guided_output.add_argument("-a", "--assemblyName", help="specifies the assembly name (e.g testAssembly) to organize output files")
guided_output.add_argument("-o", "--output", help="changes the base directory for output files (default is current working directory)")
# Mutually exclusive normalization
guided_norm = guided.add_argument_group("normalization")
guided_norm_group = guided_norm.add_mutually_exclusive_group()
guided_norm_group.add_argument("-n","--noNorm", action="store_true", help="skip digital normalization of filtered reads")
guided_norm_group.add_argument("-c", "--normCov", type=int, default=50, help="change the max kmer coverage for digital normalization - \
                            lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)")
guided_config = guided.add_argument_group("configurations")
guided_config.add_argument("-f", "--force", help="ignores the log file and runs entire pipeline")
guided_config.add_argument("-m", "--assemblers", help="file containing assemblers and their respective configurations")
guided_config.add_argument("-t", "--time", action="store_true", help="records the amount of time each stage of the pipeline takes")
guided_filter = guided.add_argument_group("quality filter")
guided_filter.add_argument("-q", "--q_value", type=int, default=20, help="changes q-value threshold for reads filtering by erne-filter \
                            (default 20)")
guided_reference = guided.add_argument_group("transcriptome reference")
guided_reference.add_argument("-r", "--reference", help="location for the reference Fasta file for read mapping")
# Mutually exclusive single reads and paired-end reds
guided_reads = guided.add_argument_group("reads")
guided_reads_mutual = guided_reads.add_mutually_exclusive_group()
guided_reads_mutual.add_argument("-s", "--single", help="specifies the file including single-end reads")
guided_reads_mutual.add_argument("-p", "--paired", nargs=2, metavar=("reads1", "reads2"), help="specifies the paired-end reads files")
# Specify function to run
guided.set_defaults(func=run_guided_assembly)

# Add analysis command
# This does not run the assembly but performs analysis on already assembled files
analysis = subparsers.add_parser("analyze", aliases=['an', 'x'], help="run analysis on already assembled transcript files")
analysis.add_argument("-a", "--assemblies", nargs='+', help="specifies a list of the assemblies and their respective files i.e \
                    SPades spades.fa Trinity trinity.fa ...")
analysis.add_argument("-n", "--name", default="result", help="specifies the analysis name (e.g testAnalysis) to organize output files")
analysis.add_argument("-w", "--write", action="store_true", help="write intermediate output to files")
analysis.add_argument("-b", "--benchmark", help="bechmark transcriptome")
analysis.add_argument("-p", "--aminoacids", action="store_true", help="specifies whether the files are amino acid sequences files \
        (default is nucleotide sequences fies)")

## Add process reads command
process_reads = subparsers.add_parser("process_reads", aliases=['pr', 'p'], help="tool to filter and normalize reads prior to \
                                        transcriptome assembly")
assembly_reads = process_reads.add_argument_group("reads")
assembly_reads_mutual = assembly_reads.add_mutually_exclusive_group()
assembly_reads_mutual.add_argument("-s", "--single", help="specifies the file containing single-end reads")
assembly_reads_mutual.add_argument("-p", "--paired", nargs=2, metavar=("reads1", "reads2"), help="specifies the \
                                    paired-end reads files")
process_reads_output = process_reads.add_argument_group("output")
process_reads_output.add_argument("-o", "--output", help="changes the base directory for output files (default is current \
                    working directory)")
process_reads_output.add_argument("-n", "--name", default="result", help="specifies the process reads name (e.g col0) to organize \
                                    output files")
process_reads_config = process_reads.add_argument_group("configurations")
process_reads_config.add_argument("-q", "--q_value", type=int, default=20, help="changes q-value threshold for reads filtering \
                                  by erne-filter (default 20)")
process_reads_config.add_argument("-c", "--normCov", type=int, default=50, help="change the max kmer coverage for digital \
                            normalization; lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)")
process_reads.set_defaults(func=run_process_reads)

# Specify function to run
analysis.set_defaults(func=run_analysis)
## Parse args
args = parser.parse_args()
try:
    args.func(args)
except AttributeError:
    parser.print_help()
    parser.exit()