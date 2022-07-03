import os
import re
import shlex
import subprocess
import argparse
import time
import warnings

from typing import Any
from lib.util import Util

cwd = os.getcwd()
utils = Util()


def execute_process(command: str, errormessage: str = "An error occured") -> subprocess.CompletedProcess:
    """Execute supprocess on host system

    Arguments:
        command {str} -- command to run
        errormessage {str} -- message to be printed if process fails

    Returns:
        subprocess.CompletedProcess -- _description_
    """
    normargs = shlex.split(command)
    try:
        proc = subprocess.run(normargs, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as err:
        warnings.warn(message=f"{errormessage}. Return code: {err.returncode}")
        warnings.warn(message=f"\n{err.stderr}")
        
    return proc

def filter_reads(reads1: str, reads2: str = None, outputdir: str = cwd, paired: bool = True) -> tuple:
    """Filter raw reads

    Arguments:
        reads1 {str} -- first file containing reads. If paired = False, this is equivalent to single reads

    Keyword Arguments:
        reads2 {str} -- second file containg paired end reads. If paired = True, this is required (default: {None})
        outputdir {str} -- directory to which the result will be saved (default: {cwd})
        paired {bool} -- specify whether paired end reads or single reads (default: {True})
    """

    # create output directory
    filteredreads = os.path.join(outputdir, "filteredreads")
    utils.create_directory(filteredreads)

    prefix = os.path.join(filteredreads, "merged")
    emsg = "A system error occured when filtering reads"
    command = ""

    if paired:
        command = f"erne-filter --query1 {reads1} --query2 {reads2} --output-prefix {prefix} --ultra-sensitive \
            --force-standard"
        proc = execute_process(command=command, errormessage=emsg)
        if proc:
            print(f"Filter reads ran successfully. Output folder is {filteredreads}")
    else:
        command = f"erne-filter --query1 {reads1} --output-prefix {prefix} --ultra-sensitive --force-standard"
        proc = execute_process(command=command, errormessage=emsg)
        if proc:
            print(f"Filter reads ran successfully. Output folder is {filteredreads}")

    return f"{prefix}_1.fastq", f"{prefix}_2.fastq"

def interleave_reads(reads1: str, reads2: str, outputdir: str = cwd) -> str:
    """Interleave paired end reads

    Arguments:
        reads1 {str} -- first file containing reads
        reads2 {str} -- second file containing reads

    Keyword Arguments:
        outputdir {str} -- directory to save results in (default: {cwd})

    Returns:
        str -- full interleaved file path
    """
    # interleave reads
    # outname = re.split('.', reads1)[0]
    interfile = os.path.join(outputdir, "merged_interleaved.fq")

    intercommand = f"interleave-reads.py  -o {interfile} --force {reads1} {reads2}"
    emsg = "A system error occured when interleaving reads"

    proc = execute_process(command=intercommand, errormessage=emsg)
    if proc:
        print(f"Interleaving reads ran successfully. Output folder is {outputdir}")

    return interfile 

def extract_paired_end_reads(normfile: str, outputdir: str) -> str:
    """Extract paired end reads

    Arguments:
        normfile {str} -- normalized reads file
        outputdir {str} -- output directory

    Returns:
        str -- extracted paired ends full path
    """
    extractedfile = f"{normfile}.pe"
    # extractedfile = "merged.fq.pe"
    command = f"extract-paired-reads.py -d {outputdir} {normfile}"
    emsg = "A system error occured when extracting paired reads"

    proc = execute_process(command=command, errormessage=emsg)

    if proc:
        print(f"Extraction of reads ran succeffully. Reads are stored in: {extractedfile}")

    return extractedfile

def split_paired_end_reads(extractedfile: str, outputdir: str) -> tuple:
    """Split paired end reads

    Arguments:
        extractedfile {str} -- _description_
        outputdir {str} -- output directory

    Returns:
        tuple -- normlized paired-ends 1 and 2
    """
    # outname = re.split('.', extractedfile)[0]
    emsg = "A system error occured when splitting paired reads. Error details"
    command = f"split-paired-reads.py -d {outputdir} --force {extractedfile}"

    proc = execute_process(command=command, errormessage=emsg)

    if proc:
        print(f"Splitting of paired end reads ran succesfully.")

    return (f"{extractedfile}.1", f"{extractedfile}.2")

def normalize_reads(reads1: str, reads2: str = None, outputdir: str = cwd, paired: bool = True, 
                    normvalue: int = 50) -> Any:
    """Normalize reads

    Arguments:
        reads1 {str} -- first file containing reads. If paired = False, this acts as single reads

    Keyword Arguments:
        reads2 {str} -- second file containing reads (default: {None})
        outputdir {str} -- directory to save normalized reads (default: {cwd})
        paired {bool} -- whether paired or single reads (default: {True})
        normvalue {int} -- normalization parameter value (default: {50})

    Returns:
        tuple -- files containing normalized reads
    """
    
    # create output directory
    normalizedreads = os.path.join(outputdir, "normalizedreads")
    utils.create_directory(normalizedreads)

    # outname = re.split('.', reads1)[0]
    
    # normalize reads
    normcommand = ""
    normreads1, normreads2 = "", ""
    emsg = "A system error occured when normalizing reads"
    normfile = os.path.join(normalizedreads, "merged_normalized.fq")

    if paired:
        # interleave reads
        interfile = interleave_reads(reads1=reads1, reads2=reads2, outputdir=normalizedreads)
        print(f"Interleaved file is {interfile}")
        normcommand = f"normalize-by-median.py -o {normfile} -x 1000000000 -C {normvalue} --n_tables 32 \
                        --force --ksize 32 {interfile}"
        
        # run process
        proc = execute_process(command=normcommand, errormessage=emsg)
        if proc:
            print(f"Normalization of reads ran successfully. Output folder is {normalizedreads}")
            extractedfile = extract_paired_end_reads(normfile=normfile, outputdir=normalizedreads)
            normreads1, normreads2 = split_paired_end_reads(extractedfile=extractedfile, outputdir=normalizedreads)
    else:
        normcommand = f"normalize-by-median.py -o {normfile} -x 1000000000 -C {normvalue} --n_tables 32 \
                        --force --ksize 32 {reads1}"

        _ = execute_process(command=normcommand, errormessage=emsg)

    return (normreads1, normreads2) if paired else normfile


if __name__ == '__main__':
    # Process commandline arguments
    parser = argparse.ArgumentParser(description="Helper tool to normalize reads prior to transcriptome assembly")
    assembly_reads = parser.add_argument_group("reads")
    assembly_reads_mutual = assembly_reads.add_mutually_exclusive_group()
    assembly_reads_mutual.add_argument("-s", "--single", help="specifies the file containing single-end reads")
    assembly_reads_mutual.add_argument("-p", "--paired", nargs=2, metavar=("reads1", "reads2"), help="specifies the \
                                        paired-end reads files")
    parser.add_argument("-o", "--output", help="changes the base directory for output files (default is current \
                        working directory)")
    args = parser.parse_args()
    
    if args.paired:
        reads1 = args.paired[0]
        reads2 = args.paired[1]
        # Filter reads
        filreads1, filreads2 = filter_reads(reads1=reads1, reads2=reads2, outputdir=args.output)
        # Normalize reads
        norm1, norm2 = normalize_reads(reads1=filreads1, reads2=filreads2, outputdir=args.output)
    else:
        # Filter reads
        filreads1, _ = filter_reads(reads1=args.single, outputdir=args.output, paired=False)
        # Normalize reads
        norm1 = normalize_reads(reads1=filreads1, outputdir=args.output, paired=False)