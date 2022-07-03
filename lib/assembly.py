import os

from typing import Any
from lib.util import Util as utillib


class ProcessReads():
    """Process reads prior to trsnscriptome assembly
    """
    def __init__(self, args, **kwargs) -> None:
        self.processname = args.name
        self.basedir = args.output if args.output else kwargs['basedir']
        self.output = os.path.join(self.basedir, f"process_reads/{self.processname}")
        self.logfile = os.path.join(self.output, "process_reads.log")
        self.paired = args.paired
        self.normvalue = args.normCov
        self.qvalue = args.q_value

        self.filtreadsdir = os.path.join(self.output, "filtered_reads")
        utillib.create_directory(self.filtreadsdir)
        self.normreadsdir = os.path.join(self.output, "normalized_reads")
        utillib.create_directory(self.normreadsdir)
        self.normreads = os.path.join(self.normreadsdir, "normalized_reads.fq")
 
        if self.paired:
            self.reads1 = args.paired[0]
            self.reads2 = args.paired[1]
            self.filtreads1 = os.path.join(self.filtreadsdir, "filtered_reads_1.fastq")
            self.filtreads2 = os.path.join(self.filtreadsdir, "filtered_reads_2.fastq")
            self.interreads = os.path.join(self.normreadsdir, "interleaved_reads.fq")
            self.extractedfile = f"{self.normreads}.pe"
            self.normreads1 = f"{self.extractedfile}.1"
            self.normreads2 = f"{self.extractedfile}.2"
        else:
            self.reads1 = args.single
            self.filtreads1 = os.path.join(self.filtreadsdir, "filtered_reads.fastq")

    def __filter_reads(self) -> None:
        """Filter raw reads
        """
        prefix = os.path.join(self.filtreadsdir, 'filtered_reads')
        emsg = "A system error occured when filtering reads"
        command = ""

        if self.paired:
            command = f"erne-filter --query1 {self.reads1} --query2 {self.reads2} --output-prefix {prefix} --ultra-sensitive \
                --force-standard"
            proc = utillib.execute_process(command=command, errormessage=emsg)
            if proc:
                print(f"Filter paired-end reads ran successfully. Output folder is {self.filtreadsdir}")
        else:
            command = f"erne-filter --query1 {self.reads1} --output-prefix {prefix} --ultra-sensitive --force-standard"
            proc = utillib.execute_process(command=command, errormessage=emsg)
            if proc:
                print(f"Filter single reads ran successfully. Output folder is {self.filtreadsdir}")

    def __interleave_reads(self) -> None:
        """Interleave paired end reads
        """
        intercommand = f"interleave-reads.py  -o {self.interreads} --force {self.filtreads1} {self.filtreads2}"
        emsg = "A system error occured when interleaving reads"

        proc = utillib.execute_process(command=intercommand, errormessage=emsg)
        if proc:
            print(f"Interleaving reads ran successfully. Output folder is {self.normreadsdir}")

    def __extract_paired_end_reads(self) -> None:
        """Extract paired end reads
        """
        command = f"extract-paired-reads.py -d {self.normreadsdir} {self.normreads}"
        emsg = "A system error occured when extracting paired reads"

        proc = utillib.execute_process(command=command, errormessage=emsg)

        if proc:
            print(f"Extraction of reads ran succeffully. Reads are stored in: {self.extractedfile}")

    def __split_paired_end_reads(self) -> None:
        """Split paired end reads
        """
        emsg = "A system error occured when splitting paired reads. Error details"
        command = f"split-paired-reads.py -d {self.normreadsdir} --force {self.extractedfile}"

        proc = utillib.execute_process(command=command, errormessage=emsg)

        if proc:
            print(f"Splitting of paired end reads ran succesfully.")

    def normalize_reads(self) -> Any:
        """Normalize reads

        Returns:
            tuple -- files containing normalized reads
        """
        
        # normalize command
        normcommand = ""
        emsg = "A system error occured when normalizing reads"

        # filter reads
        self.__filter_reads()

        if self.paired:
            # interleave reads
            self.__interleave_reads()
            print(f"Interleaved file is {self.interreads}")
            normcommand = f"normalize-by-median.py -o {self.normreads} -x 1000000000 -C {self.normvalue} --n_tables 32 \
                            --force --ksize 32 {self.interreads}"
            
            # run process
            proc = utillib.execute_process(command=normcommand, errormessage=emsg)
            if proc:
                print(f"Normalization of reads ran successfully. Output folder is {self.normreadsdir}")
                # extract paired-end reads
                self.__extract_paired_end_reads()
                # split extracted paired-end reads
                self.__split_paired_end_reads()
        else:
            normcommand = f"normalize-by-median.py -o {self.normreads} -x 1000000000 -C {self.normvalue} --n_tables 32 \
                            --force --ksize 32 {self.normreads1}"

            proc = utillib.execute_process(command=normcommand, errormessage=emsg)
            if proc:
                print(f"Normalization of reads ran successfully. Output folder is {self.normreadsdir}")