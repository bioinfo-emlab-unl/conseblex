import os
import csv
import re
import logging
import shlex, subprocess
from typing import Tuple


from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sys import maxsize
from lib.util import Util as utillib

class Analysis():
    """Performs analysis of the assemblies
    """
    def __init__(self, args, **kwargs) -> None:
        """Initialize Analysis class

        Arguments:
            args {Namespace} -- commandline arguments
        """
        self.analysisname = args.name
        self.write = args.write
        self.basedir = kwargs['basedir']
        self.output = os.path.join(self.basedir, 'analysis/{}'.format(self.analysisname))
        self.logfile = os.path.join(self.basedir, 'analysis/analysis.log')
        self.assemblynames = args.assemblies[0::2] # Get even indexes
        self.assemblyfiles = args.assemblies[1::2] # Get odd indexes
        self.assemblydict = dict(zip(self.assemblynames, self.assemblyfiles))
        self.benchmark = args.benchmark
        self.orfs = {}
        self.mergedtranscripts = {}
        self.mergedtranscripts_file = os.path.join(self.output,"mergedTranscripts.fa")
        self.transcriptsassemblers = {}
        self.longestorfs = {}
        self.longestmergedORFs_file = os.path.join(self.output, "longestMergedORFs.faa")
        self.sequenceclusters = {}
        self.contigsassemblers = {}
        self.overlapcombinations = {}

    def merge_files(self, **kwargs) -> None:
        """Merge given files into one file
        """

        for assemblyname, assemblyfile in self.assemblydict.items():
            # Read assembly file
            records = SeqIO.parse(assemblyfile, "fasta")
            addedseqs = set()
            for idx in records:
                seqid = idx.description
                sequence = idx.seq
    
                # Remove assembly name on seqid if already exists
                if(re.match(assemblyname, seqid, re.I)):
                    seqid = re.sub(assemblyname,"", seqid, flags=re.I)
                # Replace all spaces with underscore in seqid
                seqid = re.sub("[\[\]]", "", seqid)
                seqid = re.sub("[\s+,:]", "_", seqid)
                # substring seqid to 50 characters
                seqid = seqid[0:50]
                # Add to merged transcripts dictionary for later processing
                if sequence not in addedseqs:
                    self.mergedtranscripts[seqid] = sequence
                    self.transcriptsassemblers[seqid] = assemblyname.upper()
                    addedseqs.add(sequence)
        # Write to merged transcripts file  
        if self.write:
            utillib.write_id_and_sequence(data=self.mergedtranscripts, file=self.mergedtranscripts_file)

    def orf_finder(self, command: str, outfile: str, **kwargs) -> None:
        """Run ORFfinder to find open reading frames in given file

        Arguments:
            command {str} -- ORF command, with arguments, to run
            outfile {str} -- output file in ORF command
        """
        commandargs = shlex.split(command)
        if not os.path.isfile(outfile):
            try:
                # proc = subprocess.Popen(commandargs)
                # proc.wait()
                proc = True

                if proc:
                    # Use index file instead of creating dictionary
                    self.orfs = SeqIO.index(outfile, "fasta")
            except subprocess.CalledProcessError as err:
                pass
        else:
            self.orfs = SeqIO.index(outfile, "fasta")
        
    def get_longest_orfs(self, **kwargs) -> None:
        """Get longest open reading frame (ORF) for each Id in merged ORF file
        """
        if os.path.isfile(self.longestmergedORFs_file):
            records = SeqIO.index(self.longestmergedORFs_file, 'fasta')

            for idx in records:
                seqid = str(idx)
                sequence = str(records[idx].seq)

                self.longestorfs[seqid] = sequence
        else:
            # Read assembly file
            for idx in self.orfs:
                sqid = str(idx)
                sequence = str(self.orfs[idx].seq)
                # Split to get ID of sequence
                seqidinfo = re.split("ORF\d+_", sqid)
                info = re.split(":", seqidinfo[1])
                seqid = info[0] # Sequence Id after splits
                # Save to longest ORFs dictionary
                if seqid not in self.longestorfs:
                    # New entry in longest orfs
                    self.longestorfs[seqid] = sequence
                else: 
                    if len(sequence) > len(self.longestorfs[seqid]):
                        # Replace existing sequemce with longer sequence
                        self.longestorfs[seqid] = sequence
            # Write to longest merged ORFs file            
            if self.write:
                utillib.write_id_and_sequence(data=self.longestorfs, file=self.longestmergedORFs_file)

    def get_sequence_clusters(self, **kwargs) -> None:
        """Cluster sequences that have identical proteins
        """
        for seqid, sequence in self.longestorfs.items():
            # Add to sequences to match dictionary for later processing
            if sequence not in self.sequenceclusters:
                self.sequenceclusters[sequence] = [seqid]
            elif seqid not in self.sequenceclusters[sequence]:
                # if seqid not in self.sequenceclusters[sequence]:
                self.sequenceclusters[sequence].append(seqid)

    def __get_representative_sequences(self, seqids: list) -> Tuple[str, str, str, str, str, str]:
        """Get representative sequences used in the final output

        Arguments:
            seqids {list} -- sequence ids cluster

        Returns:
            Tuple[str, str, str, str, str, str] -- longest id, shortest id, longest reference id, \\
            shortest reference id, concantenated ids, reference assembly name 
        """
        longestlength = 0
        shortestlength = maxsize
        longestseqid, shortestseqid = "", ""
        # Initialize ids to use use in reference
        faaid, longfaid, shortfaid = "", "", ""

        for seqid in seqids:
            seqlength = len(self.mergedtranscripts[seqid])
            # Get Id of longest sequence
            if seqlength > longestlength:
                longestlength = seqlength
                longestseqid = seqid
            # Get Id of shortest sequence
            if seqlength < shortestlength:
                shortestlength = seqlength
                shortestseqid = seqid

            # Concantenate Ids with their respective assemblers
            longfaid = "{}-{}".format(self.transcriptsassemblers[longestseqid], longestseqid)
            shortfaid = "{}-{}".format(self.transcriptsassemblers[shortestseqid], shortestseqid)

            assemblername = self.transcriptsassemblers[seqid]

            if not faaid:
                faaid = "{}-{}".format(assemblername, seqid)
            else:
                faaid += " | {}-{}".format(assemblername, seqid)

        return longestseqid, shortestseqid, longfaid, shortfaid, faaid

    def get_consensus_sequences(self, **kwargs) -> Tuple[dict, dict, dict, dict, dict]:
        """Get overlapping sequences among assemblers

        Returns:
            Tuple[dict, dict, dict, dict, dict] -- longest nucleotides, shortest nucleotides, amino contigs, \\
            contigs reference, contigs assemblers
        """
        contignum = 1
        nucleotidelongest, nucleotideshortest, aminocontigs, contigsref = {}, {}, {}, {}
        for sequence, sequenceids in self.sequenceclusters.items():
            contiglabel = "contig-{}".format(contignum)
            # Get assembly overlap
            overlapset = set([self.transcriptsassemblers[seqid] for seqid in sequenceids])
            overlapset = sorted(overlapset)
            faaref = " | ".join(overlapset)
            overlap = len(overlapset)
            # Get representative sequences 
            longestseqid, shortestseqid, longfaid, shortfaid, faaid = self.__get_representative_sequences(sequenceids)

            ## Add contigs to respective dictionaries
            # Add longest sequence with sequence and overlap to final dictionary
            nucleotidelongest[longfaid] = {"seq": self.mergedtranscripts[longestseqid], "overlap": overlap}
            # Add shortest sequence with sequence and overlap to final dictionary
            nucleotideshortest[shortfaid] = {"seq": self.mergedtranscripts[shortestseqid], "overlap": overlap}
            # Add amono contig with sequence and overlap to final dictionary
            aminocontigs[contiglabel] = {"seq": sequence, "overlap": overlap}
            # Add contig label with refence id and overlap to final contigs dictionary
            contigsref[contiglabel] = {"id": faaid, "overlap": overlap}
            # Add contig label with assembly reference name and overlap to final assemblers dictionary
            self.contigsassemblers[contiglabel] = {"id": faaref, "overlap": overlap}
            # Add assembly overlap combination type e.g TRINITY | SPADES | SOAP to combinations dictionary
            self.overlapcombinations[faaref] = {"total": 0, "hits": 0}

            # Increment contig count
            contignum += 1

        return nucleotidelongest, nucleotideshortest, aminocontigs, contigsref

    def run_benchmark(self, contigs: dict) -> dict:
        """Compare contig set to benchmark transcriptome

        Arguments:
            contigs {dict} -- contig set

        Returns:
            dict -- hits and misses from comparison
        """
        assemblyhits = {}
        benchmarkids, benchmarkseqs = [], []
        benchmarkset = set()

        benchmark = SeqIO.index(self.benchmark, "fasta")

        # Create benchmark lists and set for comparison
        for record in benchmark:
            benchmarkids.append(record)
            benchmarkseqs.append(benchmark[record].seq)
            benchmarkset.add(benchmark[record].seq)

        # Set benchmark sequences total
        self.benchmarkseqtotal = len(benchmarkseqs)

        # Find matching sequences in contig set and benchmark
        for contig, seqitem in contigs.items():
            seq = seqitem["seq"]
            # Get target assembler overlap combination
            faaref = self.contigsassemblers[contig]["id"]
            if seq in benchmarkset:
                assemblyhits[contig] = {"id": benchmarkids[benchmarkseqs.index(seq)], "overlap": seqitem["overlap"]}
                # Increment hits count of combination in assembly overlap combination dictionary
                self.overlapcombinations[faaref]["hits"] += 1
            else:
                assemblyhits[contig] = {"id": "No hit", "overlap": seqitem["overlap"]}

            # Increment total count of combination in assembly overlap combination dictionary
            self.overlapcombinations[faaref]["total"] += 1

        return assemblyhits

    def calculate_statistics(self) -> None:
        totalhits, totalcontigs = 0, 0

        # Write to CSV file
        statisticsfile = os.path.join(self.output, "statistics.csv")
        with open(statisticsfile, 'w', newline='') as csvfile:
            statswriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_NONNUMERIC)
            statswriter.writerow(['Overlap', ' Total Contigs', 'True Positives', 'False Positives', 'False Negatives',\
                                'Total Benchmark', 'Precision', 'Recall', 'F1-Score'])

            for faaref, refitem in self.overlapcombinations.items():
                totalhits = refitem['hits']
                totalcontigs = refitem['total']

                # False negatives
                falseneg = self.benchmarkseqtotal - totalhits
                # False positives
                falsepos = totalcontigs - totalhits
                # Precision
                precision = round(totalhits / (totalhits + falsepos), 2)
                # Recall
                recall = round(totalhits / (totalhits + falseneg), 2)
                # F1 score
                fscore = round(2 * totalhits / ((2 * totalhits) + falsepos + falseneg), 2)


                statswriter.writerow([faaref, totalcontigs, totalhits, falsepos, falseneg, self.benchmarkseqtotal, \
                                    precision, recall, fscore])