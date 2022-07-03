import os
import re
import logging
import shlex
import subprocess
import warnings

from io import TextIOWrapper


class Util():
    def __init__(self) -> None:
        self.openfiles = []

    @staticmethod
    def write_id_and_sequence(data: dict, file: str, final: bool=False, overlap: int=0, plusonly: bool=True ) -> None:
        """Write sequence Id and sequence to file

        Arguments:
            data {dict} -- data to be written
            file {str} -- file handler
            final {bool} -- is final output file
            overlap {int} -- number of assembly overlaps
            plusonly {int} -- whether the output is overlap+ (i.e 1+, 2+, 3+ e.t.c) only
                          or includes overlap (i.e 1, 2, 3 e.t.c).
        """
        if not plusonly:
            # Get path
            loc, filename = os.path.split(file)
            # Get file extension
            ext = re.split('\.', filename)[-1]

            # Write to consemble single file
            consemble_single = os.path.join(loc, f"consemble{overlap}.{ext}")
            with open(consemble_single, 'w') as outfile:
                for seqid, sequence in data.items():
                    if sequence["overlap"] == overlap:
                        outfile.write(">{}\n".format(seqid))
                        outfile.write("{}\n".format(sequence["seq"]))
        
        # write to overlap+ file               
        with open(file, 'w') as outfile:
            for seqid, sequence in data.items():
                if final:
                    if sequence["overlap"] >= overlap:
                        outfile.write(">{}\n".format(seqid))
                        outfile.write("{}\n".format(sequence["seq"]))
                else:
                    outfile.write(">{}\n".format(seqid))
                    outfile.write("{}\n".format(sequence)) 

    @staticmethod
    def write_ref(data: dict, names: dict, file: str, overlap: int, plusonly: bool=True) -> None:
        """Write reference for contig Id against coresponding Ids and assembler names

        Arguments:
            data {dict} -- contig id and concatenated sequence id pairs
            names {dict} -- contig id and concatenated assembler names pairs
            file {str} -- file to be written to
            overlap {int} -- number of assembly overlaps
            plusonly {int} -- whether the output is overlap+ (i.e 1+, 2+, 3+ e.t.c) only
                          or includes overlap (i.e 1, 2, 3 e.t.c)
        """
        if not plusonly:
            # Get path
            loc, filename = os.path.split(file)
            # Get file extension
            ext = re.split('\.', filename)[-1]

            # Write to consemble single file
            ref_single = os.path.join(loc, f"consemble{overlap}_ref.{ext}")
            with open(ref_single, 'w') as outfile:
                for key, value in data.items():
                    if value["overlap"] == overlap:
                        outfile.write("{}\t".format(key))
                        outfile.write("\t\t{}\t".format(names[key]["id"]))
                        outfile.write("\t\t{}\n".format(value["id"]))

        # Write to overlap+ file                
        with open(file, 'w') as outfile:
            for key, value in data.items():
                if value["overlap"] >= overlap:
                    outfile.write("{}\t".format(key))
                    outfile.write("\t\t{}\t".format(names[key]["id"]))
                    outfile.write("\t\t{}\n".format(value["id"]))

    @staticmethod
    def write_benchmark_ref(data: dict, file: str, overlap: int=0, plusonly: bool=True) -> None:
        """Write benchmark reference to file

        Arguments:
            data {dict} -- contig and matching id in benckmark
            file {str} -- file to write to

        Keyword Arguments:
            overlap {int} -- number of assembly overlaps (default: {0})
            plusonly {int} -- whether the output is overlap+ (i.e 1+, 2+, 3+ e.t.c) only
                          or includes overlap (i.e 1, 2, 3 e.t.c) (default: {True})
            
        """
        if not plusonly:
            # Get path
            loc, filename = os.path.split(file)
            # Get file extension
            ext = re.split('\.', filename)[-1]

            # Write to consemble single file
            benchmark_single = os.path.join(loc, f"consemble{overlap}_hits.{ext}")
            with open(benchmark_single, 'w') as outfile:
                for id, value in data.items():
                    if value["overlap"] == overlap:
                        outfile.write("{}\t".format(id))
                        outfile.write("\t\t{}\n".format(value["id"]))

        # write to overlap+ file                
        with open(file, 'w') as outfile:
            for id, value in data.items():
                if value["overlap"] >= overlap:
                    outfile.write("{}\t".format(id))
                    outfile.write("\t\t{}\n".format(value["id"]))

    @staticmethod
    def directory_exits(directorypath: str) -> bool:
        """Check if directory exists

        Arguments:
            directorypath {str} -- directory path

        Returns:
            bool -- true is directory exists otherwise false
        """
        exists = os.path.exists(directorypath)
        return exists

    @staticmethod
    def create_directory(directorypath: str) -> None:
        """Create directory

        Arguments:
            directorypath {str} -- [description]
        """
        if not Util.directory_exits(directorypath):
            os.makedirs(directorypath)

    @staticmethod
    def configure_logging(logfile: str, level: int) -> None:
        """Configure looging

        Arguments:
            logfile {str} -- file to write logs to
            level {int} -- level of logging
        """
        logging.basicConfig(filename=logfile, level=level, 
                format='%(asctime)s.%(msecs)03d %(levelname)s {%(module)s} [%(funcName)s] %(message)s',
                datefmt='%Y-%m-%d,%H:%M:%S')

    def open_file(self, file: str, mode: str) -> TextIOWrapper:
        """Open given file

        Arguments:
            file {str} -- file name
            mode {str} -- mode in which file is opened e.g 'w', 'r', 'a'

        Returns:
            TextIOWrapper -- file handler
        """
        f = open(file, mode)
        self.openfiles.append(f)

        return f

    @staticmethod
    def execute_process(command: str, errormessage: str = "An error occured") -> subprocess.CompletedProcess:
        """Execute supprocess on host system

        Arguments:
            command {str} -- command to run
            errormessage {str} -- message to be printed if process fails

        Returns:
            subprocess.CompletedProcess -- _description_
        """
        normargs = shlex.split(command)
        print(f"Current command: {command}")
        try:
            proc = subprocess.Popen(normargs, stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate() # interact with sub process, wait for it to terminate
        except subprocess.CalledProcessError as err:
            warnings.warn(message=f"{errormessage}. Return code: {err.returncode}")
            warnings.warn(message=f"\n{err.stderr}")
            
        return proc