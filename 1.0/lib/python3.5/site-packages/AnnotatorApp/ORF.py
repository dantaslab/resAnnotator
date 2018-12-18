from AnnotatorApp.settings import os, SeqIO, logger

class ORF(object):
    """Class to find open reading frames from nucleotide sequence."""
    def __init__(self,input_file, clean=True, temp_dir=None, working_directory=None, low_quality=False):
        """Creates ORF object for finding open reading frames."""
        self.input_file = input_file
        self.clean = clean
        self.working_directory = working_directory
        self.temp_dir = temp_dir
        self.low_quality = low_quality

    def __repr__(self):
        """Returns ORF class full object."""
        return "ORF({}".format(self.__dict__)

    def contig_to_orf(self):
        """Converts contigs to open reading frames."""
        self.orf_prodigal()

    def min_max_sequence_length(self):
        """Returns minimum and maximun sequence length in multi-fasta inputs""" 
        sequences = []
        for record in SeqIO.parse(self.input_file, "fasta"):
            sequences.append(len(record.seq))
        return min(sequences), max(sequences)

    def orf_prodigal(self):
        """Runs PRODIGAL to find open reading frames."""
        quality = "-n -p single"

        _min, _max = self.min_max_sequence_length()
        logger.info("minimum sequence length: {}, maximun sequence length {}".format(_min,_max))

        if self.low_quality == True or _min < 20000:
            quality = "-p meta"

        filename = os.path.splitext(os.path.basename(self.input_file))[0]

        stdout = "2> /dev/null"

        cmd = "prodigal -q -m -a {trans_file} -i {input_file} -o {output_file} -d {nuc_file} -s {potential_genes} {quality} {stdout}" \
           .format(
                trans_file=os.path.join(self.temp_dir, "{}.contig.fsa".format(filename)),
                input_file=self.input_file,
                output_file=os.path.join(self.temp_dir, "{}.draft".format(filename)),
                quality=quality,
                stdout=stdout,
                nuc_file= os.path.join(self.temp_dir, "{}.contigToORF.fsa".format(filename)),
                potential_genes= os.path.join(self.temp_dir, "{}.potentialGenes".format(filename))
            )

        logger.info(cmd)
        os.system(cmd)

        # format the contig file headers to remove after first space
        trans_file = os.path.join(self.temp_dir, "{}.contig.fsa".format(filename))
        self.format_fastaheaders(trans_file)
        nuc_file = os.path.join(self.temp_dir, "{}.contigToORF.fsa".format(filename))
        self.format_fastaheaders(nuc_file)

        if self.clean == True:
            os.remove(os.path.join(self.working_directory, "{}.draft".format(filename)))


    def get_character_len(self,file_path):
        """Returns character count in a file."""
        chars = words = lines = 0
        with open(file_path, 'r') as in_file:
            for line in in_file:
                if line[0] == '>':
                    pass
                else:
                    lines += 1
                    words += len(line.split())
                    chars += len(line)
        # logger.info("chars count: {}".format(chars))
        return chars

    def format_fastaheaders(self, filepath):
        fpath, fname = os.path.split(filepath)
        fname_prefix, fname_ext = os.path.splitext(fname)
        outfile = os.path.join(fpath, "{}_fixedHeader{}".format(fname_prefix, fname_ext))
 
        if os.stat(filepath).st_size > 0:
            with open(outfile, 'w') as fout:
               for record in SeqIO.parse(filepath, 'fasta'):
                    header = record.id.rstrip()
                    seq = record.seq
                    fout.write(">{}\n{}\n".format(header, seq))
        else:
            logger.error('Contig file is empty!!')
