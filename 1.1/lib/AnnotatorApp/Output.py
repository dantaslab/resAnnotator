from AnnotatorApp.settings import *

class Output(object):
    """Class to format RGI output."""
    def __init__(self, output_dir, temp_dir, output_dict, db_params, prefix='resAnnotator'):
        """Creates Output object for formatting outputs."""
        self.output_dir = os.path.abspath(output_dir)
        self.temp_dir = os.path.abspath(temp_dir)
        self.prefix = prefix
        self.output_dict = output_dict
        self.db_params = db_params
        self.db_name = str(os.path.splitext(db_params['dbname'])[0])
        self.file_prefix = '{}_{}'.format(str(self.prefix), self.db_name)
        #self.output_filename = os.path.join(self.output_dir,"{}_{}_results_tab.txt".format(str(self.prefix),db_params['dbname']))
        #self.output_filepath = os.path.join(self.output_dir,self.file_prefix)
        #logger.debug('Created Output object')
        #logger.info(repr(self))

    def __repr__(self):
        """Returns Output class full object."""
        return "Output({}".format(self.__dict__)

    def run(self):
        """ Run output """
        #logger.info('Write annotated and unannotated sequences in fasta format')
        #self.write_fasta_files()

        logger.info("Create tab-delimited file.")
        if self.db_params['dbtype'].upper() == 'BLAST':
            self.write_blast_output()
        elif self.db_params['dbtype'].upper() == 'HMM':
            self.write_hmm_output()
        else:
            logger.error('invalid dbtype value for {}'.format(str(db_params['dbtype'])))


    def write_fasta_files(self):
        """ Write annotated sequences in fasta format"""
        file_list = []
        annot_fpath = os.path.join(self.output_dir, '{}_final_resistanceSeqs.fasta'.format(self.file_prefix))
        unannot_fpath = os.path.join(self.output_dir, '{}_final_unannotatedSeqs.fasta'.format(self.file_prefix))
        
        for filename in os.listdir(self.temp_dir):
            filepath = os.path.join(self.temp_dir, filename)
            if os.path.isfile(filepath) and filename.endswith('resistanceSeqs.fasta'):
                file_list.append(filepath)
            elif os.path.isfile(filepath) and filename.endswith('unannotatedSeqs.fasta'):
                #file_list.append(filepath) 
                pass
            else:
                pass 
        
        if len(file_list) > 0:
            with open(annot_fpath, 'w') as w_file:
                for filen in file_list:
                    with open(filen, 'rU') as o_file:
                        seq_records = SeqIO.parse(o_file, 'fasta')
                        SeqIO.write(seq_records, w_file, 'fasta')
        else:
            logger.warning("No annotated sequences found")
 

    def write_blast_output(self):
        """ Write blast output in tabular format."""
        if self.output_dict is not None:
            output_file = os.path.join(self.output_dir,'{}_blastout.txt'.format(self.file_prefix))
            with open(output_file, 'w') as OUT:
                OUT.write("QueryName\tQueryLength\tSubjectName\tSubjectSource\tSubjectLength\tAlignmentLength\tQueryStart\tQueryEnd\tSubjectStart\tSubjectEnd\tHspScore\tHspExpect\tPercentIdentity\tQueryCoverage\tSubjectCoverage\n")
            
                for qid, records in self.output_dict.items():
                    for hitid, result in records.items():
                        OUT.write("{query}\t{qlen}\t{sbjct}\t{source}\t{slen}\t{aln_len}\t{qstart}\t{qend}\t{sstart}\t{send}\t{score}\t{evalue}\t{ident}\t{qcov}\t{scov}\n".format(
                            query = str(result["qname"]),
                            qlen = str(result["qlen"]),
                            sbjct = str(result["sname"]),
                            source = str(result["source"]),
                            slen = str(result["slen"]),
                            aln_len = str(result["alignment_len"]),
                            qstart = str(result["qstart"]),
                            qend = str(result["qend"]),
                            sstart = str(result["sstart"]),
                            send = str(result["send"]),
                            score = str(result["bit_score"]),
                            evalue = str(result["evalue"]),
                            ident = str(result["perc_identity"]),
                            qcov = str(result["query_cov"]),
                            scov = str(result["sbjct_cov"])
                        ))



    def write_hmm_output(self):
        """ Write HMM output in tabula format."""
        if self.output_dict is not None:
            output_file = os.path.join(self.output_dir,'{}_hmmout.txt'.format(self.file_prefix))
            with open(output_file, 'w') as OUT:
                #header = ["target_name", "target_accession", "query_name", "query_accession", "seq_eval", 
                #"seq_score", "seq_bias", "dom_eval", "dom_score","dom_bias",
                #"exp_dom","num_reg", "num_clu", "num_overlap", "num_env", 
                #"num_dom", "num_rep", "num_inc", "description_of_target"]
 
                OUT.write("QueryName\tQueryAccession\tSubjectName\tSubjectSource\tSubjectAccession\tSubjectDescription\tSeqEvalue\tSeqScore\n")
                for qid, records in self.output_dict.items():
                    for hitid, result in records.items():
                        OUT.write("{qname}\t{qacc}\t{sname}\t{source}\t{sacc}\t{sdesc}\t{evalue}\t{score}\n".format(
                            qname = str(result["query_name"]),
                            qacc = str(result["query_accession"]),
                            sname = str(result["target_name"]),
                            source = str(result["source"]),
                            sacc = str(result["target_accession"]),
                            sdesc = str(result["description_of_target"]),
                            evalue = str(result["seq_eval"]),
                            score = str(result["seq_score"]),
                        ))
        else:
            logger.error("No result for {}".format(str(self.db_params['dbname'])))

    def checkKeyExisted(self,key, my_dict):
        try:
            nonNone = my_dict[key] is not None
        except KeyError:
            nonNone = False
        return nonNone


# obj = Output("./_tests/protein.fasta.blast.output.json")
# print(repr(obj))
# obj.run()

# obj = Output("./_tests/nucleotide.fa.diamond.output.json")
# print(repr(obj))
# obj.run()
