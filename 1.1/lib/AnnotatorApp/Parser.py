
from AnnotatorApp.settings import *
from Bio import SeqIO
from Bio.Blast import NCBIXML
import collections
from collections import OrderedDict

class Parser(object):
    """Class to create a parser for blast and hmm output"""

    def __init__(self, input_sequence, result_file, input_type=None, temp_dir=None, db_params=None): 
        
        self.input_sequence = input_sequence
        self.result_file = result_file
        self.db_params = db_params
        self.input_type = input_type

        if temp_dir == None:
            self.temp_dir = os.path.split(inputfile)[0]
        else:
            self.temp_dir = temp_dir
        
        self.dbname = self.db_params['dbname']

        #self.pass_bitscore = 5
        #self.pass_evalue = float(0.001)
        #self.pass_query_cov = 30.0
        #self.pass_sbjct_cov = 30.0
        #self.pass_identity = 50.0
        #self.max_target_seq = 1

       
    def __repr__(self):
        """Returns Parser class full object."""
        return "Parser({})".format(self.__dict__)

    def run(self):
        """Parse blast/hmm output file"""
        output_dict = None
        final_outdict = None
        
        hit_cnts = self.db_params.get('max_target_seq', 1)
        logger.debug('hit_cnts = {}'.format(str(hit_cnts)))
        sortby_col = None

        if self.db_params['dbtype'].upper() == 'BLAST':
            output_dict = self.parse_xml()
            sortby_col = 'bit_score'
        elif self.db_params['dbtype'].upper() == 'HMM':
            output_dict = self.parse_hmm()
            sortby_col = 'seq_score'
        else:
            logger.error('invalid dbtype value for {}'.format(str(db_params['dbtype'])))
       
        if output_dict:  
            final_outdict = self.filter_by_hitCounts(output_dict, sortby_col, hit_cnts)

        return final_outdict

    def parse_xml(self):
        """Parse blast output in xml file"""
        logger.debug("Parser({})".format(self.__dict__))

        blastOutput = {}
        xml_file = self.result_file

        for key, val in self.db_params.items():
            if key == 'query_cov' and 0 < val <= 100:
                pass_query_cov = float(val)
            elif key == 'subject_cov' and 0 < val <= 100:
                pass_sbjct_cov = float(val)
            elif key == 'identity_cutoff' and 0 < val <= 100:
                pass_identity = float(val)
            elif key == 'evalue' and 0 <= float(val) <= 0.001:
                pass_evalue = float(val)
            elif key == 'max_target_seqs' and isinstance(val, int) == True:
                max_tagets_seqs = int(val)

        if not os.path.exists(xml_file) or os.stat(xml_file).st_size == 0:
            logger.error('{} does not exists'.format(xml_file))
            return None
        
        with open(xml_file, 'r') as result_handle:
            
            blast_records = NCBIXML.parse(result_handle)
            init = 0
            
            for blast_record in blast_records:
                blast_results = {}
                
                for alignment in blast_record.alignments:
                    alignTitle = alignment.title
                    spacepos = alignTitle.index(' ')
                    hitid = alignTitle[0:spacepos]
                    hitid = hitid.encode('ascii','replace')
                    
                    for hsp in alignment.hsps:
                        hspdict = {}
                        querySeq = hsp.query.replace('-', '')
                        sbjctSeq = hsp.sbjct.replace('-', '')
                        
                        if self.input_type.lower() == 'protein' and self.db_params['seqtype'].lower() == 'nucl':
                            realSubjectLength = int(len(sbjctSeq)*3)
                        else:
                            realSubjectLength = len(sbjctSeq)
                        
                        if self.input_type.lower() == 'contig' and self.db_params['seqtype'].lower() == 'prot':
                            realQueryLength = int(len(querySeq)*3) # this includes matches and mis-matches
                        else:
                            realQueryLength = int(len(querySeq)) # this includes matches and mis-matches
                        
                        hspdict["alignment_len"] = hsp.align_length
                        hspdict["source"] = str(os.path.splitext(self.dbname)[0])
                        hspdict["evalue"] = hsp.expect
                        hspdict["bit_score"] = hsp.bits
                        hspdict["max_identities"] = hsp.identities
                        #hspdict["match"] = hsp.match
                        #hspdict["alignment_len"] = hsp.align_length
                        # percent identity is calculated w.r.t aligned region
                        hspdict["perc_identity"] = float(format(float(hspdict["max_identities"]*100) / int(hspdict["alignment_len"]), '.2f'))
                        hspdict["qname"] = blast_record.query
                        #hspdict["qseq"] = hsp.query
                        hspdict["qstart"] = hsp.query_start
                        hspdict["qend"] = hsp.query_end #hsp.query_start + realQueryLength
                        hspdict["qlen"] = blast_record.query_length
                        # query coverage is the num of bases in the query sequence aligned with subject seq
                        hspdict["query_cov"] = float(format(float(realQueryLength*100) / int(hspdict["qlen"]), '.2f'))        
                        #hspdict["query_cov"] = float(format(float(hspdict["alignment_len"]*100) / int(hspdict["qlen"]), '.2f'))     
                        hspdict["sname"] = str(alignment.title)                           
                        hspdict["sseq"] = str(hsp.sbjct)
                        hspdict["sstart"] = str(hsp.sbjct_start)
                        hspdict["send"] = str(hsp.sbjct_end)
                        hspdict["slen"] = str(alignment.length)
                        #hspdict["sbjct_cov"] = float(format(float(hspdict["alignment_len"]*100) / int(hspdict["slen"]), '.2f')) 
                        hspdict["sbjct_cov"] = float(format(float(realSubjectLength*100) / int(hspdict["slen"]), '.2f'))        
                        
                        if float(hspdict["evalue"]) <= float(pass_evalue) \
                            and float(hspdict["perc_identity"]) >= float(pass_identity) \
                            and float(hspdict["query_cov"]) >= float(pass_query_cov) \
                            and float(hspdict["sbjct_cov"]) >= float(pass_sbjct_cov):
                            
                            blast_results["{}".format(hitid.decode())] = hspdict
                                        
                if blast_results:
                    blastOutput[blast_record.query] = blast_results
                
        return blastOutput
    

    def parse_hmm(self):
        """ Parse hmm tblout format """
        logger.debug("Parser({})".format(self.__dict__))
        hmm_outfile = self.result_file
        output_dict = collections.defaultdict(dict)
        
#        pass_bitscore = 5
#        pass_evalue = float(0.001)
#        pass_query_cov = 30.0
#        pass_sbjct_cov = 30.0
#        pass_identity = 50.0


        header = ["target_name", "target_accession", "query_name", "query_accession", "seq_eval", 
                "seq_score", "seq_bias", "dom_eval", "dom_score","dom_bias",
                "exp_dom","num_reg", "num_clu", "num_overlap", "num_env", 
                "num_dom", "num_rep", "num_inc", "description_of_target"]
    
        with open(hmm_outfile, 'r') as fhandle:
            for line in fhandle:
                record = {}
                
                if line.startswith("#") or line is "":
                    continue
                line = line.rstrip()
                #line = split("\s+", line, 18)
                itemlist = line.split(None, 18)
                for idx, item in enumerate(itemlist):
                    record[header[idx]] = item
               
                record["source"] = str(os.path.splitext(self.dbname)[0])
                output_dict[itemlist[2]][itemlist[1]] = record
        
        return output_dict


    def filter_sequences(self, output_dict):
        """ filter fasta sequences based on blast results"""
        logger.debug("Parser => {}".format(str(output_dict)))
        result_fasta = {}
        unannotated_fasta = {}
        seq_filepath = os.path.abspath(self.input_sequence)
    
        if os.stat(seq_filepath).st_size != 0 and len(output_dict) > 0:
    
            for record in SeqIO.parse(seq_filepath, 'fasta'):
                if str(record.description) in output_dict:
                    result_fasta[record.id] = record
                else:
                    unannotated_fasta[record.id] = record
     
        return result_fasta, unannotated_fasta
    

    def filter_by_hitCounts(self, output_dict, sortby_col, hit_cnts=1):
        """ The function sorts the dictionary by bits_score and keep hits based
            on 'max_target_seqs' value set in the configuration file.
            By default, function will keep only the top hit for each query,
            But, if max_target_seqs is set to 0, then it will return all hits 
            sorted by bits_score
        """ 
        final_dict = OrderedDict() 
        hitCount = int(hit_cnts)

        for qid in output_dict:
            count = 1
            final_dict[qid] = OrderedDict()
            #print("{}---->".format(str(qid)))
            for hitid in sorted(output_dict[qid], key=lambda x: float(output_dict[qid][x][sortby_col]), reverse=True):
                if (count <= hitCount) or (hitCount <= 0):
                    #print("\t\t{} => {}".format(str(hitid), str(output_dict[qid][hitid]['seq_score'])))
                    final_dict[qid][hitid] = output_dict[qid][hitid]
                    count += 1
                else:
                    break

        return final_dict 

   
#    def write_results_to_file(outfile, blast_results = None):
#            
#        with open(outfile, 'w') as OUT:
#            OUT.write("Query Name\tQuery Length\tSubject Name\tSubject Length\tAlignment Length \
#                        \tQuery Start\tQuery End\tSubject Start\tSubject End\tHsp Score\tHsp Expect \
#                        \tPercent Identity\tQuery Coverage\tSubject Coverage\n")
#            
#            if blast_results is not None:
#                for qid, records in blast_results.items():
#                    for hitid, result in records.items():
#                        OUT.write("{query}\t{qlen}\t{sbjct}\t{slen}\t{aln_len}\
#                                    \t{qstart}\t{qend}\t{sstart}\t{send}\t{score}\t{evalue}\
#                                    \t{ident}\t{qcov}\t{scov}\n"\
#                            .format(
#                                    query = str(result["qname"]),
#                                    qlen = str(result["qlen"]),
#                                    sbjct = str(result["sname"]),
#                                    slen = str(result["slen"]),
#                                    aln_len = str(result["alignment_len"]),
#                                    qstart = str(result["qstart"]),
#                                    qend = str(result["qend"]),
#                                    sstart = str(result["sstart"]),
#                                    send = str(result["send"]),
#                                    score = str(result["bit_score"]),
#                                    evalue = str(result["evalue"]),
#                                    ident = str(result["perc_identity"]),
#                                    qcov = str(result["query_cov"]),
#                                    scov = str(result["sbjct_cov"])
#                            ))
#
#
#with open(xml_file, 'r') as result_handle:
#    blast_records = NCBIXML.parse(result_handle)
#    for blast_record in blast_records:  
#        for alignment in blast_record.alignments:
#            for hsp in alignment.hsps:
#                OUT.write("{query}\t{qlen}\t{sbjct}\t{slen}\t{aln_len}\t{qstart}\t{qend}\t{sstart}\t{send}\t{score}\t{evalue}\t{ident}\n"\
#                    .format(
#                        query = str(blast_record.query),
#                        qlen = str(blast_record.query_length),
#                        sbjct = str(alignment.title),
#                        slen = str(alignment.length),
#                        aln_len = str(hsp.align_length),
#                        qstart = str(hsp.query_start),
#                        qend = str(hsp.query_end),
#                        sstart = str(hsp.sbjct_start),
#                        send = str(hsp.sbjct_end),
#                        score = str(hsp.score),
#                        evalue = str(hsp.expect),
#                        ident = str(hsp.identities)
#
#                    ))
#                cnt += 1
