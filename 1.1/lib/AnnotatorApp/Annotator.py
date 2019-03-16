from Bio import SeqIO
import glob
import time
import yaml
import shutil
import filetype
import tempfile

from AnnotatorApp.settings import *
from AnnotatorApp.Database import Database
from AnnotatorApp.Blast import Blast
from AnnotatorApp.Hmm import Hmm
from AnnotatorApp.Diamond import Diamond
from AnnotatorApp.Parser import Parser
from AnnotatorApp.Output import Output
from AnnotatorApp.ORF import ORF

class Annotate:
    """Class to predict resistome(s) from protein or nucleotide data based on BLAST or HMM based AR databases such as CARD and Resfams"""

    def __init__(self, input_type='contig', input_sequence=None, config_file=None, prefix=None, output_dir=None, \
                clean=True, num_threads=4, alignment_tool=None, log_level='DEBUG', general=None, indexes=None, order=None):

        self.input_type = input_type.lower()
        self.input_sequence = os.path.abspath(input_sequence)
        self.config_file = os.path.abspath(config_file)
        self.general = general
        self.indexes = indexes
        self.order = order
        self.aligner = alignment_tool.lower()
        
        self.num_sequences = 1
        self.output_path = os.path.abspath(output_dir)
        self.prefix = str(prefix)
        self.clean = clean
        self.num_threads = num_threads

        self.working_directory = self.output_path
        self.results_file = ''        
        self.output_dict = {} 
        
        self.idx = 1
        self.db_params = None
        self.log_level = log_level.upper()

    
    def __repr__(self):
        """Returns RGI class full object."""
        return "RGI({}".format(self.__dict__)


    def run(self):
        """Runs ResAnnotator Pipeline"""
        t0 = time.time()
        logger.info("Started the ResAnnotator Pipeline v{}".format(str(SOFTWARE_VERSION)))
        self.validate_inputs()

        for idx, dbname in enumerate(self.order, 1):
            logger.info("=====Working with database {}=====".format(str(dbname)))
            print("{} Working with database => {}".format(str(idx), str(dbname)))
            logger.debug('dbname {} => {}'.format(str(idx), str(dbname)))
            logger.debug('input_file => {}'.format(str(self.input_sequence)))
            
            db_params = self.indexes.get(dbname, None)
            flag_setdb = self.set_current_db(idx, db_params)    
            
            if flag_setdb:
                self.create_databases()
                self.align_sequences()
                self.filter_process()
                self.filter_input_sequences()
                self.write_output_file()

        logger.info('Total running time {}s'.format(round(time.time() - t0, 3)))
        logger.info('ResAnnotator pipeline ended!!')
        print("ResAnnotator pipleine completed successfully!!!")

    def validate_inputs(self):
        """Validate inputs."""
        #print("Output files will be written to: " + self.output_path)        
        self.temp_dir = tempfile.mkdtemp(prefix="temp_", dir=os.path.abspath(self.output_path))
        print("Temporary files will be written to:{}".format(str(self.temp_dir)))

        if not os.path.exists(self.input_sequence):
            logger.error("input file does not exist: {}".format(self.input_sequence))
            exit()

        logger.info("{} => {}".format(self.input_sequence, filetype.guess(self.input_sequence)))
        kind = filetype.guess(self.input_sequence)
	#print("File Type kind = {}".format(str(kind)))
        header_flag = False
        if kind is None:
            if is_fasta(self.input_sequence):
                logger.debug("Input sequence=> {} is in Fasta format".format(os.path.basename(self.input_sequence)))
                logger.info('Checking fasta headers for white spaces')
                header_flag = check_headers(self.input_sequence)
                
                dst_fname = str(os.path.splitext(os.path.basename(self.input_sequence))[0]) + '.' + str(self.idx) + '.fasta'
                dst_fpath = os.path.join(self.temp_dir, dst_fname)

                if not header_flag:
                    logger.info("Removing white spaces from fasta headers")
                    count = fix_headers(self.input_sequence, dst_fpath)
                    logger.info("Copied {} to temp dir {}".format(str(os.path.basename(self.input_sequence)), str(dst_fpath)))
                    logger.info("We have found {} sequences in {}".format(str(count), str(os.path.basename(dst_fname))))
                else:
                    shutil.copy(self.input_sequence, dst_fpath)
                    logger.info("Copied {} to temp dir {}".format(str(os.path.basename(self.input_sequence)), str(dst_fpath)))
                
                self.set_input_filepath(dst_fpath)
        
            logger.info("Validating inputs and sequence files....Completed")
            print("Validating inputs and sequence files....Completed\n")

            print("Identifying Open Reading Frames in contig file using prodigal")
            input_seqfile = self.get_input_filepath()
            f_path, f_name = os.path.split(input_seqfile)
            f_name = os.path.splitext(f_name)[0]
            if self.input_type == 'contig':
                orf_obj = ORF(input_file=input_seqfile, clean=self.clean, temp_dir=self.temp_dir, working_directory=self.working_directory)
                orf_obj.contig_to_orf()
                contig_fsa_file = os.path.join(self.temp_dir,"{}.contigToORF_fixedHeader.fsa".format(f_name))	
                print("Input orf file = {}".format(str(os.path.basename(contig_fsa_file))))
                self.set_input_filepath(contig_fsa_file)

        else:
            logger.error(kind.extension)
            logger.error(kind.mime)
            logger.warning("Sorry, no support for this format.")
            exit()
   

    def create_databases(self):
        """Creates databases."""
        #logger.info("----------Build database index fies----------")
        idx, db_params = self.get_current_db()
        logger.info("-----Step {}.1 Building Database => {}-----".format(str(idx), db_params['dbname']))
        print("{}.1 Building Database => {}".format(str(idx), db_params['dbname']))
        db_obj = Database(db_params, self.temp_dir, self.aligner)
        db_obj.build_databases()
      
 
    def align_sequences(self):
        """Align input sequences to list of databases"""        
        idx, db_params = self.get_current_db()
        logger.info("-----Step {}.2 Aligning sequences againsts {}-----".format(str(idx), db_params['dbname']))
        print("{}.2 Aligning sequences againsts {}".format(str(idx), db_params['dbname']))
        if db_params['dbtype'].upper() == 'BLAST':
            self.run_blast(db_params)
        
        elif db_params['dbtype'].upper() == 'HMM':
            self.run_hmmscan(db_params)
 
        else:
            logger.error('invalid dbtype value for {}'.format(str(db_params['dbtype'])))
           

    def run_blast(self, db_params):
        """Runs blast."""
        input_seqfile = self.get_input_filepath()
        input_seqtype = self.input_type
        f_path, f_name = os.path.split(input_seqfile)
        f_name = os.path.splitext(f_name)[0] 
        db_name = os.path.splitext(db_params['dbname'])[0]
        #if input_seqtype == 'contig':
        #    orf_obj = ORF(input_file=input_seqfile, clean=self.clean, temp_dir=self.temp_dir, working_directory=self.working_directory)
        #    orf_obj.contig_to_orf()
             #if db_params['seqtype'].lower() == 'nucl':
             #    contig_fsa_file = os.path.join(self.temp_dir,"{}.contigToORF_fixedHeader.fsa".format(f_name))
             #    self.set_input_filepath(contig_fsa_file)
             #else:
             #    contig_fsa_file = os.path.join(self.temp_dir,"{}.contig_fixedHeader.fsa".format(f_name))
             #    self.set_input_filepath(contig_fsa_file)
             #    input_seqtype = 'protein'
            
        #input_seqfile = self.get_input_filepath()
        #input_seqfile = contig_fsa_file
        xml_file = os.path.join(self.temp_dir,"{}_{}.xml".format(f_name, str(db_name)))
        logger.debug("blast xml_file = {}".format(xml_file))
        self.set_results_filepath(xml_file)

        try:
            if os.stat(input_seqfile).st_size > 0:
                if self.aligner == "diamond" and db_params['seqtype'].lower()=='prot':
                    logger.debug('Align sequences using DIAMOND')
                    diamond_obj = Diamond(input_seqfile, input_seqtype, self.results_file, self.temp_dir, db_params, self.num_threads)
                    diamond_obj.run()
                else:
                    logger.debug('Align sequences using BLAST')
                    blast_obj = Blast(input_seqfile, input_seqtype, self.results_file, self.temp_dir, db_params, self.num_threads)
                    blast_obj.run()
            else:
                self.write_stub_output_file()
                logger.error("no open reading frames (orfs) found.")
        except Exception as e:
            self.write_stub_output_file()
            logger.exception("failed to write orf file")
        else:
            # logger.info("success procession orf file")
            pass


    def run_hmmscan(self,db_params):
        """Run hmmscan"""
        #logger.debug("Running Hmmer")
        input_seqfile = self.get_input_filepath()
        f_path, f_name = os.path.split(input_seqfile)
        f_name = os.path.splitext(f_name)[0] 
        db_name = os.path.splitext(db_params['dbname'])[0]
        hmm_out = os.path.join(self.temp_dir,"{}_{}_hmmout.txt".format(f_name, str(db_name)))
        logger.debug("results_hmm_file = {}".format(hmm_out))
        self.set_results_filepath(hmm_out)

        if self.input_type == 'protein': 
            hmm_obj = Hmm(input_seqfile, self.results_file, self.temp_dir, db_params, self.num_threads)
            hmm_obj.run()
        elif self.input_type == 'contig':
            orf_obj = ORF(input_file=input_seqfile, clean=self.clean, temp_dir=self.temp_dir, working_directory=self.working_directory)
            orf_obj.contig_to_orf()
            contig_fsa_file = os.path.join(self.temp_dir,"{}.contig_fixedHeader.fsa".format(f_name))
            logger.debug("Inside HMM function Input orf file = {}".format(str(os.path.basename(contig_fsa_file))))	
            #print("Inside HMM function Input orf file = {}".format(str(os.path.basename(contig_fsa_file))))
            self.set_input_filepath(contig_fsa_file)
            input_seqfile = self.get_input_filepath()
            logger.debug("input file name= {}".format(str(input_seqfile)))
            #orf_obj = ORF(input_file=input_seqfile, clean=self.clean, working_directory=self.working_directory)
            #orf_obj.contig_to_orf()
            #contig_fsa_file = os.path.join(self.temp_dir,"{}.contig_fixedHeader.fsa".format(f_name))
            #self.set_input_filepath(contig_fsa_file)
            hmm_obj = Hmm(input_seqfile, self.results_file, self.temp_dir, db_params, self.num_threads)
            hmm_obj.run()
        else:
            logger.error("Invalid input_type!!")
            sys.exit(1)

    
    def filter_process(self):
        """ Parse output files and then filter the input sequences"""
        
        idx, db_params = self.get_current_db()
        output_dict = None
        input_seqfile = self.get_input_filepath()
        logger.debug('input_file => {}'.format(str(input_seqfile)))
        
        logger.info("-----Step {}.3 parsing xml output {}-----".format(str(idx),str(os.path.basename(self.results_file))))
        print("{}.3 parsing xml output {}".format(str(idx),str(os.path.basename(self.results_file))))
        parser_obj = Parser(input_seqfile, self.results_file, self.input_type, self.temp_dir, db_params)
        output_dict = parser_obj.run()
 
        self.set_output_dict(output_dict)


    def filter_input_sequences(self):
        """ Write annotated and un-annotated sequences in fasta format """ 
        input_seqfile = self.get_input_filepath()
        f_path, f_name = os.path.split(input_seqfile)
        current_db_idx, db_params = self.get_current_db()
        output_dict = self.get_output_dict()
        sep_idx = f_name.index('.')
        unannotated_fname = '{}.{}.unannotatedSeqs.fasta'.format(f_name[:sep_idx], str(current_db_idx + 1))
        annotated_fname = '{}.{}.annotatedSeqs.fasta'.format(f_name[:sep_idx], str(current_db_idx))
        logger.info("-----Step {}.4 Filtering sequences => {}-----".format(str(current_db_idx), str(os.path.basename(input_seqfile))))
        print("{}.4 Filtering sequences => {}".format(str(current_db_idx), str(os.path.basename(input_seqfile))))
        parser_obj = Parser(input_seqfile, self.results_file, self.input_type, self.temp_dir, db_params)

        if output_dict is not None:
            filteredseqs, unannotatedseqs = parser_obj.filter_sequences(self.output_dict)
            annotated_fpath = os.path.join(self.output_path, annotated_fname)
            unannotated_fpath = os.path.join(self.temp_dir, unannotated_fname)

            with open(annotated_fpath, 'w') as faa_handle:
                for id, seq in filteredseqs.items():
                    SeqIO.write(seq, faa_handle, "fasta")

            with open(unannotated_fpath, 'w') as faa_handle:
                for id, seq in unannotatedseqs.items():
                    SeqIO.write(seq, faa_handle, "fasta")
                
            if os.path.exists(unannotated_fpath) and os.stat(unannotated_fpath).st_size > 0:
                #logger.debug('unannotated_fasta exists => {}'.format(str(unannotated_fasta)))
                logger.debug('setting input_file to => {}'.format(str(unannotated_fpath)))
                self.set_input_filepath(unannotated_fpath)
            else:
                logger.debug('unannotated_fasta file does not exists')


    def write_output_file(self): 
        #logger.info("----------Write output----------")
        output_dict = self.get_output_dict()
        idx, db_params = self.get_current_db()
        logger.info("-----Step {}.5 Writing output to a file-----".format(str(idx)))
        print("{}.5 Writing output to a file".format(str(idx)))
        out_obj = Output(self.output_path, self.temp_dir, output_dict, db_params, prefix=self.prefix)
        out_obj.run()


    def clean_files(self):
        """Cleans temporary files."""
        if self.clean == True:
            basename_output_file = os.path.splitext(os.path.basename(self.output_file))[0]
            logger.info("Cleaning up temporary files...{}".format(basename_output_file))
            # clean working_directory
            self.clean_directory(self.working_directory, basename_output_file)
            d_name, f_name = os.path.split(self.output_file)
            # clean destination_directory
            self.clean_directory(d_name, basename_output_file)
            #shutil.rmtree(self.tempdir)
        else:
            logger.info("Clean up skipped.")


    def clean_directory(self, directory, basename_output_file):
        """Cleans files in directory."""
        logger.info(directory)
        files = glob.glob(os.path.join(directory, "*"))
        for f in files:
            if os.path.basename(self.input_sequence) + ".temp" in f and os.path.isfile(f):
                self.remove_file(f)

    def remove_file(self, f):
        """Removes file."""
        if os.path.exists(f):
            logger.info("Removed file: {}".format(f))
            os.remove(f)
        else:
            logger.warning("Missing file: {}".format(f))

    def out(self):
        """Writes tab-delimited, output files."""
        tab_obj = Output(self.output_file)
        tab_obj.run()

    def set_results_filepath(self, fp):
        logger.info("setting result file to : [{}]".format(fp))
        self.results_file = fp
    
    def set_input_filepath(self,fp):
        logger.debug("setting input file to {}".format(str(os.path.abspath(fp))))
        self.input_sequence = os.path.abspath(fp)

    def set_current_db(self, idx, db_params):
        flag_setdb = False
        if db_params and 'dbtype' in db_params:
            self.idx = idx
            self.db_params = db_params
            flag_setdb = True
        else:
            logger.warning("Database {} has no values in configuration file..Skipping..".format(str(dbname)))

        return flag_setdb

    def set_output_dict(self, output_dict):
        self.output_dict = output_dict

    def get_xml_filepath(self,fp):
        return self.results_xml_file

    def get_hmm_filepath(self,fp):
        return self.results_hmm_file

    def get_input_filepath(self):
        return self.input_sequence

    def get_current_db(self):
        return self.idx, self.db_params
 
    def get_output_dict(self):
        return self.output_dict
   
#    def write_stub_output_file(self):
#        # write empty output file if there are no open reading frames
#        with open(os.path.join(self.output_file), 'w') as fout:
#            fout.write(json.dumps({}))
