from AnnotatorApp.settings import *
from AnnotatorApp.Database import Database

class Diamond(Database):
    """Class to create Diamond object and align for protein and translated DNA searches."""
    src_cnt = 0

    def __init__(self, input_file, input_type, output_file=None, temp_dir=None, db_params=None, num_threads=4):
        """Creates Diamond object for running DIAMOND algorithm."""
        self.input_file = input_file
        
        if output_file == None:
            f_path, f_name = os.path.split(input_file)
            self.output_file = os.path.join(f_path,"{}.blastRes.xml".format(f_name))
        else:
            self.output_file = output_file
        
        self.db_params = db_params
        self.dbloc = db_params['dbpath']
        self.dbname = db_params['dbname']
        self.stdout = "2>&1 >> /dev/null" #"2> /dev/null"
        self.num_threads = num_threads
        self.outfmt = 5 # xml format
        #self.outfmt = "6 qseqid sseqid qlen slen qstart qend sstart send evalue pident ppos bitscore length nident positive mismatch"
        self.seqtype = db_params.get('seqtype','prot')
        self.input_type = input_type.lower()
        
        if self.input_type.lower() == 'protein':
            self.program = 'blastp'
        elif self.input_type.lower() == 'contig':
            self.program = 'blastx'
        else:
            logger.error('Invalid paramater!!')
            sys.exit()

        if temp_dir == None:
            self.temp_dir = os.path.split(inputfile)[0]
        else:
            self.temp_dir = temp_dir
        
    def __repr__(self):
        """Returns Diamond class full object."""
        return "Diamond({})".format(self.__dict__)

    def run(self):
        """Runs DIAMOND algorithm."""
        logger.debug("Diamond => ({})".format(str(self.__dict__)))
        logger.info('Running Diamond against {}'.format(self.dbname))
        db_path = os.path.join(self.dbloc, self.dbname)

        if os.path.exists(db_path + ".dmnd") == True:
            if os.access(self.input_file, os.R_OK) and is_fasta(self.input_file):
                logger.debug('Aligning input seq {}'.format(self.input_file))

                if not os.path.exists(self.temp_dir):
                    os.makedirs(self.temp_dir)

                Diamond.src_cnt += 1

                script_file = os.path.join(self.temp_dir, "s1" + str(Diamond.src_cnt) + "_run_diamond.sh")
                fh1 = open(script_file, 'w')

                extra_params = ''

                for key, val in self.db_params.items():
                    logger.debug("Key = {} and Value = {}".format(str(key),str(val)))
                    if key == 'query_cov' and 0 < val <= 100:
                        extra_params += '--query-cover {} '.format(str(val))
                    elif key == 'subject_cov' and 0 < val <= 100:
                        extra_params += '--subject-cover {} '.format(str(val))
                    elif key == 'identity_cutoff' and 0 < val <= 100:
                        extra_params += '--id {} '.format(str(val))
                    elif key == 'evalue' and 0 <= float(val) <= 0.001:
                        extra_params += '--evalue {} '.format(str(val))
                    elif key == 'max_target_seqs' and isinstance(val, int) == True:
                        extra_params += '--max-target-seqs {} '.format(str(val))
                    elif key == 'params':
                        extra_params += ' {} '.format(str(val))
                
                logger.debug("Extra Params = {}".format(str(extra_params)))
                
                cmd = 'diamond {program} --in {in_ref} --db {db} --query {input} --out {output_file} --outfmt {outfmt} {params} --threads {num_threads} {stdout}' \
                        .format(
                                program = self.program, 
                                in_ref = db_path + ".fasta",
                                db = db_path + ".dmnd",
                                input = self.input_file, 
                                output_file = os.path.join(self.temp_dir,self.output_file), 
                                outfmt = self.outfmt,
                                params = extra_params, 
                                num_threads = self.num_threads,
                                stdout = self.stdout
                        )

                logger.debug("Cmd = {}".format(str(cmd)))
                os.system(cmd)
                fh1.write(cmd + '\n')
                fh1.close()

                #if os.stat(script_file).st_size > 0:
                #    sbatch_params = '--mem=8000M'
                #    run_sbatch_script(script_file, 1, self.temp_dir, sbatch_params)
                #else:
                #    logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
                #    sys.exit(1)

            else:
                logger.error("input file {} is empty or not in fasta format".format(self.input_file))

        else:
            logger.error('database index file does not exist for {}'.format(self.dbname))



#       os.system('diamond {program} --in {in_ref} --db {db} \
#                  --query {input} --outfmt {outfmt} --out {output_file}  \
#                  --threads {num_threads}  --index-chunks {index_chunks} \
#                  --block-size {block_size}  \
#                  --salltitles  --quiet --more-sensitive' \
#                   .format(
#                       program=self.program,
#                       in_ref=os.path.join(self.db,"proteindb.fsa"),
#                       db=os.path.join(self.dbloc,self.dbname),
#                       input=self.input_file,
#                       output_file=self.output_file,
#                       #num_threads=self.num_threads, 
#                       #index_chunks=self.index_chunks,
#                       #block_size=self.block_size,
#                       outfmt=self.outfmt
#                   )
#               )
