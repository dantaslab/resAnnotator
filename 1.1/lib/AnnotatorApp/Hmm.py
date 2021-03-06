from AnnotatorApp.settings import *
from AnnotatorApp.Database import Database

class Hmm(Database):
    """Class to align protein and translated DNA against HMM database using Hmmer3."""
    src_cnt = 0

    def __init__(self, input_file, output_file=None, temp_dir=None, db_params=None, num_threads=4):
        """Creates HMM object for running HMMER algorithm."""
        self.input_file = input_file
        
        if output_file == None:
            f_path, f_name = os.path.split(input_file)
            self.output_file = os.path.join(f_path,"{}.hmmRes_tab.txt".format(f_name))
        else:
            self.output_file = output_file
        
        self.db_params = db_params
        self.dbloc = db_params['dbpath']
        self.dbname = db_params['dbname']
        self.stdout = "/dev/null" #"2> /dev/null"
        self.num_threads = num_threads
        self.outfmt = 'tblout'

        if temp_dir == None:
            self.temp_dir = os.path.split(inputfile)[0]
        else:
            self.temp_dir = temp_dir

    def __repr__(self):
        """Returns Hmm class full object."""
        return "Hmm({}".format(self.__dict__)

    def run(self):
        """Runs Hmmer3""" 
        logger.info('Running HMM against {}'.format(self.dbname))
        logger.debug("run HMM => {}".format(str(self.__dict__)))
 
        db_path = os.path.join(self.dbloc, self.dbname)

        if os.path.exists(db_path + ".h3m") == True and os.path.exists(db_path + ".h3i") == True \
            and os.path.exists(db_path + ".h3f") == True and os.path.exists(db_path + ".h3p") == True:
            
            if os.access(self.input_file, os.R_OK) and is_fasta(self.input_file):
                logger.debug('Aligning input seq {}'.format(self.input_file))
                
                if not os.path.exists(self.temp_dir):
                    os.makedirs(self.temp_dir)

                Hmm.src_cnt += 1

                script_file = os.path.join(self.temp_dir, "s1" + str(Hmm.src_cnt) + "_run_hmm.sh")
                fh1 = open(script_file, 'w')
                extra_params = ''

                if 'gathering_threshold' in self.db_params and self.db_params.get('gathering_threshold', None) is not None:
                    extra_params += '--cut_ga'
                    
                if 'params' in self.db_params and self.db_params.get('params', None) is not None:
                    extra_params += '{} '.format(str(self.db_params.get('params')))

                logger.debug("Extra Params = {}".format(str(extra_params)))
                
                cmd = 'hmmscan --cpu {num_threads} -o {stdout} --{outfmt} {output_file} {params} {hmmdb} {seqfile}' \
                        .format(
                                num_threads = self.num_threads,
                                seqfile = self.input_file, 
                                hmmdb = db_path,
                                output_file = os.path.join(self.temp_dir,self.output_file), 
                                stdout=self.stdout,
                                outfmt = self.outfmt,
                                params = extra_params 
                        )

                logger.debug("Cmd = {}".format(str(cmd)))
                fh1.write(cmd + '\n')
                fh1.close()
                os.system(cmd)
                
                #if os.stat(script_file).st_size > 0:
                #    sbatch_params = '--mem=2000M'
                #    run_sbatch_script(script_file, 1, self.temp_dir, sbatch_params)
                #else:
                #    logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
                #    sys.exit(1)

            else:
                logger.error("input file {} is empty or not in fasta format".format(self.input_file))    
        else:
            logger.error('database index file does not exist for {}'.format(self.dbname))





        #os.system('{program} -query {input} -db {path} \
        #            -outfmt {outfmt} -out {output_file}' \
        #            .format(
        #                program=self.program, 
        #                outfmt=self.outfmt,
        #                input=self.input_file,
        #                path=os.path.join(self.dbloc,self.dbname), 
        #                output_file=self.output_file
        #            )
        #        )
