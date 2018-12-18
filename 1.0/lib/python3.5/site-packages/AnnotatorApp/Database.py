from AnnotatorApp.settings import *
import shutil

#from AnnotatorApp.Annotator import Annotate

class Database(object):
    """Class to create BLAST databases"""
    seqtype = 'prot'
    src_cnt = 0
    
    def __init__(self, db_params, temp_dir, aligner='BLAST'):
        """Creates Database object."""
        self.db_params = db_params
        self.stdout = "2>&1 >> /dev/null" #"2> /dev/null"
        self.aligner = aligner
        self.temp_dir = temp_dir

    def __repr__(self):
        """Returns Database class full object."""
        return "Database({}".format(self.__dict__)

    def build_databases(self):
        """Build databases."""
        if self.db_params['dbtype'].upper() == 'BLAST':
            if self.aligner == 'diamond' and self.db_params['seqtype'].lower()  == 'prot':
                self.make_diamond_database()
            else:
                self.make_blastdb()
        elif self.db_params['dbtype'].upper() == 'HMM':
            self.make_hmmdb()
        else:
            logger.exception("Invalid database type specified!..should be either BLAST|HMM")
            sys.exit(1)

        #logger.debug("building databases = {}".format(str(self.__dict__)))
        #for (name, params) in self.db.items():
        #    logger.debug("{} => {}".format(name, str(params)))
        #    if params['dbtype'] == 'BLAST':
        #        if self.aligner == 'diamond':
        #            self.make_diamond_database(name, params)
        #        else:
        #            self.make_blastdb(name, params)
        #    elif params['dbtype'] == 'HMM':
        #        self.make_hmmdb(name, params)
        #    else:
        #        logger.exception("Invalid database type specified!")
        #        sys.exit(1)

    def make_blastdb(self):
        """Build blast database from fasta file"""
        index_params = self.db_params
        db_path = os.path.join(index_params['dbpath'])
        db_name = os.path.basename(index_params['dbname'])
        db_seqtype = index_params['seqtype']
        db_status = self.check_db_exists(db_path, db_seqtype)
        
        if db_status:
            #os.path.exists(os.path.join(db_path, db_name + ".phr")) == True \
            #and os.path.exists(os.path.join(db_path, db_name + ".pin")) == True \
            #and os.path.exists(os.path.join(db_path, db_name + ".psq")) == True:
            logger.info("Index files exist for database {}".format(db_name))
            pass
        else:
            logger.warning("Index files does not exist for database {}".format(db_name))
            fastafile = os.path.abspath(os.path.join(index_params['dbpath'],index_params['dbname']))              
            logger.debug("Checking if fasta file {} exists".format(fastafile))
            logger.debug("before setting dbpath path variable {}".format(str(self.db_params)))
            
            if os.access(fastafile, os.R_OK) and is_fasta(fastafile):
                logger.info("Creating Blast database using fasta file")
                self.set_db_indexpath(db_name, os.path.abspath(self.temp_dir))
                Database.src_cnt += 1

                if not os.path.exists(self.temp_dir):
                    os.makedirs(self.temp_dir)
                
                script_file = os.path.join(self.temp_dir, "s0" + str(Database.src_cnt) + "_makeblastdb.sh")
                fh1 = open(script_file, 'w')
                cmd = 'makeblastdb -in {input_file} -dbtype {db_type} -out {output_dir} {stdout}' \
                        .format(
                                input_file = fastafile, 
                                db_type = index_params['seqtype'], 
                                output_dir = os.path.join(os.path.abspath(self.temp_dir), db_name), 
                                stdout=self.stdout
                        )

                logger.debug("cmd = {}".format(str(cmd)))
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
                logger.error("Either fasta file does not exists or is not in fasta format. Please check if filename ends with .fasta")

    def make_diamond_database(self):
        """Build DIAMOND database from a FASTA file."""
        index_params = self.db_params
        db_path = index_params['dbpath']
        db_name = os.path.basename(index_params['dbname'])
        db_file = os.path.join(db_path, db_name + '.dmnd')

        logger.info("db_path for {} = {}".format(str(db_name), str(db_path)))
        
        if os.path.exists(db_file) == True:
            logger.info("diamond DB exists")
            pass
        else:
            logger.info('Index files do not exist for {}. So, we will now try to create using diamond'.format(db_name))
            fastafile = os.path.abspath(os.path.join(index_params['dbpath'],index_params['dbname']))              
            logger.info("Checking if fasta file {} exists".format(fastafile))

            if os.access(fastafile, os.R_OK) and is_fasta(fastafile):
                logger.info("Creating Diamond database using fasta file")

                self.set_db_indexpath(db_name, os.path.abspath(self.temp_dir))
                Database.src_cnt += 1

                if not os.path.exists(self.temp_dir):
                    os.makedirs(self.temp_dir)

                script_file = os.path.join(self.temp_dir, "s0" + str(Database.src_cnt) + "_makediamonddb.sh")
                fh1 = open(script_file, 'w')
                cmd = 'diamond makedb --in {input_file} --db {db_path}' \
                            .format(
                                    input_file = fastafile, 
                                    db_path = os.path.join(os.path.abspath(self.temp_dir), db_name) 
                                    #stdout = self.stdout
                            )
                
                logger.debug("cmd = {}".format(str(cmd)))
                os.system(cmd)
                fh1.write(cmd + '\n')
                fh1.close()

                #if os.stat(script_file).st_size > 0:
                #    sbatch_params = '--mem=2000M'
                #    run_sbatch_script(script_file, 1, self.temp_dir, sbatch_params)
                #else:
                #    logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
                #    sys.exit(1)
            else:
                logger.error("Either fasta file does not exists or is not in fasta format. Please check if filename ends with .fasta")
            
            #os.system('diamond makedb --quiet --in {} --db {} {stdout}'.format(os.path.join(self.db,"proteindb.fsa"),os.path.join(self.db,"protein.db"),stdout=self.stdout))

    def make_hmmdb(self):
        """Build HMM database if not exists"""
        logger.debug('HMM database')
        index_params = self.db_params
        db_path = index_params['dbpath']
        db_name = index_params['dbname']
        db_file = os.path.join(index_params['dbpath'], index_params['dbname'])
        
        if os.path.exists(db_file + ".h3m") == True \
            and os.path.exists(db_file + ".h3i") == True \
            and os.path.exists(db_file + ".h3f") == True \
            and os.path.exists(db_file + ".h3p") == True:
            logger.info("Index files exist for database {}".format(db_name))
            pass
        else:
            logger.warning("Index files are not present {}".format(db_name))
            logger.info("Checking if hmm file {} exists".format(os.path.abspath(db_file)))
            if os.access(os.path.abspath(db_file), os.R_OK) and is_hmm(os.path.abspath(db_file)):
                logger.info("Creating HMM index database using hmm file")
                hmmfile = os.path.join(self.temp_dir, index_params['dbname'])
                if not os.path.exists(hmmfile):
                    shutil.copyfile(os.path.abspath(db_file), os.path.abspath(hmmfile))
       
                self.set_db_indexpath(db_name, os.path.abspath(self.temp_dir))
                Database.src_cnt += 1
            
                if not os.path.exists(self.temp_dir):
                    os.makedirs(self.temp_dir)
                
                script_file = os.path.join(self.temp_dir, "s0" + str(Database.src_cnt) + "_hmmpress.sh")
                #print(script_file)
                
                fh1 = open(script_file, 'w')
                cmd = 'hmmpress {}'.format(hmmfile)
                logger.debug("cmd = {}".format(str(cmd)))
                os.system(cmd)
                fh1.write(cmd + '\n')   
                fh1.close()

                #if os.stat(script_file).st_size > 0:
                #    sbatch_params = '--mem=2000M'
                #    run_sbatch_script(script_file, 1, self.temp_dir, sbatch_params)
                #else:
                #    logger.error("The script file {} is empty. Exiting....".format(script_file), exc_info=True)
                #    sys.exit(1)
            else:
                logger.error("Either fasta file does not exists or is not in fasta format. Please check if filename ends with .fasta")
    

    def set_db_indexpath(self, dbname, db_dir):
        """Sets blast xml filepath."""
        logger.info("setting path for {}".format(dbname))
        self.db_params['dbpath'] = os.path.abspath(db_dir)
        #self.db_params[dbname]['dbpath'] = db_dir
        logger.debug("After setting database path = {}".format(str(self.db_params['dbpath'])))


    def check_db_exists(self, db_path, db_seqtype):
        if db_seqtype == 'prot':
            if os.path.exists(db_path + ".phr") == True and os.path.exists(db_path + ".pin") == True \
            and os.path.exists(db_path + ".psq") == True:
                return True
            else:
                logger.warning("Index files does not exists")
                return False
        elif db_seqtype == 'nucl':
            if os.path.exists(db_path + ".nhr") == True and os.path.exists(db_path + ".nin") == True \
            and os.path.exists(db_path + ".nsq") == True:
                return True
            else:
                logger.warning("Index files does not exists")
                return False
        else:
            logger.error('Invalid database seq type!')
            return False

