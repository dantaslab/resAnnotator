from AnnotatorApp.settings import *
from AnnotatorApp.Database import Database

class Blast(Database):
    """Class to create Blast object and align for protein and translated DNA searches."""
    src_cnt = 0

    def __init__(self, input_file, input_type, output_file=None, temp_dir=None, db_params=None, num_threads=4):
        """Creates Blast object for running NCBI BLAST algorithm."""
        self.input_file = input_file
        
        if output_file == None:
            f_path, f_name = os.path.split(input_file)
            self.output_file = os.path.join(f_path,"{}.blastRes.xml".format(f_name))
        else:
            self.output_file = output_file
        
        self.db_params = db_params
        self.dbloc = db_params['dbpath']
        self.dbname = db_params['dbname']
        self.dbseqtype = db_params.get('seqtype', None)
        #self.program = program
        self.stdout = "2>&1 >> /dev/null" #"2> /dev/null"
        self.num_threads = num_threads
        self.input_type = input_type
        self.outfmt = 5
        
        if temp_dir == None:
            self.temp_dir = os.path.split(inputfile)[0]
        else:
            self.temp_dir = temp_dir
    
        if self.input_type.lower() == 'protein' and self.dbseqtype.lower() == 'prot':
            self.program = 'blastp'
        elif self.input_type.lower() == 'protein' and self.dbseqtype.lower() == 'nucl':
            self.program = 'tblastn'
        elif self.input_type.lower() == 'contig' and self.dbseqtype.lower() == 'prot':
            self.program = 'blastx'
        elif self.input_type.lower() == 'contig' and self.dbseqtype.lower() == 'nucl':
            self.program = 'blastn'
        else:
            logger.error("Invalid input seqtype or database seqtype!")
            sys.exit(1)

    def __repr__(self):
        """Returns Blast class full object."""
        return "Blast({}".format(self.__dict__)

    def run(self):
        """Runs BLAST algorithm.""" 
        logger.info('Running blast against {}'.format(self.dbname))
        db_path = os.path.join(self.dbloc, self.dbname)
        db_status = self.check_db_exists(db_path, self.dbseqtype)

        if db_status:
            #print("Blast object values =>\n {}".format(self.__dict__))
            if os.access(self.input_file, os.R_OK) and is_fasta(self.input_file):
                logger.debug('Aligning input seq {}'.format(self.input_file))
                
                if not os.path.exists(self.temp_dir):
                    os.makedirs(self.temp_dir)

                Blast.src_cnt += 1

                script_file = os.path.join(self.temp_dir, "s1" + str(Blast.src_cnt) + "_run_blast.sh")
                fh1 = open(script_file, 'w')
                cmd = "{program} -query {input_file} -db {dbpath} -out {output_file} -outfmt {outfmt}"\
                        .format(
                                program = self.program,
                                input_file = self.input_file,
                                dbpath = db_path,
                                output_file = os.path.join(self.temp_dir, self.output_file),
                                outfmt = self.outfmt
                        )

#                cmd = "time cat {input_file}| parallel --gnu --plain --progress -j+0 --recstart '>' --pipe \
#                        {program} -query - -db {dbpath} -out {output_file} -outfmt {outfmt}" \
#                        .format(
#                                program = self.program, 
#                                input_file = self.input_file, 
#                                dbpath = db_path, 
#                                output_file = os.path.join(self.temp_dir, self.output_file), 
#                                outfmt = self.outfmt 
#                                #stdout=self.stdout
#                        )
                
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

#time cat D-e-a_filtered_old.faa | parallel --gnu --plain --progress -j+0 --recstart '>' --pipe blastp -query - -db _db/card_protein_homolog -outfmt 5 -out parallel_blast_out.xml



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

