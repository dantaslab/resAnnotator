#!/usr/bin/env python3

from AnnotatorApp.settings import *
from AnnotatorApp.Annotator import Annotate

import argparse
import yaml

class MainBase(object):
    def __init__(self, api=False):

        parser = self.main_args()
        self.args = vars(parser.parse_args(sys.argv[1:]))
        
        #if not os.path.exists(self.args['output_dir']):
        #    os.makedirs(self.args['output_dir'])

        if not os.path.isdir(self.args['output_dir']):
            try:
                os.makedirs(self.args['output_dir'])
            except EnvironmentError:
                sys.exit("CRITICAL ERROR: Unable to create the directory")

        if not os.access(self.args['output_dir'], os.W_OK):
            sys.exit("CRITICAL ERROR: The output directory is not " +
                    "writeable. This software needs to write files to this directory.\n" +
                    "Please select another directory.")

        f_name = os.path.splitext(os.path.basename(self.args['input_sequence']))[0]     
        log_filename = 'zz_{}.log'.format(str(os.path.basename(f_name)))

        self.logger = logger or logging.getLogger(__name__)
        self.file_handler = logging.FileHandler(os.path.join(self.args['output_dir'],log_filename))
        self.file_handler.setFormatter(formatter)
        self.logger.addHandler(self.file_handler)

        self.main_run()     

    def main_args(self):
        parser = argparse.ArgumentParser(prog="resAnnotator", description="{} - {} - The program assigns antibiotic resistance function using blast and HMM based databases".format(APP_NAME,SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY)
        parser.add_argument(
            '-i','--input_sequence', 
            dest = "input_sequence", 
            required = True,
            help = 'input file must be in either FASTA (contig and protein)')
        parser.add_argument(
            '-o','--output_dir', 
            dest = "output_dir", 
            required = True, 
            help = "specify output folder")
        parser.add_argument(
            '-c','--config', 
            dest = "config_file", 
            required = True,
            help = "path to the config file")
        parser.add_argument(
            '--clean', 
            dest = "clean", 
            action = "store_true", 
            help = "removes temporary files")
        parser.add_argument(
            '-t', '--num_threads',
            dest = "num_threads",
            default = 4,
            help = "number of threads to use(default=4)")
        parser.add_argument(
            '-v','--version',
            action = 'version', 
            version = "{}".format(SOFTWARE_VERSION), 
            help = "prints software version number")
        
        return parser

    def main_run(self):
        print("\nStarted the ResAnnotator Pipeline v{}".format(str(SOFTWARE_VERSION)))
        #print("Before parsing =>\n {}".format(str(self.__dict__)))
        configuration = yaml.load(open(self.args['config_file']).read())
        
        if configuration['log_level'] in ['DEBUG','INFO','WARNING','ERROR','CRITICAL']:
            self.logger.setLevel(configuration['log_level'])
        else:
            print("Error Occurred!..Invalid logger level value")
            self.logger.setLevel(level)
            sys.exit(1)

        #logger.debug("parameters in config file =>\n{}".format(str(configuration)))
        args = self.args.update(configuration)
        res_obj = Annotate(**self.args)
        #print("After parsing =>\n{}".format(str(res_obj.__dict__)))
        res_obj.run()


if __name__ == '__main__':
    MainBase()
