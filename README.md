# resAnnotator.py 
Description: searches amino acid sequence against known ARG specific databases sequentially and assigns annotation using BLAST and HMM approach.

## Usage
      resAnnotator [-h] -i INPUT_SEQUENCE -o OUTPUT_DIR -c CONFIG_FILE
                    [--clean] [-t NUM_THREADS] [-p PREFIX] [-v]

      ResAnnotator - 0.5.0 - The program assigns antibiotic resistance function
      using blast and HMM based databases

      optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_SEQUENCE, --input_sequence INPUT_SEQUENCE
                        input file must be in FASTA format (contig and
                        protein)
      -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        specify output folder
      -c CONFIG_FILE, --config CONFIG_FILE
                        path to the config file
      --clean               removes temporary files
      -t NUM_THREADS, --num_threads NUM_THREADS
                        number of threads to use(default=4)
      -p PREFIX, --prefix PREFIX
                        prefix for the output files
      -v, --version         prints software version number

    A program to annotate antibiotic resistance genes/proteins in contigs/protein
    sequences using BLAST(Uniprot, CARD, Resfinder) and/orHMM (Resfams, PFam-A)
    databases
    
## Config file template: example_config_file.txt
