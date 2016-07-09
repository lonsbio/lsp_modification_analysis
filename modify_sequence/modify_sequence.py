#!/bin/env python

# coding: utf-8

import logging
import time
import datetime
import re
import os
import sys
import protein_variation_functions as seqmod



from subprocess import Popen, PIPE, check_output, STDOUT
import pandas as pd
from io import StringIO
class MorphedSequence:


    def __init__(self,header,description,seq,bootstrap,signal_function):


        self.seq = seq
        self.header = header
        self.description = description
        self.signal_function = signal_function
        self.rev = seqmod.protein_reverse(seq)

        self.random = []
        for _ in range(0,bootstrap):
            self.random.append(seqmod.protein_random(seq))

        drange = signal_function(seq)

        if drange:
            self.positive_data = True
            self.drange = drange

        else:
            self.positive_data = False
            self.drange = '1-30'

        self.sp_cterm = seqmod.protein_sp_cterm(self.seq,self.drange)
        self.sp_remove = seqmod.protein_sp_remove(self.seq,self.drange)
        self.sp_random = []
        for _ in range(0,bootstrap):
            self.sp_random.append(seqmod.protein_sp_random(self.seq,self.drange))





def call_signalP(string_seq):
    """
    Call signalP and return the range of signal peptide predicted - none for no predictions.
    """
    # works for single sequence only
    exec_path = ['/home/ubuntu/signalp4/signalp-4.1/signalp']
    exec_options = [ '-t','euk','-f','summary']
   
    string_seq = ">test\n" + string_seq
   
    seq = str.encode(string_seq)
   
    p=Popen(exec_path + exec_options, stdin=PIPE, stdout=PIPE,stderr=PIPE,shell=True)
    (exec_stdout,exec_stderr) = p.communicate(input=seq)
    exec_stdout 

    stdout_csv="#"
    signalp_headers=''
    for line in exec_stdout.decode().split('\n'):
        stdout_csv = stdout_csv + '\n' + ','.join(line.split())
    result = pd.read_csv(StringIO(stdout_csv), sep=",",comment='#',index_col=0,warn_bad_lines=True,error_bad_lines=True,header=None)
    result.columns = ['Cmax','Cpos','Ymax','Ypos','Smax','Spos','Smean','D','?','Dmaxcut','Networks-used']
    result.index.name="Accession"
   
    sp_found = result["?"] == "Y"


    if sp_found.bool():
        return '1-'+str(result[result["?"] == "Y"]['Ypos'][0]-1)
    else:
        return None



def all(args):
    """
    run all sequence modification changes
    """


    import random as rand
    import quantumrandom as qrand
    import platform

    if (args.seed):
        RANDOM_SEED=args.seed
        logging.info("User supplied random seed")
    else:
        RANDOM_SEED= int(qrand.randint(0, 20000000000000))

    logging.info("Random seed is: %s",RANDOM_SEED)
    rand.seed(RANDOM_SEED)


    logging.info("==Python Details==")
    logging.info("Version: %s",platform.python_version())
    logging.info("Implementation: %s",platform.python_implementation())
    logging.info("Platform: %s",platform.platform())



    import os
    morph_dir=args.output+"/morph"

    if not os.path.exists(morph_dir):
        os.makedirs(morph_dir)
    import skbio

    bootstrap = args.bootstraps

    for f in args.files:
        
        exp_sequences = []
        # since SecretomeP will remove duplciates in input, implement a ful maping from original header to an incrementing ID
        morph_id=1
        logging.info("Morph_ID initialised at: %d", morph_id)

        sk_seqs = skbio.io.read(f, format='fasta')

        for s in sk_seqs:
            print(s.metadata)
            exp_sequences.append(MorphedSequence(s.metadata['id'],s.metadata['description'],str(s),bootstrap,call_signalP))




        output_stub = morph_dir + "/" + os.path.basename(f.name)
        exclude_less_than=41
        with open(output_stub+"_original.fasta", 'w') as original, open(output_stub+"_rev.fasta", 'w') as reverse,open(output_stub+"_spremove.fasta", 'w') as sp_remove,open(output_stub+"_spcterm.fasta", 'w') as sp_cterm,open(output_stub+"_sprandom.fasta", 'w') as sp_random,open(output_stub+"_random.fasta", 'w') as random:
         for m in exp_sequences:
          if (len(m.seq) < exclude_less_than):
            continue
          # use annotations for headers
          # O)riginal with description intact
          # R)everse
          # r(A)ndom with iteration
          # S)p_remove
          # sp_ra(N)dom with iteration
          # sp_(C)term
          # i.e O1, R2,  A3 
          #
          skbio.sequence.Protein(sequence = m.seq , metadata= {'id': "O"+str(morph_id), 'description': m.header +":"+m.description } ).write(original)
          morph_id = morph_id+1 

          skbio.sequence.Protein(sequence = m.rev , metadata= {'id': "R"+str(morph_id), 'description':m.header +":" + 'sequence reversed'} ).write(reverse)
          morph_id = morph_id+1 

          if (m.positive_data):
             skbio.sequence.Protein(sequence = m.sp_remove , metadata= {'id': "S"+str(morph_id), 
                                                                       'description': m.header+": "+'SP at positions' + m.drange + ' removed using '+str(m.signal_function.__name__)} ).write(sp_remove)
             morph_id = morph_id+1 
             skbio.sequence.Protein(sequence = m.sp_cterm , metadata= {'id': "C"+str(morph_id), 
                                                                      'description': m.header +": "+'SP at positions' + m.drange + ' placed at  C-terminus using '+ str(m.signal_function.__name__)} ).write(sp_cterm)
             morph_id = morph_id+1 
          else:
             skbio.sequence.Protein(sequence = m.sp_remove , metadata= {'id': "S"+str(morph_id),
                                                                       'description': m.header+": "+ 'False SP at positions' + m.drange + ' removed (negative data)'} ).write(sp_remove)
             morph_id = morph_id+1 
             skbio.sequence.Protein(sequence = m.sp_cterm , metadata= {'id': "C"+str(morph_id),
                                                                      'description': m.header +": "+'False SP at positions' + m.drange + ' placed at  C-terminus (negativa data)'} ).write(sp_cterm)
             morph_id = morph_id+1 


          randoms = m.random
          iteration=1
          for r in randoms:
            skbio.sequence.Protein(sequence = r , metadata= {'id': "A"+str(morph_id), 'description':  m.header+":"+  ' randomised iteration '+str(iteration)} ).write(random)
            iteration = iteration +1
            morph_id = morph_id+1 

          sp_randoms = m.sp_random
          iteration=1
          for spr in sp_randoms:
                skbio.sequence.Protein(sequence = spr , metadata= {'id': "N"+str(morph_id), 'description': m.header+": "+' SP randomised' +' using '+ str(m.signal_function.__name__) + 'iteration '+str(iteration)} ).write(sp_random)
                iteration = iteration +1
                morph_id = morph_id+1 


          logging.info("Morph_ID finished at: %d", morph_id)



def echo(args):
    """
    echo
    """
    print("No subcommands given. use -h or --help for a list of subcommands")

def info(args):
    """
    Function to print to STDERR to ask for a subcommand, unless doi or citation information is required
    """

    doi_text = "No DOI"

    citation_text = "Unpublished"

    if(args.doi):
        print("{} {}".format(sys.argv[0] , doi_text))
    elif(args.citation):
        print("{} \n {}".format(sys.argv[0],citation_text))
    else:
        print("A subcommand is required. Use -h or --help for more information.")

def main():


    import argparse
    parser = argparse.ArgumentParser()
    logginggroup = parser.add_mutually_exclusive_group()
    logginggroup.add_argument('--debug', action='store_true', help='include debug information in log')
    logginggroup.add_argument('--verbose', action='store_true',help='include extra information in log')
    logginggroup.add_argument('--quiet', action='store_true',help='include only errors in log')

    parser.add_argument('--logfile', type=str, help='log filename')
    parser.add_argument('-v','--version', action='version', version='%(prog)s 1.0')
    parser.set_defaults(func=echo)


    subparsers = parser.add_subparsers(help='subcommands of %(prog)s',                                       title='subcommands')

    #citation information
    parser_info = subparsers.add_parser('info', help='get citation info')
    parser_info.add_argument('--citation', action='store_true', default=True, help='display citation only')
    parser_info.add_argument('--doi', action='store_true',help='display doi only')
    parser_info.set_defaults(func=info)

    #subcommmand 1
    parser_all = subparsers.add_parser('all', help='do tasks required for all')
    parser_all.add_argument('-s','--seed',type=int,help='Integer to use as random seed. A random one will be used otherwise')
    parser_all.add_argument('-o','--output',type=str,help='Directory to store output. Default current directoy',default=".")
    parser_all.add_argument('-b','--bootstraps',type=int,help='Number of bootstraps, default is 500',default=500)
    parser_all.add_argument('files',nargs='*', type=argparse.FileType('r'), default=sys.stdin,help="input files")
    parser_all.set_defaults(func=all)

    args = parser.parse_args()
    print(args)
    if (args.debug):
        loglevel=logging.DEBUG
    elif(args.verbose):
        loglevel=logging.INFO
    elif(args.quiet):
            loglevel=logging.ERROR
    else:
        loglevel=logging.WARNING

    if (args.logfile):
        logging.basicConfig(filename=args.logfile,level=loglevel,
                            format='%(asctime)s %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(filename='{}_{}.log'.format(sys.argv[0],
                            datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d%H%M%S')),level=loglevel,format='%(asctime)s %(message)s',datefmt='%m/%d/%Y %I:%M:%S %p')

    args.func(args)

if __name__ == "__main__":
    main()

