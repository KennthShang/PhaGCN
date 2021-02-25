import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
import pandas as pd
import subprocess
import argparse
import re

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default = 'contigs.fa')
parser.add_argument('--len', type=int, default=8000)
args = parser.parse_args()


if not os.path.exists("input"):
    _ = os.makedirs("input")
else:
    print("folder {0} exist... cleaning dictionary".format("input"))
    if os.listdir("input"):
        try:
            _ = subprocess.check_call("rm -rf {0}".format("input"), shell=True)
            _ = os.makedirs("input")
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)

if not os.path.exists("pred"):
    _ = os.makedirs("pred")
else:
    print("folder {0} exist... cleaning dictionary".format("pred"))
    if os.listdir("pred"):
        try:
            _ = subprocess.check_call("rm -rf {0}".format("pred"), shell=True)
            _ = os.makedirs("pred")
            print("Dictionary cleaned")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)



if not os.path.exists("Split_files"):
    _ = os.makedirs("Split_files")
else:
    print("folder {0} exist... cleaning dictionary".format("Split_files"))
    if os.listdir("Split_files"):
        try:
            _ = subprocess.check_call("rm -rf {0}".format("Split_files"), shell=True)
            _ = os.makedirs("Split_files")
        except:
            print("Cannot clean your folder... permission denied")
            exit(1)
try:
    make_diamond_cmd = 'diamond makedb --threads 8 --in database/Caudovirales_protein.fasta -d database/database.dmnd'
    print("Creating Diamond database...")
    res = subprocess.run(make_diamond_cmd, check=True, stdout=subprocess.PIPE)
    if res.returncode != 0:
        print('Error creating Diamond database')
        exit(1)
    diamond_cmd = 'diamond blastp --threads 8 --sensitive -d database/database.dmnd -q database/Caudovirales_protein.fasta -o database/database.self-diamond.tab'
    print("Running Diamond...")
    _ = subprocess.check_call(diamond_cmd, shell=True)
    diamond_out_fp = "database/database.self-diamond.tab"
    database_abc_fp = "database/database.self-diamond.tab.abc"
    _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
except:
    print("create database failed")
    exit(1)


#####################################################################
##########################    Start Program  ########################
#####################################################################

def special_match(strg, search=re.compile(r'[^ACGT]').search):
    return not bool(search(strg))


cnt = 0
file_id = 0
records = []
for record in SeqIO.parse(args.contigs, 'fasta'):
    if cnt !=0 and cnt%1000 == 0:
        SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta") 
        records = []
        file_id+=1
        cnt = 0
    seq = str(record.seq)
    seq = seq.upper()
    if special_match(seq):
        if len(record.seq) > args.len:
            records.append(record)
            cnt+=1

SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta")
file_id+=1

for i in range(file_id):
    cmd = "mv Split_files/contig_"+str(i)+".fasta input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Moving file Error for file {0}".format("contig_"+str(i)))
        continue

    cmd = "python run_CNN.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Pre-trained CNN Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue
        

    cmd = "python run_KnowledgeGraph.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Knowledge Graph Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    cmd = "python run_GCN.py"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("GCN Error for file {0}".format("contig_"+str(i)))
        cmd = "rm input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    # Clean files
    cmd = "rm input/*"
    out = subprocess.check_call(cmd, shell=True)

    name_list = pd.read_csv("name_list.csv")
    prediction = pd.read_csv("prediction.csv")
    prediction = prediction.rename(columns={'contig_names':'idx'})
    contig_to_pred = pd.merge(name_list, prediction, on='idx')
    contig_to_pred.to_csv("pred/contig_"+str(i)+".csv", index = None)

    cmd = "rm name_list.csv prediction.csv"
    out = subprocess.check_call(cmd, shell=True)


cmd = "cat pred/* > final_prediction.csv"
out = subprocess.check_call(cmd, shell=True)

