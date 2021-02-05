import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
import pandas as pd
import subprocess


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

#try:
#    _ = subprocess.check_call("mkdir {0}".format("input"), shell=True)
#except:
#    print("folder {0} exist... cleaning dictionary".format("input"))
#    try:
#        _ = subprocess.check_call("rm -rf {0}".format("input"), shell=True)
#        _ = subprocess.check_call("mkdir {0}".format("input"), shell=True)
#    except:
#        print("Cannot clean your folder... permission denied")
#        exit(1)

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

#try:
#    _ = subprocess.check_call("mkdir {0}".format("Split_files"), shell=True)
#except:
#    print("folder {0} exist... cleaning dictionary".format("Split_files"))
#    try:
#        _ = subprocess.check_call("rm -rf {0}".format("Split_files"), shell=True)
#        _ = subprocess.check_call("mkdir {0}".format("Split_files"), shell=True)
#    except:
#        print("Cannot clean your folder... permission denied")
#        exit(1)

cnt = 0
file_id = 0
records = []
for record in SeqIO.parse("contigs.fa", 'fasta'):
   if cnt !=0 and cnt%1000 == 0:
       SeqIO.write(records, "Split_files/contig_"+str(file_id)+".fasta","fasta")
       records = []
       file_id+=1
   if "N" not in record.seq or "n" not in record.seq:
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

