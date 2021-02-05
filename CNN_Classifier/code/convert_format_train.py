import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def create_reads(file_name):
    with open("stride50_train/"+file_name, 'w') as file:
        for record in SeqIO.parse("filtered_train/"+file_name, "fasta"):
            file.write(str(record.seq) +'\n')
        file.close()

if __name__ == "__main__":
    path = "filtered_train/"
    name_list = os.listdir(path)
    for name in name_list:
        create_reads(name)
        #print(name + " finished")    