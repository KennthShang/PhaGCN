import numpy as np
import os 

k_list = ["A", "C", "G", "T"]
nucl_list = ["A", "C", "G", "T"]


for i in range(2):
    tmp = []
    for item in nucl_list:
        for nucl in k_list:
            tmp.append(nucl+item)
    k_list = tmp


int_to_vocab = {ii: word for ii, word in enumerate(k_list)}
vocab_to_int = {word: ii for ii, word in int_to_vocab.items()}

def encode(file_name):
    file = open("stride50_train/"+file_name) 
    data = file.readlines() 
    feature = []
    for read in data:
        read = read[:-1]
        int_read = []
        for i in range(len(read)):
            if i + 3 > len(read):
                break
                
            int_read.append(vocab_to_int[read[i:i+3]])
        
        
        if len(int_read) != 1998:
            print("error length")
        
        feature.append(int_read)
    name = file_name.split(".")[0]
    np.savetxt("int_train/"+name+".csv", feature, delimiter=",", fmt='%d')


if __name__ == "__main__":
    Load_path = "stride50_train/"
    name_list = os.listdir(Load_path)
    for name in name_list:
        encode(name)
        #print(name + " finished")
