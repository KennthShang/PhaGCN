# PhaGCN

PhaGCN is a GCN based model, which can learn the species masking feature via deep learning classifier, for new Phage taxonomy classification. To use PhaGCN, you only need to input your contigs to the program.


# Required Dependencies
* Python 3.x
* Numpy
* Pytorch
* cuda 10.1 (only for gpu accelerate mode)
* Networkx
* Pandas

## An easiler way to install
We recommend you to install all the package with [anaconda](https://anaconda.org/)

After cloning this respository, you can use anaconda to install the **environment.yaml**. This will install all packages you need with gpu mode.

# Usage (example)
Here we present an example to show how to run PhaGCN. We support a file named "contigs.fa" in the Github folder and it contain contigs simulated from E. coli phage. The only command that you need to run is `python run_Speed_up.py --contigs contigs.fa -len 8000`. 

There are two parameters for the program: 1. `--contigs` is the path of your contigs file. 2. `--len` is the length of the contigs you want to predict. As shown in our paper, with the length of contigs increases, the recall and precision increase. We recommend you to choose a proper length according to your needs. The default length is 8000bp. 

# References
how to cite this tool:
Jiayu Shang, Jingzhe Jiang and Yanni Sun, Bacteriophage classification for assembled contigs using Graph Convolutional Network, submitted to ISMB 2021
