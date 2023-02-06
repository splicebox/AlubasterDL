# eXAlu
The eXAlu is a Deep Learning model to predict *Alu* Exonization events in the human genome based on sequence alone.

# Installation
The eXAlu works on Linux (tested), Windows and macOS. It requires Python 3.9+, CUDA 11.2+, and PyTorch 1.10+.

We recommend users to install and use this tool in a [conda](https://www.anaconda.com/) envirnoment. Please follow these steps to configure the proper envirnoment.

1. To create a conda envirnoment and activate it,
```
conda create -n exalu python=3.9
conda activate exalu
```
2. Install PyTorch, scikit-learn, tensorboard, matplotlib, pybedtools, seaborn
```
conda install pytorch torchvision cudatoolkit=11.3 -c pytorch
conda install scikit-learn tensorboard matplotlib -c conda-forge
conda install pybedtools -c bioconda
conda install seaborn -c anaconda
```
3. To install eXAlu in developing mode, enter the project root directory, then, 
```
pip install -e .
```

# Usage

## Inference
The program can take bed file or fasta file as input, that user need to specify input mode,
```
python infer.py {bed,fasta} ...

positional arguments:
  {bed,fasta}  select bed or fasta input mode
    bed        infer with bed file
    fasta      infer with fasta file
```

To input a bed file,
```
python infer.py bed -b ALU_BED_FILE -r REF_GENOME_FILE -m MODEL_WEIGHTS_FILE -o OUTPUT_DIR

optional arguments:
  -b ALU_BED_FILE       the input Alu bed file
  -r REF_GENOME_FILE    the reference genome file, it should be hg38c.fa
  -m MODEL_WEIGHTS_FILE the trained model weights file
  -o OUTPUT_DIR         the directory contains temp files and final output file, default ./out
```

To input a fasta file,
```
python infer.py fasta -f ALU_FASTA_FILE -m MODEL_WEIGHTS_FILE -o OUTPUT_DIR

optional arguments:
  -f ALU_FASTA_FILE     the input Alu fa file
  -m MODEL_WEIGHTS_FILE the trained model weights file
  -o OUTPUT_DIR         the directory contains temp files and final output file, default ./out
```

### Example
Below is an example which shows you how to perform inference on a small *Alu* bed file using the trained network weights,

```
conda activate exalu
cd test/inference
python infer.py bed -b example_alu.bed -r ~/data_lflorea1/shared/genomes/hg38/hg38c.fa -m ../models/model_weights.pt -o ./demo_out
python infer.py fasta -f example_alu.fa -m ../models/model_weights.pt -o ./demo_out
```

## Mutagenesis
To show the effect that muatations on a sequence, we developedd the mutatgenesis program. Within an *Alu* sequence and its 25 bp surrounding context regions, we mutate each base three times. We calculate and plot the differences between mutated predicted scores and baseline predicted scores.

```
python mutagenesis.py [-h] {bed,fasta} ...

positional arguments:
  {bed,fasta}  select bed or fasta input mode
    bed        infer with bed file
    fasta      infer with fasta file

optional arguments:
  -h, --help   show this help message and exit
```
The bed file input contains *Alu* sequences only, while the fasta file input contains the contexted *Alu* sequence.

To input a bed file,
```
python mutagenesis.py bed [-h] -b ALU_BED_FILE -r REF_GENOME_FILE -m MODEL_WEIGHTS_FILE -o OUTPUT_DIR [--yaxis Y_AXIS_MODE]

optional arguments:
  -h, --help            show this help message and exit
  -b ALU_BED_FILE       the input Alu bed file
  -r REF_GENOME_FILE    the reference genome file, it should be hg38c.fa
  -m MODEL_WEIGHTS_FILE
                        the trained model weights file
  -o OUTPUT_DIR         the directory contains temp files and final output file
  --yaxis Y_AXIS_MODE   limits of y-axis is fixed to +/-0.3 or adaptive. The default is fixed mode
```
Since we need the formatted description lines (start with ">") to lable the sequences and plot text information inside the output images, you may want to format the description lines in the fasta input file like below,
```
>h38_mk_AluY::chr12:70285190-70285525(-)
```
To input a fasta file,
```
python mutagenesis.py fasta [-h] -f ALU_FASTA_FILE -m MODEL_WEIGHTS_FILE [-o OUTPUT_DIR] [--yaxis Y_AXIS_MODE]

optional arguments:
  -h, --help            show this help message and exit
  -f ALU_FASTA_FILE     the input Alu fa file
  -m MODEL_WEIGHTS_FILE
                        the trained model weights file
  -o OUTPUT_DIR         the directory contains temp files and final output file, default ./out
  --yaxis Y_AXIS_MODE   limits of y-axis is fixed to +/-0.3 or adaptive. The default is fixed mode
```

The output images are located in ./demo_out/imgs, and the text files having all of the changes data are located in ./demo_out/tables.

### Example
Below is an example shows you how to plot the mutagenesis graphs giving bed or fasta file input,
```
conda activate exalu
cd test/analysis/mutagenesis
python mutagenesis.py bed -b ./example_alu.bed -r ~/data_lflorea1/shared/genomes/hg38/hg38c.fa -m ../../models/model_weights.pt -o ./demo_out --yaxis fixed
python mutagenesis.py fasta -f ./example_alu.fa -m ../../models/model_weights.pt -o ./demo_out --yaxis adaptive
```

