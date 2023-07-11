# eXAlu
A Deep Learning model to predict *Alu* exonization events in the human genome based on sequence alone.

```
Copyright (C) 2022-2023, and GNU GPL v3.0, by Zitong He, Liliana Florea
```
This software is beta version.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## What is eXAlu
*Alu* repeats make up ~11% of the human genome, with more than one million copies. While most *Alu* elements map to nonfunctional intergenic and intronic sequences, some become incorporated into genes. *Alu* exonization, in which an intronic *Alu* sequence is recruited into a gene transcript (see below), can disrupt gene function and potentially lead to disease. *Alu* exons are estimated to represent 5% of alternatively spliced exons in the human genome, but the extent of *Alu* exonization in humans is not yet known.

![transcript_alu](images/transcript_alu.png)

We developed a Deep Learning (DL) model eXAlu that predicts *Alu* exonization events from genomic sequences. The model takes as input *Alu* sequences surrounded by 25 bp context sequences, and the output is the probability that the *Alu* will cause exonization. We provide model weights trained on the human *Alu* sequences. The training dataset consists of *Alu* exonization events in 28 human tissues, extracted from RNA-seq data generated by the GTEx project. We used the [STAR] (https://github.com/alexdobin/STAR) and [CLASS2] (https://sourceforge.net/p/splicebox/wiki/CLASS/) tools to map and assemble reads into transcripts, then extracted internal exons overlapping *Alu" annotations in the antisense to the gene.

The model's network has four convolutional layers, and corresponding batch-norm layers and pooling layers, followed by fully-connected layers that map the features extracted by convolutional layers to the output probabilities.

<!-- ![model_network](images/model_network.png) -->

This repository provides *inference* and *mutagenesis plotting* functions. The mutagenesis graphs show the difference in the model's score when mutating single nucleotides at each position on given *Alu* sequences.

## Installation
eXAlu works on Linux (tested), Windows and macOS. It requires Python 3.9+, CUDA 11.2+, and PyTorch 1.10+.

We recommend users to install and use this tool in a [conda](https://www.anaconda.com/) envirnoment. Please follow these steps to configure the proper envirnoment.

1. To create a conda environment and activate it,
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

## Usage

### Inference
The program can take as input a bed file, containing genomic *Alu* intervals, or a fasta file, containing the *Alu* plus context sequences, and the user needs to specify the input mode,
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
  -b ALU_BED_FILE       the input *Alu* bed file
  -r REF_GENOME_FILE    the path to the reference genome file (which should be hg38c.fa)
  -m MODEL_WEIGHTS_FILE the trained model weights file
  -o OUTPUT_DIR         the directory containing temp files and final output file, default ./out
```

To input a fasta file,
```
python infer.py fasta -f ALU_FASTA_FILE -m MODEL_WEIGHTS_FILE -o OUTPUT_DIR

optional arguments:
  -f ALU_FASTA_FILE     the input *Alu* fasta file
  -m MODEL_WEIGHTS_FILE the trained model weights file
  -o OUTPUT_DIR         the directory containing temp files and final output file, default ./out
```

#### Example
Below is an example showing how to perform inference on a small *Alu* bed file using the trained network weights,

```
conda activate exalu
cd test/inference
python infer.py bed -b example_alu.bed -r ~/data_lflorea1/shared/genomes/hg38/hg38c.fa -m ../models/model_weights.pt -o ./demo_out
python infer.py fasta -f example_alu.fa -m ../models/model_weights.pt -o ./demo_out
```

### Mutagenesis
To show the effects that sequence mutations have on the model's prediction, we developed a mutagenesis plotting program. Within an *Alu* sequence and its 25 bp surrounding context regions, it mutates each base into each of the three alternate bases, then plots the difference in scores between the mutated and original sequencxes. Lastly, positive and negative peaks are determined with a sliding window algorithm. 

```
python mutagenesis.py [-h] {bed,fasta} ...

positional arguments:
  {bed,fasta}  select bed or fasta input mode
    bed        infer with bed file
    fasta      infer with fasta file

optional arguments:
  -h, --help   show this help message and exit
```
Note that the bed file input contains *Alu* sequences only, while the fasta file input contains the *Alu* sequence plus context.

To input a bed file,
```
python mutagenesis.py bed [-h] -b ALU_BED_FILE -r REF_GENOME_FILE -m MODEL_WEIGHTS_FILE -o OUTPUT_DIR [--yaxis Y_AXIS_MODE]

optional arguments:
  -h, --help            show this help message and exit
  -b ALU_BED_FILE       the input *Alu* bed file
  -r REF_GENOME_FILE    the reference genome file, it should be hg38c.fa
  -m MODEL_WEIGHTS_FILE
                        the trained model weights file
  -o OUTPUT_DIR         the directory contains temp files and final output file
  --yaxis Y_AXIS_MODE   limits of y-axis is fixed to +/-0.3 or adaptive. The default is fixed mode
```
Since we need the formatted description lines (start with ">") to labee the sequences and plot text information within the output images, you may want to format the description lines in the fasta input file like below,
```
>h38_mk_AluY::chr12:70285190-70285525(-)
```
To input a fasta file,
```
python mutagenesis.py fasta [-h] -f ALU_FASTA_FILE -m MODEL_WEIGHTS_FILE [-o OUTPUT_DIR] [--yaxis Y_AXIS_MODE]

optional arguments:
  -h, --help            show this help message and exit
  -f ALU_FASTA_FILE     the input *Alu* fa file
  -m MODEL_WEIGHTS_FILE
                        the trained model weights file
  -o OUTPUT_DIR         the directory contains temp files and final output file, default ./out
  --yaxis Y_AXIS_MODE   limits of y-axis is fixed to +/-0.3 or adaptive. The default is fixed mode
```

The output images are located in ./demo_out/imgs, and the text files with the score change data are located in ./demo_out/tables.

#### Example
Below is an example that shows how to plot the mutagenesis graphs giving bed or fasta file input,
```
conda activate exalu
cd test/analysis/mutagenesis
python mutagenesis.py bed -b ./example_alu.bed -r ~/data_lflorea1/shared/genomes/hg38/hg38c.fa -m ../../models/model_weights.pt -o ./demo_out --yaxis fixed
python mutagenesis.py fasta -f ./example_alu.fa -m ../../models/model_weights.pt -o ./demo_out --yaxis adaptive
```
## Support
Contact: Zitong He, hezt@jhu.edu

## License information
See the file LICENSE for information on the history of this software, terms & conditions for usage, and a DISCLAIMER OF ALL WARRANTIES.
