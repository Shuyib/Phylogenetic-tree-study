[![Phylogenetic tree study Continuous Integration/Continuos Delivery](https://github.com/Shuyib/Phylogenetic-tree-study/actions/workflows/devops.yml/badge.svg)](https://github.com/Shuyib/Phylogenetic-tree-study/actions/workflows/devops.yml)
# Phylogenetic-tree-study

Estimating Phylogenetic trees using 30 microorganisms (previously 6 organisms: review [data folder](https://github.com/Shuyib/Phylogenetic-tree-study/tree/master/data) and [notebook](https://github.com/Shuyib/Phylogenetic-tree-study/blob/master/Phylogenetic_trees_unsupervised_learning.ipynb). Looking at the 16S rRNA gene with Unsupervised Learning, web based tools and Molecular Evolutionary Genetics Analysis MEGA7. Further we are looking at motifs and finding out what they do.

It is important to know these regions since they can potentially give use clues about the regions we can target for targeted DNA studies.

Setting up your environment
---
* Download Anaconda for your operating system for Python 3 [anaconda](https://www.anaconda.com/download/)
* Create a conda environment like mine:

  `conda env create -f environment.yml`

  This creates an environment called py35. Activate it with this command in your terminal

  `source activate py35`

* In your terminal, in the directory where you cloned this repository. Run this command

  `jupyter notebook Phylogenetic_trees_unsupervised_learning.ipynb`

I've not provided a codebook/data description file since one of the headings cover that in the notebook.
Otherwise, you can check out the [notebook](https://nbviewer.jupyter.org/github/Shuyib/Phylogenetic-tree-study/blob/master/Phylogenetic_trees_unsupervised_learning.ipynb) or the HTML file 
i've provided in the repository. 

# How to run in a virtual environment   
Make the virtual environment. When working in your own system   

```bash
python3 -m venv phylo-env   
```     
Activate the virtual environment.   

```bash
source phylo-env/bin/activate   
```   

Install packages.  

```bash
make install run_script
```   

# How to run the docker   

Build DockerFile  
```bash
sudo docker build -t phylo-exp .
```   

Run the Docker image   
```bash
sudo docker run -it -p 8888:8888 --rm phylo-exp:latest
```  

# Add [data version control](https://dvc.org/doc/install)

Initialize dvc to the folder. To allow us to use dvc functionality to be used in the repo. NB. files will be created in the directory.    

```bash
dvc init
```  

Track changes to the different data files. The reason why we are doing this is because these files will change during the experiment especially if the investigator wants to try other experiments with more data.   

```bash
dvc add updated_data/data.md
dvc add updated_data/sequences.fasta
dvc add updated_data/sequence_metrics.csv
```  

Commit changes to save the changes that have occured in the repository.     

```bash
dvc commit
```

# FAQ
Q: What is a phylogenetic tree (aka phylogeny)? 
A: According to [Baum in Nature](https://www.nature.com/scitable/topicpage/reading-a-phylogenetic-tree-the-meaning-of-41956/#), this is a diagram that shows lines of evolutionary descent of different species, organisms or genes from a common ancestor. It is useful for organizing knowledge of biological diversity, for structuring classifications, or understanding evolutionary events.   

Q: What is antimicrobial resistance?
A: This is when a microorganism is no longer sensitive to action of antibiotics. Antibiotics are compounds produced by the natural metabollic processes of microbes for instance fungi that kill or inhibit other microbes. Key people who discovered the first named them are Alexander Fleming & Selman Waksman. Now, these microbes are not affected by these antibiotics and they can lose/gain function via misspelling their genetic material to also bypass the mechanism of action of antibiotics.Some microbes can also share genes and may get this resistance via horizontal gene transfer.   

Q: What are nitrogenous bases?
A: Adenine, Guanine, Cytosine and Thymine (Uracil, if in messanger ribonucleic acid form). These are the building blocks of DNA (Deoxyribonucleic acid) or the variant (Ribonucleic acid). Notice the sugar is how distinguish them deoxyribose and ribose. 

Q: What is 16S rRNA (ribosomal ribonucleic acid)?
A: It is component of the 30S subunit of the prokaryotic ribosome. It consists about 1500 nucleotides. It possesses a conserved and variable region which is important in Phylogenetic studies due to [slow evolution of the area.](https://en.wikipedia.org/wiki/16S_ribosomal_RNA) Similar regions of interest in other organisms are Internal transcribed spacer for fungi, 18S for microbial eukaryotes and for you + me 28S and 18S fragments.   

Q: What are GC and GA content?
A: These are common bioinformatic metrics that can enable us to find out the origin of a DNA we don't know the origin of. Higher GC content could mean contamination of a sample with a microbe. For most organisms its about 50 percent with some exceptions in some areas. Read [more.](https://rosalind.info/glossary/gc-content/)  

Q: What are Kmers?
A: A nucleotide sequence of a certain length. Read more [here](https://www.biostars.org/p/286438/) & [here.](https://rosalind.info/glossary/k-mer-composition/)
