[![Phylogenetic tree study Continuous Integration/Continuos Delivery](https://github.com/Shuyib/Phylogenetic-tree-study/actions/workflows/devops.yml/badge.svg)](https://github.com/Shuyib/Phylogenetic-tree-study/actions/workflows/devops.yml)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Shuyib/Phylogenetic-tree-study/HEAD)

# Phylogenetic-tree-study

Estimating Phylogenetic trees using 30 microorganisms (previously 6 organisms: review [data folder](https://github.com/Shuyib/Phylogenetic-tree-study/tree/master/data) and [notebook](https://github.com/Shuyib/Phylogenetic-tree-study/blob/master/Phylogenetic_trees_unsupervised_learning.ipynb). Looking at the 16S rRNA gene with Unsupervised Learning, web based tools and Molecular Evolutionary Genetics Analysis MEGA7. Further we are looking at motifs and finding out what they do.

It is important to know these regions since they can potentially give use clues about the regions we can target for targeted DNA therapies.

# How to run in a virtual environment   
Make the virtual environment. When working in your own system   

```bash
python3 -m venv phylo-env   
```     
Activate the virtual environment.   

```bash
source phylo-env/bin/activate   
```   

Install packages. You need an email to be in your .bashrc file to run biopython.   

```bash
make install run_script
```   

In your terminal, in the directory where you cloned this repository. Run this command to run notebooks.   

```bash
jupyter notebook Phylogenetic_trees_unsupervised_learning.ipynb
```

Previously, we've not provided a codebook/data description file since one of the headings cover that in the notebook.
Otherwise, you can check out the [notebook](https://nbviewer.jupyter.org/github/Shuyib/Phylogenetic-tree-study/blob/master/Phylogenetic_trees_unsupervised_learning.ipynb) or the HTML file 
i've provided in the repository.    

# How to run project in Docker

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
