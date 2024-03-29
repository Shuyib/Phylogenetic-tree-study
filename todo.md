- [ ] Update project with more microbes to improve the grouping
- [ ] Update table created with the old microbes and descriptions
- [ ] Add todo list with updated kanban
- [ ] Add project workflow (last steps)
- [x] Add Makefile
- [x] Update the packages to the only ones useful to project
- [ ] Update environment.yml to have all the packages based on requirements.txt (use miniconda)
- [x] Add Dockerfile
- [x] Add License
- [x] Add nine potential organisms to the 16S rRNA (*Neisseria gonorrhea*(https://www.ncbi.nlm.nih.gov/nuccore/X07714.1/ & paper normally streptomycin resistant https://www.ncbi.nlm.nih.gov/pmc/articles/PMC336862/?page=1)[x], *Bordtella species*(https://www.ncbi.nlm.nih.gov/nuccore/AF177667 [x] & https://www.ncbi.nlm.nih.gov/nuccore/AF177666 [x] paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC86205/), *Mycobacterium tuberculosis* (https://www.ncbi.nlm.nih.gov/nuccore/MN166772.1 & paper unpublished), *Salmonella Typhi*(https://www.ncbi.nlm.nih.gov/nuccore/Z47544 & paper DOI: 10.13057/biodiv/d120101)[], *Mycoplasma genitalium*(potential https://www.ncbi.nlm.nih.gov/nuccore/X77334 and paper https://pubmed.ncbi.nlm.nih.gov/12517858/)[x], *Clostridium difficile* (https://www.ncbi.nlm.nih.gov/nuccore/MW798268 & paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8708398/)(Non-human) replace with this https://www.ncbi.nlm.nih.gov/nuccore/NR_112172.1[x], *Pseudomonas aeruginosa* (https://www.ncbi.nlm.nih.gov/nuccore/LT599799.2 & paper unpublished)[x] and Carbapenem-resistant Enterobacterales Raoultella ornithinolytica/ Klebsiella ornithinolytica(https://www.ncbi.nlm.nih.gov/nuccore/KT378441.1 & unpublished)[x]. Easier to go through search engine first to narrow down better and the nucleotide database in ncbi is easier.    
- [ ] Add CI/CD step with github actions.  
- [x] Use biopython or other libraries to try automate steps of getting sequences.   
      - [x] fix the fastadataloader to allow getting data for all the accession numbers and append them to a file.   
      - [x] Write a function to automate the getting of sequences with the accession numbers.  
      - [x] Merge fasta files to multifasta.  
- [ ] Use [MEME](https://rosalind.info/glossary/meme/) to find motifs and improve clustering done before. NB: Might want to check the [anatomy of a gene](https://learn.genetics.utah.edu/content/basics/geneanatomy/) and if you run into duplicates errors. Choose a few sequences at a time.  
- [ ] Add more different versions of python from 3.8, 3.9 and 3.10 in your tests. Just to circumvent the end of life python things to github actions.
- [x] Update Dockerfile to Python 3.8-slim-buster.
- [ ] Edit the requirements.txt on another branch for easy comparison.
- [ ] Consider using DictVectorizer and semi-supervised learning to see if any generalizations arise from using a neural network. Review contrastive loss and ideas [here](https://towardsdatascience.com/contrastive-loss-explaned-159f2d4a87ec)   
- [ ] Assign 2 microbes per person to speed up sequence discovery and replacement.  
- [x] Test the readFasta function from this [notebook](https://github.com/BioinfoNet/TeachingJupyterNotebooks).  
- [ ] Make a dashboard with [Pyscript](https://pyscript.net/) or [streamlit](streamlit.io) to explain the results better.  
- [ ] Use KEGG or something else to find out the function of the motifs found by MEME. Refer to biopython. 
