# Description

This GitHub repository contains the scripts and pipelines used to generate bioinformatic data related to the analysis of *Drosophila* homologs of 16p12.1 genes.

There are two directories in this repository: 

1. Pipelines for identifying differentially-expressed fly genes from RNA-Sequencing data of knockdown of homologs of 16p12.1 genes. 
  * This directory contains a batch script for aligning and quantifying read counts using [TopHat2 v.2.1.1](https://ccb.jhu.edu/software/tophat/index.shtml) and [HTSeq v.0.6.1](https://htseq.readthedocs.io/en/release_0.11.1/), and an R pipeline for identifying differentially-expressed genes using [edgeR v.3.20.1](https://bioconductor.org/packages/release/bioc/html/edgeR.html).  
  * The raw RNA-Seq reads and quantified read counts for biological replicates are available at [NCBI GEO accession number GSE151330](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE151330).

2. Network analysis of CNV genes and simulated random gene sets.
 * This directory contains a Python script for analyzing the connectivity of genes within a human brain-specific interaction network. 
 * The script (`nearest_neighbor_weighted_allgenes.py`) takes an individual gene as standard input (gene name and Entrez ID), and outputs the shortest distance to every other gene in the network, as well as the genes located in the shortest path between the two target genes. For example, to find the connectivity of POLR3E (Entrez ID 55718) to all other genes in the genome, one would run:
 `python nearest_neighbor_weighted_allgenes.py POLR3E_55718`
 * The network file used in this analysis (`brain.degnorm-ge2.prob-gept02.dat`) is described in [Greene *et al, Nat. Genet.* 2015](https://www.ncbi.nlm.nih.gov/pubmed/25915600) and [Krishnan *et al, Nat. Neurosci.* 2016](https://www.ncbi.nlm.nih.gov/pubmed/27479844); we generated a sub-network that only contained edges with weights >2.0 (the top 0.5% of interactions in the network).

Bash pipeline scripts can be run in any Unix environment. R scripts can be run using any R version (scripts were generated using R v.3.4.2.). Python scripts for network analysis can be run in Python2 (scripts were generated using Python v.2.7.16) and require the [NetworkX package v.2.4](https://networkx.github.io/). 

# Copyright/License
The code in this repository is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the code.  If not, see <https://www.gnu.org/licenses/>.

# Contact
For questions or comments, please contact Matthew Jensen (mpj5142@psu.edu) or Santhosh Girirajan (sxg47@psu.edu).