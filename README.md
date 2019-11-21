eBay (empirical bayes normalization method for microbiome data)
# Introduction
The characteristic of microbiome data were complicated: varied sequencing depth, lots of zeros, over-dispersion,
a phylogenetic tree representing the relationship among taxa and compositionality.
In order to deal with several problems, we proposed a normalization method based on a bayes method called eBay.
[1] shows the statistical model in detail.

# Installation

install.packages("devtools")  
devtools::install_github("liudoubletian/eBay")  
library(eBay)  

# Basic Usage
eBay(otu.data=ntree_table, group=group, test.method="t", cutf=0.05)  
* otu.data : otu table which was a n*m matrix including n samples and m taxa
* group : the group for each sample  
* test.method : a t-test or wilcoxon ramk sum test  
* cutf : a significance level  
it returns a list of results:
* final.p : the adjusted p values 
* dif.otus : the detected differential abundance otus





