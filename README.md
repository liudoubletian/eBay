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
*`otu.data` : otu table which was a n*m matrix including n samples and m taxa
* `group` : the group for each sample  
* `test.method` : a t-test or wilcoxon ramk sum test  
* `cutf` : a significance level  
it returns a list of results:
* `final.p` : the adjusted p values 
* `dif.otus` : the detected differential abundance otus  

# Example

The following function shows how to simulate data from a dirichlet multinomial distribution.  
set.seed(1)  
rand_pi <- runif(20)   
control_pi = case_pi = rand_pi/sum(rand_pi)   
control_theta = case_theta = 0.1  
group <- rep(c(0,1),each =20)  
ntree_table <- simulation_dm(p=20,seed=1, N=20,control_pi, case_pi,control_theta,case_theta)  

We can run the eBay function to normalize the simulated data and return the detected differential abundance taxa.  

ebay.res <- eBay(otu.data=ntree_table, group=group, test.method="t", cutf=0.05)  
ebay.res  








