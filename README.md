 eBay (An empirical Bayes normalization method for microbiome data)

# Introduction
The characteristic of microbiome data were complicated: varied sequencing depth, lots of zeros, over-dispersion,
a phylogenetic tree representing the relationship among taxa and compositionality.

The characteristics of microbiome data are complicated: varied sequencing depth, lots of zeros, over-dispersion, a phylogenetic tree representing the relationship among taxa, and compositionality. In order to deal with these problems, we proposed an empirical Bayes normalization method called eBay. [1] shows the statistical model in detail.

# Installation
```r
install.packages("devtools")  
devtools::install_github("liudoubletian/eBay")  
library(eBay)  
```
# Basic Usage
## normalization in the absence of a phylogenetic tree
```r
eBay(otu.data=ntree_table, group=group, test.method="t", cutf=0.05)
```
* `otu.data` : an OTU table with n rows (samples) and m columns (taxa)
* `group` : a n-vector of group indicators
* `test.method` : t-test or Wilcoxon rank sum test
* `cutf` : level of significance

it returns a list of results:  
* `final.p` : the adjusted p values 
* `dif.otus` : the set of differentially abundant OTUs  

## phylogeny-aware empirical Bayes normalization
```r
eBay_tree(otu.data=tree_table,tree=tree,group=group,test.method="t",cutf=0.05)
```
* `otu.data` : an OTU table with n rows (samples) and m columns (taxa)
* `tree` : a phylogenetic tree among m taxa
* `group` : a n-vector of group indicators
* `test.method` : t-test or Wilcoxon rank sum test
* `cutf` : level of significance

it returns a list of results:  
* `final.p` : the adjusted p values 
* `dif.otus` : the set of differentially abundant OTUs  
# Example
## simulation from DM
The following function shows how to simulate data from a Dirichlet-multinomial distribution.  
```r
set.seed(1)  
rand_pi <- runif(20)   
control_pi = case_pi = rand_pi/sum(rand_pi)   
control_theta = case_theta = 0.1  
group <- rep(c(0,1),each =20)  
ntree_table <- simulation_dm(p=20,seed=1, N=20,control_pi, case_pi,control_theta,case_theta)  
```

Run the eBay function to normalize the data and return a set of differentially abundant taxa.

```r
ebay.res <- eBay(otu.data=ntree_table, group=group, test.method="t", cutf=0.05)  
ebay.res  
```
## simulation from DTM
First, generate a tree randomly.
```r
p <- 40
set.seed(1)
tree <- simulate_tree(p)
```

The following function shows how to simulate data from the Dirichlet-tree multinomial model.
```r
set.seed(1)    
control_pi = case_pi = c()
for(j in (p+1):(p+tree$Nnode)){
   set.seed(j)
   random_pi <- runif(1,0.2,0.4)
   control_pi[which(tree$edge[,1]==j)] <- c(random_pi, 1-random_pi)
   case_pi[which(tree$edge[,1]==j)] <- c(random_pi, 1-random_pi)
}
control_theta = case_theta = rep(0.1, tree$Nnode) 
group <- rep(c(0,1),each =20)  
tree_table <- simulation_dtm(p=40,tree, seed=1, N=20,control_pi, case_pi,control_theta,case_theta)  
```
Run the eBay_tree function to normalize the data and return a set of differentially abundant taxa.
```r
ebay_tree.res <- eBay_tree(otu.data=tree_table, tree=tree, group=group, test.method="t", cutf=0.05)  
ebay_tree.res  
```




