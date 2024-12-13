# Little Toads, Big Diversity

Frog Species Delimitation using _delimtools_
Delimitation of frogs for Herp lab meeting


### Install software

```bash
########################
### Install software ###
########################

# install code
git clone https://github.com/boopsboops/frog-delimitation-using-delimtools.git
cd frog-delimitation-using-delimtools
Rscript -e "renv::restore()"
mkdir temp temp-local
# SYMLINK a software dir
ln -s ~/Software/RSoftware software
```

```bash
# install ncbi-supermatrix 
git clone https://github.com/boopsboops/ncbi-supermatrix.git
cd ncbi-supermatrix
Rscript -e "renv::restore()"
```


### Get GenBank data

```bash
########################
### NCBI SUPERMATRIX ###
########################

# run ncbi-supermatrix (cd ncbi-supermatrix first)
scripts/download-sequences.R -c Amazophrynella -n 500 -x 600 -b 200 -a false -d false
scripts/download-sequences.R -c Melanophryniscus_stelzneri -n 500 -x 600 -b 10 -a true -d false
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2 -d false
# choose cluster
scripts/pick-clusters.R -c 5 -g 16s
scripts/annotate-ncbi.R -t 1 -c ncbi
scripts/filter-species.R -n 2 -i true
scripts/align-trim-concatenate.R -p 0.05 -t 2 -i true
scripts/tree-search.R -m TN93+G -v false -e 0.1 -t 4
scripts/tree-plot.R -w 0.5 -h 3 -s 1.5
scripts/tidy-results-directory.R
cd ..
```


```r
########################
### SETUP DIRS/FILES ###
########################

# load libs
source(here::here("scripts/load-libs.R"))

# get latest dir and load paths to files
latest.dir <- sort(list.dirs(here::here("ncbi-supermatrix/temp"),full.names=FALSE,recursive=FALSE),decreasing=TRUE)[1]
tr.path <- here::here("ncbi-supermatrix/temp",latest.dir,"trees/16s.aligned.trimmed.fasta.raxml.bestTree")
beast.path <- here::here("assets/16s.aligned.trimmed.haps.tre")
#beast.trees.path <- here::here("assets/16s.aligned.trimmed.haps.run1.run2.trees")
fas.path <- here::here("ncbi-supermatrix/temp",latest.dir,"alignments/16s.aligned.trimmed.fasta")
tab.path <- here::here("ncbi-supermatrix/temp",latest.dir,"metadata/ncbi-clean.csv")
fas.path.haps <- glue::glue("{fas.path}.haps.fasta")
tr.path.root <- glue::glue("{tr.path}.rooted.haps.nwk")

# load trees and fasta, data
fas <- ape::read.FASTA(fas.path)
tr <- ape::read.tree(tr.path)
beast.tr <- treeio::read.beast(beast.path)
ncbi.tab <- readr::read_csv(tab.path,show_col_types=FALSE)


########################
#### COLLAPSE HAPS  ####
########################

# collapse haplotypes and write out fasta
fas.haps <- delimtools::hap_collapse(fas,collapseSubstrings=TRUE,clean=TRUE)
fas.haps |> ape::write.FASTA(fas.path.haps)

# root and drop haplotypes
tr.rooted.haps <- tr |> 
    castor::root_in_edge(root_edge=which.max(tr$edge.length)) |> 
    ape::ladderize() |>
    ape::keep.tip(tip=names(fas.haps))
tr.rooted.haps |> ape::ladderize() |> plot()
tr.rooted.haps |> ape::write.tree(tr.path.root)


#####################################
#### SPECIES DELIMITATION - GMYC ####
#####################################

# gmyc
ape::is.binary(beast.tr@phylo)
set.seed(42)
gmyc.res <- splits::gmyc(beast.tr@phylo,method="single",interval=c(0,5),quiet=FALSE)
summary(gmyc.res)
gmyc.df <- delimtools::gmyc_tbl(gmyc.res)
gmyc.df |> delimtools::report_delim()


######################################
#### SPECIES DELIMITATION - bGMYC ####
######################################

# bgmyc
set.seed(42)
bgmyc.res.single <- bGMYC::bgmyc.singlephy(beast.tr@phylo,mcmc=11000,burnin=1000,thinning=100,t1=2,t2=length(beast.tr@phylo$tip.label),start=c(1,0.5,50))
bgmyc.df <- delimtools::bgmyc_tbl(bgmyc.res.single,ppcutoff=0.05)
bgmyc.df |> delimtools::report_delim()


#####################################
#### SPECIES DELIMITATION - mPTP ####
#####################################

# mptp
delimtools::min_brlen(tree=tr.path.root,n=10)
mptp.df <- delimtools::mptp_tbl(infile=tr.path.root,exe=here::here("software/mptp/bin/mptp"),method="single",minbrlen=0.0019)#minbrlen=0.01
mptp.df <- mptp.df |> rename(mptp=ptp)
mptp.df |> delimtools::report_delim()


#####################################
#### SPECIES DELIMITATION - ASAP ####
#####################################

# asap
asap.df <- delimtools::asap_tbl(infile=fas.path.haps,exe=here::here("software/ASAP/bin/asap"),model=3)
asap.df |> delimtools::report_delim()


#####################################
#### SPECIES DELIMITATION - ABGD ####
#####################################

# abdg
abgd.df <- delimtools::abgd_tbl(infile=fas.path.haps,slope=0.5,exe=here::here("software/Abgd/bin/abgd"),model=3)
abgd.df |> delimtools::report_delim()


#######################################
#### SPECIES DELIMITATION - LOCMIN ####
#######################################

# locmin
mat <- ape::dist.dna(fas.haps,model="raw",pairwise.deletion=TRUE)
lmin <- spider::localMinima(as.matrix(mat))
plot(lmin); abline(v=lmin$localMinima[1],col="red")
locmin.df <- delimtools::locmin_tbl(mat,threshold=lmin$localMinima[1])
locmin.df |> delimtools::report_delim()

# 2%
fixed.df <- delimtools::locmin_tbl(mat,threshold=0.02)
fixed.df <- fixed.df |> rename(fixed=locmin)
fixed.df |> delimtools::report_delim()


#######################################
#### SPECIES DELIMITATION - MORPHO ####
#######################################

# morph
ncbi.tab.sub <- ncbi.tab |> filter(gbAccession %in% names(fas.haps))
morph.df <- delimtools::morph_tbl(labels=dplyr::pull(ncbi.tab.sub,gbAccession),sppVector=dplyr::pull(ncbi.tab.sub,scientificName))
morph.df
morph.df |> delimtools::report_delim()


##########################################
#### SPECIES DELIMITATION - PARSIMNET ####
##########################################

# collapse haps with parsimnet
pnet <- haplotypes::parsimnet(haplotypes::as.dna(as.matrix(fas)),indels="missing",prob=0.95)
pnet.df <- delimtools:::parsimnet_tbl(dna=fas, parsimnet=pnet)
pnet.df |> delimtools::report_delim()
#nrow(pnet.df) == nrow(as.matrix(fas))
#names(fas)[-which(names(fas) %in% pull(pnet.df,labels))]
#haps <- pull(pnet.df,labels)
#fas.haps <- fas[haps]


############################
#### JOIN AND SUMMARISE ####
############################

# join delims
delimtools::delim_join(list(mptp.df,gmyc.df,bgmyc.df,locmin.df,fixed.df,asap.df,abgd.df,morph.df,pnet.df))
all.delims.df <- delimtools::delim_join(list(mptp.df,gmyc.df,bgmyc.df,locmin.df,fixed.df,asap.df,abgd.df,morph.df))
all.delims.df |> delimtools::report_delim()
all.delims.df |> delim_consensus(n_match=3)

# match ratio congruence
all.delims.df |> delimtools::match_ratio() |> dplyr::arrange(desc(match_ratio)) |> print(n=Inf)


#######################
#### PLOT THE TREE ####
#######################

# make tip labels
tip.tab <- ncbi.tab.sub |> 
    dplyr::mutate(labs=glue::glue("{gbAccession} | {scientificName}")) |> 
    dplyr::select(gbAccession,labs,scientificName)


# get colour palettes
cols <- delimtools::delim_brewer(delim=all.delims.df,package="viridisLite",palette="viridis",seed=42)
cols <- delimtools::delim_brewer(delim=all.delims.df,package="RColorBrewer",palette="Set2",seed=42)
cols <- delimtools::delim_brewer(delim=all.delims.df,package="randomcoloR",seed=42)

# plot
p <- delimtools::delim_autoplot(delim=all.delims.df,tr=beast.tr,tbl_labs=tip.tab,col_vec=cols,hexpand=0.3,widths=c(0.4,0.1),n_match=3,delim_order=c("asap","abgd","locmin","fixed","gmyc","bgmyc","mptp","morph"),consensus=TRUE)
ggplot2::ggsave(here::here("temp/amazophrynella-delimitation.pdf"),plot=p,height=500,width=400,units="mm")


# autoplot2
p <- delimtools::delim_autoplot2(delim=all.delims.df,tr=beast.tr, consensus=TRUE, n_match=3, tbl_labs=tip.tab, species="scientificName",hexpand= 0.1, widths= c(0.5, 0.2), delim_order=c("asap","abgd","locmin","fixed","gmyc","bgmyc","mptp","morph"))
ggplot2::ggsave(here::here("temp/amazophrynella-delimitation-bw.pdf"),plot=p,height=500,width=400,units="mm")


# confidence intervals
#source("https://raw.githubusercontent.com/legalLab/delimtools/refs/heads/main/R/confidence_intervals.R")
#beast.trees <- ape::read.nexus(beast.trees.path)
#confidence_intervals(method=, tr=beast.tr@phylo, posterior=beast.trees[1:10], method="single", interval=c(0,5))

```