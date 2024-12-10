# Little Toads, Big Diversity

Frog Species Delimitation using _delimtools_
Delimitation of frogs for Herp lab meeting


### Install software

```bash
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
# run ncbi-supermatrix (cd ncbi-supermatrix first)
scripts/download-sequences.R -c Amazophrynella -n 500 -x 600 -b 200 -a false -d false
scripts/download-sequences.R -c Melanophryniscus_stelzneri -n 500 -x 600 -b 10 -a true -d false
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2 -d false
# choose cluster
scripts/pick-clusters.R -c 5 -g 16s
scripts/annotate-ncbi.R -t 1 -c ncbi
scripts/filter-species.R -n 2 -i true
scripts/align-trim-concatenate.R -p 0.05 -t 2 -i true
scripts/tree-search.R -m GTR+G -v false -e 0.1 -t 4
scripts/tree-plot.R -w 1 -h 4 -s 1.5
scripts/tidy-results-directory.R
```


```r
# load libs
source(here::here("scripts/load-libs.R"))

# get latest dir and load paths to files
latest.dir <- sort(list.dirs(here::here("ncbi-supermatrix/temp"),full.names=FALSE,recursive=FALSE),decreasing=TRUE)[1]
tr.path <- here::here("ncbi-supermatrix/temp",latest.dir,"trees/16s.aligned.trimmed.fasta.raxml.bestTree")
beast.path <- here::here("assets/16s.aligned.trimmed.haps.tre")
fas.path <- here::here("ncbi-supermatrix/temp",latest.dir,"alignments/16s.aligned.trimmed.fasta")
tab.path <- here::here("ncbi-supermatrix/temp",latest.dir,"metadata/ncbi-clean.csv")
fas.path.haps <- glue::glue("{fas.path}.haps.fasta")
tr.path.root <- glue::glue("{tr.path}.rooted.haps.nwk")


# load trees and fasta, data
fas <- ape::read.FASTA(fas.path)
tr <- ape::read.tree(tr.path)
beast.tr <- treeio::read.beast(beast.path)
ncbi.tab <- readr::read_csv(tab.path,show_col_types=FALSE)


# collapse haplotypes and write out fasta
fas.haps <- delimtools::hap_collapse(fas,collapseSubstrings=TRUE,clean=TRUE)
fas.haps |> ape::write.FASTA(fas.path.haps)

# root and drop haplotypes
tr.rooted.haps <- tr |> 
    castor::root_in_edge(root_edge=which.max(tr$edge.length)) |> 
    ape::keep.tip(tip=names(fas.haps))
# tr.rooted.haps |> ape::ladderize() |> plot()
tr.rooted.haps |> ape::write.tree(tr.path.root)


##### SPECIES DELIMITATION #####

# gmyc
ape::is.binary(beast.tr@phylo)
set.seed(42)
gmyc.res <- splits::gmyc(beast.tr@phylo,method="single",interval=c(0,5),quiet=FALSE)
summary(gmyc.res)
gmyc.df <- delimtools::gmyc_tbl(gmyc.res)
gmyc.df |> delimtools::report_delim()


# bgmyc
set.seed(42)
bgmyc.res.single <- bGMYC::bgmyc.singlephy(beast.tr@phylo,mcmc=11000,burnin=1000,thinning=100,t1=2,t2=length(beast.tr@phylo$tip.label),start=c(1,0.5,50))
bgmyc.df <- delimtools::bgmyc_tbl(bgmyc.res.single,ppcutoff=0.05)
bgmyc.df |> delimtools::report_delim()


# mptp
delimtools::min_brlen(tree=tr.path.root,n=10)
mptp.df <- delimtools::mptp_tbl(infile=tr.path.root,exe=here::here("software/mptp/bin/mptp"),method="single",minbrlen=0.0019)#minbrlen=0.01
mptp.df <- mptp.df |> rename(mptp=ptp)
mptp.df |> delimtools::report_delim()


# asap
asap.df <- delimtools::asap_tbl(infile=fas.path.haps,exe=here::here("software/ASAP/bin/asap"),model=3)
asap.df |> delimtools::report_delim()


# abdg
abgd.df <- delimtools::abgd_tbl(infile=fas.path.haps,slope=0.5,exe=here::here("software/Abgd/bin/abgd"),model=3)
abgd.df |> delimtools::report_delim()


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


# parsimnet
pnet <- haplotypes::parsimnet(haplotypes::as.dna(as.matrix(fas.haps)),indels="5th",prob=.95)
pnet.df <- delimtools:::parsimnet_tbl(dna=fas.haps, parsimnet=pnet)
pnet.df |> delimtools::report_delim()
nrow(pnet.df) == nrow(as.matrix(fas.haps))
names(fas.haps)[-which(names(fas.haps) %in% pull(pnet.df,labels))]


# morph
ncbi.tab.sub <- ncbi.tab |> filter(gbAccession %in% names(fas.haps))
morph.df <- delimtools::morph_tbl(labels=dplyr::pull(ncbi.tab.sub,gbAccession),sppVector=dplyr::pull(ncbi.tab.sub,scientificName))
morph.df
morph.df |> delimtools::report_delim()


# join delims
all.delims.df <- delimtools::delim_join(list(mptp.df,gmyc.df,bgmyc.df,locmin.df,fixed.df,asap.df,abgd.df,morph.df))#pnet.df,
all.delims.df |> delimtools::report_delim()
all.delims.df |> delim_consensus(n_match=3)


tip.tab <- ncbi.tab.sub |> 
    dplyr::mutate(labs=glue::glue("{gbAccession} | {scientificName}")) |> 
    dplyr::select(gbAccession,labs,scientificName)

cols <- delimtools::delim_brewer(delim=all.delims.df,package="viridisLite",palette="viridis",seed=42)

p <- delimtools::delim_autoplot(delim=all.delims.df,tr=beast.tr,tbl_labs=tip.tab,col_vec=cols,hexpand=0.3,widths=c(0.4,0.1),n_match=3,delim_order=c("asap","abgd","locmin","fixed","gmyc","bgmyc","mptp","morph"),consensus=TRUE)
ggplot2::ggsave(here::here("temp/amazophrynella-delimitation.pdf"),plot=p,height=500,width=400,units="mm")


# autoplot2
p <- delimtools::delim_autoplot2(delim=all.delims.df,tr=beast.tr, consensus=TRUE, n_match=3, tbl_labs=tip.tab, species="scientificName",hexpand= 0.1, widths= c(0.5, 0.2), delim_order=c("asap","abgd","locmin","fixed","gmyc","bgmyc","mptp","morph"))
ggplot2::ggsave(here::here("temp/amazophrynella-delimitation-bw.pdf"),plot=p,height=500,width=400,units="mm")



```