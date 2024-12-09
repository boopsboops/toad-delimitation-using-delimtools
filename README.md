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
scripts/download-sequences.R -c Amazophrynella -n 400 -x 1000 -b 200 -a false -d false
scripts/download-sequences.R -c Dendrophryniscus_proboscideus -n 400 -x 1000 -b 5 -a true -d false
scripts/clean-and-cluster.R -n 10 -c 0.6 -m 2 -d false
scripts/pick-clusters.R -c 12 -g 16S
scripts/annotate-ncbi.R -t 1 -c ncbi
scripts/filter-species.R -n 2 -i true
scripts/align-trim-concatenate.R -p 0.2 -t 2 -i true
scripts/tree-search.R -m TN93+G -v false -e 0.1 -t 4
scripts/tree-plot.R -w 1 -h 4 -s 1.5
scripts/tidy-results-directory.R
```


```r
library("here")
library("glue")
library("ape")
library("tidyverse")
library("delimtools")
library("haplotypes")
library("spider")

# get latest dir
latest.dir <- sort(list.dirs(here("ncbi-supermatrix/temp"),full.names=FALSE,recursive=FALSE),decreasing=TRUE)[1]
# path
tr.path <- here::here("ncbi-supermatrix/temp",latest.dir,"trees/16S.aligned.trimmed.fasta.raxml.bestTree")
fas.path <- here::here("ncbi-supermatrix/temp",latest.dir,"alignments/16S.aligned.trimmed.fasta")
tab.path <- here::here("ncbi-supermatrix/temp",latest.dir,"metadata/ncbi-clean.csv")



# collapse haps


tr <- ape::read.tree(tr.path)
fas <- ape::read.FASTA(fas.path)

outgroup.node <- ape::getMRCA(tr,c("JN867566.1","KU495200.1"))
outgroup.tips <- extract.clade(tr,node=outgroup.node)$tip.label


fas.haps <- delimtools::hap_collapse(fas,collapseSubstrings=TRUE,clean=TRUE)
fas.haps <- fas.haps[!names(fas.haps) %in% outgroup.tips]
fas.path.haps <- glue::glue("{fas.path}.haps.fasta")
fas.haps |> ape::write.FASTA(fas.path.haps)


tr.rooted.haps <- tr |> 
    ape::root(node=outgroup.node,resolve.root=TRUE) |> 
    ape::keep.tip(tip=names(fas.haps))

#tr.rooted |> plot()

tr.path.root <- glue("{tr.path}.rooted.haps.nwk")
tr.rooted.haps |> write.tree(tr.path.root)


# mptp
delimtools::min_brlen(tree=tr.path.root,n=10)
mptp.df <- delimtools::mptp_tbl(infile=tr.path.root,exe=here::here("software/mptp/bin/mptp"),method="single",minbrlen=0.002)#minbrlen=0.01
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


# parsimnet
pnet <- haplotypes::parsimnet(haplotypes::as.dna(as.matrix(fas.haps)),indels="5th",prob=.95)
pnet.df <- delimtools:::parsimnet_tbl(dna=fas.haps, parsimnet=pnet)
pnet.df |> delimtools::report_delim()

```