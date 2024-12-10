#!/usr/bin/env sh

# need to be in 'delimtools-testing'
# find most recent tempdir
tmpdir=$(ls -d ncbi-supermatrix/temp/*/ | sort -r | head -n 1)
echo $tmpdir
cd ${tmpdir}/alignments

# beauti
# beauti 16S.aligned.trimmed.fasta.haps.fasta
# set up beauti for two .xml files using strict clock, constant coalescent, TN93+G model, chain length 10000000 (10 million), log every 16000

# run beast
beast -beagle_auto -overwrite -seed 42 16s.aligned.trimmed.fasta.haps.fasta.run1.xml
beast -beagle_auto -overwrite -seed 24 16s.aligned.trimmed.fasta.haps.fasta.run2.xml
tracer 16s.aligned.trimmed.fasta.haps.fasta.run1.log 16s.aligned.trimmed.fasta.haps.fasta.run2.log 
logcombiner -trees -burnin 2016000 16s.aligned.trimmed.fasta.haps.fasta.run1.trees 16s.aligned.trimmed.fasta.haps.fasta.run2.trees 16s.aligned.trimmed.fasta.haps.fasta.run1.run2.trees
grep -c "tree STATE_" 16s.aligned.trimmed.fasta.haps.fasta.run1.run2.trees
treeannotator -burninTrees 0 -heights ca 16s.aligned.trimmed.fasta.haps.fasta.run1.run2.trees 16s.aligned.trimmed.fasta.haps.fasta.tre
figtree 16s.aligned.trimmed.fasta.haps.fasta.tre

# copy files
cp 16s.aligned.trimmed.fasta.haps.fasta.run1.run2.trees ../../../../assets/16s.aligned.trimmed.haps.run1.run2.trees
cp 16s.aligned.trimmed.fasta.haps.fasta.tre ../../../../assets/16s.aligned.trimmed.haps.tre

# run raxml on beast tree
# run in R first to convert tree to newick
#ape::read.nexus(here(today.path,"coi.geophagus.haps.beast.tre")) |> write.tree(here(today.path,"coi.geophagus.haps.beast.tre.nwk"))
#raxml-ng --evaluate --threads auto{} --tree coi.geophagus.haps.beast.tre.nwk --lh-epsilon 0.1 --redo --seed 42 --outgroup MH538063.1 --model TN93+G --msa coi.geophagus.haps.fasta