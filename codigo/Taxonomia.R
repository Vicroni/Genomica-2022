library(BiocManager)
library(dada2)
library(ShortRead)
library(Biostrings)
library(here)

packageVersion("BiocManager")
packageVersion("dada2")
packageVersion("ShortRead")
packageVersion("Biostrings")

path <-  file.path("~", "Desktop", "Genomica") 
path.cut <- file.path(path, "Zhang2022-seqs-ITS")

load(file = file.path(path.cut, "dada", "seqtab.nochim.Rdata"))
path.cut.ref <- file.path(path.cut, "ref", "sh_general_release_dynamic_s_01.12.2017.fasta")
taxa <- assignTaxonomy(seqtab.nochim, path.cut.ref, multithread = TRUE)
save(taxa, file =  file.path(path.cut, "dada", "taxa.Rdata"))

