#Intalar RTools

#Es importante instalar lo siguiente con permisos de administrado
#install.packages("tidyverse")
#install.packages("BiocManager")
#BiocManager::install("dada2")
#BiocManager::install("Biostrings")
#BiocManager::install("ShortRead")


library(BiocManager)
library(dada2)
library(ShortRead)
library(Biostrings)
library(here)

packageVersion("BiocManager")
packageVersion("dada2")
packageVersion("ShortRead")
packageVersion("Biostrings")

path <-  file.path("C:/Users/sgonz/Desktop/Escuela","Genomica", "codigo") 
path.cut <- file.path(path, "Zhang2022-seqs-ITS")
list.files(path)

cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

# Extraemos los nombres los cortes
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#plotQualityProfile(cutFs[1:2])
#plotQualityProfile(cutRs[1:2])


filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#Ser pacientes con esta linea
# out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
#                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, 
#                      multithread = TRUE)  
# save(out, file = file.path(path.cut, "dada", "out.Rdata"))
load(file = file.path(path.cut, "dada", "out.Rdata"))
head(out)

#Estimamos loa errores y los graficamos
# errF <- learnErrors(filtFs, multithread = TRUE)
# save(errF, file = file.path(path.cut, "dada", "errF.Rdata"))
load(file = file.path(path.cut, "error", "errF.Rdata"))
# errR <- learnErrors(filtRs, multithread = TRUE)
# save(errR, file = file.path(path.cut, "dada", "errR.Rdata"))
load(file = file.path(path.cut, "error", "errR.Rdata"))
#plotErrors(errR, nominalQ = TRUE)


#Dereplicar lecturas idénticas
# derepFs <- derepFastq(filtFs, verbose = TRUE)

# save(derepFs, file = file.path(path.cut, "dada", "derepFs.Rdata"))
load(file = file.path(path.cut, "dada", "derepFs.Rdata"))
# derepRs <- derepFastq(filtRs, verbose = TRUE)
# save(derepRs, file = file.path(path.cut, "dada", "derepRs.Rdata"))
load(file = file.path(path.cut, "dada", "derepRs.Rdata"))
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Inferencia de muestra
# dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
# save(dadaFs, file = file.path(path.cut, "dada", "dadaFs.Rdata"))
load(file = file.path(path.cut, "dada", "dadaFs.Rdata"))
# dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
# save(dadaRs, file = file.path(path.cut, "dada", "dadaRs.Rdata"))
load(file = file.path(path.cut, "dada", "dadaRs.Rdata"))

#Combinamos las lecturas emparejadas
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# save(mergers, file = file.path(path.cut, "dada", "mergers.Rdata"))
load(file = ile.path(path.cut, "dada", "mergers.Rdata"))


#Construir tabla de secuencia
# seqtab <- makeSequenceTable(mergers)
# save(seqtab, file.path(path.cut, "dada", "seqtab.Rdata"))
load(file = file.path(path.cut, "dada", "seqtab.Rdata"))
dim(seqtab)

#Quitamos las quimeras
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#save(seqtab.nochim, file = file.path(path.cut, "dada", "seqtab.nochim.Rdata"))
load(file = file.path(path.cut,  "dada", "seqtab.nochim.Rdata"))



#Track reads through the pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#Asignar taxonomia
load(file = file.path(path.cut, "dada", "taxa.Rdata"))
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
