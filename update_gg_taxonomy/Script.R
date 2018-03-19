source("http://bioconductor.org/biocLite.R")

biocLite("Biostrings")

library("Biostrings")

install.packages("devtools")
library(devtools)

install_github("mhahsler/rBLAST")

library("rBLAST")

Sys.which("blastn")




untar("C:/Users/marce/Desktop/16SMicrobial.tar.gz", exdir = "16SMicrobialDB")

microbial_database <- blast(db = "./16SMicrobialDB/16SMicrobial")

microbial_database


data_fasta <- readDNAStringSet("C:/Users/marce/Desktop/dna-sequences.fasta")

data_fasta

length(data_fasta)

data_predicted <- predict(microbial_database, data_fasta[1,])

data_predicted["QueryID"]

head(data_df)

data_predicted[1:5,]





Sys.setenv(ENTREZ_KEY = "ed4870836e8f61529227d9176a7c4a994c07")

getkey(service = "entrez")

genbank2uid(id = 'NR_041263.1')



ncbi_get_taxon_summary("304207")

clas <- classification("304207", db = "ncbi")

clas[[1]][2]

