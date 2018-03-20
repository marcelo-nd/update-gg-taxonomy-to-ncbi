############ Install required packages ############

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("devtools", "dplyr", "taxize")
library("devtools")
install_github("mhahsler/rBLAST")




############ Load required packages ############

library("Biostrings")
library("rBLAST")
# Check if blast is in PATH.
Sys.which("blastn")
library("dplyr", lib.loc = "C:/Program Files/R/R-3.4.4/library")
library("taxize", lib.loc = "C:/Program Files/R/R-3.4.4/library")
# set NCBI's entrez api key
Sys.setenv(ENTREZ_KEY = "ed4870836e8f61529227d9176a7c4a994c07")
# Check that key variable is in path.
getkey(service = "entrez")





############ Load required data ############

# Extract reference database if not done already.
untar("C:/Users/marce/Desktop/16SMicrobial.tar.gz", exdir = "16SMicrobialDB")
# Load reference database.
microbial_database <- blast(db = "./16SMicrobialDB/16SMicrobial")
# Check database
microbial_database
# Load taxonomy table to be updated with NCBI's taxonomy.
tax_table <- read.table("C:/Users/marce/Desktop/taxonomy.tsv", sep = "\t", header = TRUE)
tax_table
# Load sequences to be analyzed. Ids corresponding to taxonomy OTUs.
data_fasta <- readDNAStringSet("C:/Users/marce/Desktop/dna-sequences.fasta")
# Check fasta sequences
data_fasta

data_small <- head(data_fasta, 3)
data_small

tax2 <- head(tax_table, 3)
tax2

new_tax2 <- tax2
new_tax2$Taxon <- as.character(new_tax2$Taxon)
new_tax2

############ Function for parsing the NCBI's taxonomy into greegngenes format ############

parse_ncbi_to_gg <- function(ncbi_tax) {
    return(sprintf("k__%s; p__%s; c__%s; o__%s; f__%s; g__%s; s__%s",
           # Kingdom
           filter(ncbi_tax[[1]], rank == "superkingdom")[1, 1],
           # Phylum
           filter(ncbi_tax[[1]], rank == "phylum")[1, 1],
           # Class
           filter(ncbi_tax[[1]], rank == "class")[1, 1],
           # Order
           filter(ncbi_tax[[1]], rank == "order")[1, 1],
           # Family
           filter(ncbi_tax[[1]], rank == "family")[1, 1],
           # Genus
           filter(ncbi_tax[[1]], rank == "genus")[1, 1],
           # Species
           strsplit(filter(ncbi_tax[[1]], rank == "species")[1, 1], " ")[[1]][2]))
}

############ Function for parsing the NCBI's taxonomy into greegngenes format ############

is_taxonomy_complete <- function(taxonomy) {
    print("true")
    return(TRUE)
}

############ Function for blasting sequences, retrieving taxonomy from NCBI taxonomy database and replacing an old taxonomy with the NCBI's taxonomy ############

blast_and_replace <- function(seq) {
    seq_id <- seq@ranges@NAMES[1]
    seq_taxonomy <- as.character((new_tax2 %>% filter(Feature.ID == seq_id) %>% select(Taxon))[1, 1])
    print(seq_id)
    print(seq_taxonomy)
    # If original taxonomy is incomplete.
    if (is_taxonomy_complete(seq_taxonomy)) {
        print("blasting")
        blast_result <- arrange(predict(microbial_database, seq), desc(Perc.Ident), desc(Bits), E)[1,]["SubjectID"]
        print(blast_result)
        ncbi_tax_classification <- classification(genbank2uid(id = blast_result[1, 1]), db = "ncbi")
        print(ncbi_tax_classification)
        print(parse_ncbi_to_gg(ncbi_tax_classification))
        new_tax2 <- new_tax2 %>% mutate(Taxon = replace(Taxon, which(Feature.ID == seq_id), parse_ncbi_to_gg(ncbi_tax_classification)))
    }
}

data_small[1,]

blast_and_replace(data_small[1,])


new_tax2

tax2$Taxon <- as.character(tax2$Taxon)
new_tax2 <- new_tax2 %>% mutate(Taxon = replace(Taxon, which(Feature.ID == "4fdb872697ff4712d1408c2a31c881ef"), "caca_caca"))

new_tax2

blast_result2 <- arrange(predict(microbial_database, data_small[1,]), desc(Perc.Ident), desc(Bits), E)[1,]["SubjectID"]

predict(microbial_database, data_small[1,])


s









for (i in 1:length(data_small)) {
    blast_and_replace(data_small[i,])
}