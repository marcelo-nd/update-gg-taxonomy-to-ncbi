############ Install required packages ############

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("devtools")
install.packages("dplyr")
install.packages("taxize")
library("devtools")
install_github("mhahsler/rBLAST")




############ Load required packages ############

library("Biostrings")
library("rBLAST")
# Check if blast is in PATH.
Sys.which("blastn")
library("dplyr")
library("taxize")
# set NCBI's entrez api key
Sys.setenv(ENTREZ_KEY = "ed4870836e8f61529227d9176a7c4a994c07")
# Check that key variable is in path.
getkey(service = "entrez")





############ Load required data ############

# Extract reference database if not done already.
download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz", "16SMicrobial.tar.gz", mode = 'wb')
untar("16SMicrobial.tar.gz", exdir = "16SMicrobialDB")
# Load reference database.
microbial_database <- blast(db = "./16SMicrobialDB/16SMicrobial")
# Check database
microbial_database
# Load taxonomy table to be updated with NCBI's taxonomy.
tax_table <- read.table("C:/Users/marce/OneDrive/update-tax/taxonomy.tsv", sep = "\t", header = TRUE)
tax_table
# Load sequences to be analyzed. Ids corresponding to taxonomy OTUs.
data_fasta <- readDNAStringSet("C:/Users/marce/OneDrive/update-tax/dna-sequences.fasta")
# Check fasta sequences
data_fasta

data_small <- head(data_fasta, 3)
data_small

tax2 <- head(tax_table, 3)
tax2

new_tax2 <- tax2
new_tax2$Taxon <- as.character(new_tax2$Taxon)
new_tax2

############ Function for determining if original greengenes taxonomy is incomplete. Returns boolean TRUE if taxonomy is incomplete. ############

is_taxonomy_incomplete <- function(taxonomy) {
    return(TRUE)
}

############ Function for parsing the NCBI's taxonomy into greegngenes format. Returns string of taxonomy in greengenes format. ############

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


############ Function for blasting sequences agains microbial refseq 16s database and returning the taxonomy from NCBI's taxonomy database. ############
############ Returns vector of size 3 including: bool if % indentity is at least given number, def. 97; float of blast result's % identity; taxonomy in ncbi's format ############

blast_n_get_ncbi_tax <- function(seq, perc_ident = 97) {
    blast_result <- arrange(predict(microbial_database, seq), desc(Bits), desc(Perc.Ident))[1,]
    return(c((blast_result["Perc.Ident"] >= perc_ident), (blast_result["Perc.Ident"]), (classification(genbank2uid(id = blast_result["SubjectID"][1, 1]), db = "ncbi"))))
}

for (i in 1:length(data_small)) {
    current_sequence <- data_small[i,]
    seq_id <- current_sequence@ranges@NAMES[1]
    current_seq_taxonomy <- as.character((new_tax2 %>% filter(Feature.ID == seq_id) %>% select(Taxon))[1, 1])
    if (is_taxonomy_incomplete(current_seq_taxonomy)) {
        new_ncbi_taxonomy = blast_n_get_ncbi_tax(current_sequence, 97)
        print(new_ncbi_taxonomy)

        new_tax2 <- new_tax2 %>% mutate(Taxon = replace(Taxon, which(Feature.ID == seq_id), parse_ncbi_to_gg(new_ncbi_taxonomy[3])))

        print("add")
        new_tax2 <- new_tax2 %>% mutate(Confidence = replace(Confidence, which(Feature.ID == seq_id), (new_ncbi_taxonomy[2])))
    }  
}