############ Install required packages ############

source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
install.packages("devtools")
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
# Load sequences to be analyzed. Ids corresponding to taxonomy OTUs.
data_fasta <- readDNAStringSet("C:/Users/marce/Desktop/dna-sequences.fasta")
# Check fasta sequences
data_fasta




############ Function for parsing the NCBI's taxonomy into greegngenes format ############

parse_ncbi_to_gg <- function(ncbi_tax) {
    return(sprintf("k__%s; p__%s; c__%s; o__%s; f__%s; g__%s; s__%s",
           # Kingdom
           ncbi_tax[[1]][2, 1],
           # Phylum
           ncbi_tax[[1]][4, 1],
           # Class
           ncbi_tax[[1]][5, 1],
           # Order
           ncbi_tax[[1]][6, 1],
           # Family
           ncbi_tax[[1]][7, 1],
           # Genera
           ncbi_tax[[1]][8, 1],
           # Species
           strsplit(tax_classification[[1]][9, 1], " ")[[1]][2]))
}

############ Function for blasting sequences, retrieving taxonomy from NCBI taxonomy database and replacing an old taxonomy with the NCBI's taxonomy ############

blast_and_replace <- function(seq) {
    blast_result <- arrange(predict(microbial_database, seq), desc(Perc.Ident), desc(Bits), E)[1,]["SubjectID"]
    tax_classification <- classification(genbank2uid(id = blast_result[1, 1]), db = "ncbi")
    parse_ncbi_to_gg(tax_classification)
}

blast_and_replace(data_fasta[1,])






for (i in 1:length(data_fasta)) {
    print(i)
    print(data_fasta[i,])
    #
}







tax_id <- genbank2uid(id = blast_result[1,]["SubjectID"][1, 1])

tax_id[[1]][1]

tax_classification <- classification(tax_id[[1]][1], db = "ncbi")
tax_classification

tax_classification[[1]][2, 1]
tax_classification[[1]][4, 1]
tax_classification[[1]][5, 1]
tax_classification[[1]][6, 1]
tax_classification[[1]][7, 1]
tax_classification[[1]][8, 1]
tax_classification[[1]][9, 1]

typeof(tax_classification[[1]][9, 1])

strsplit(tax_classification[[1]][9, 1], " ")[[1]][2]


final_string <- sprintf("k__%s; p__%s; c__%s; o__%s; f__%s; g__%s; s__%s", tax_classification[[1]][2, 1], tax_classification[[1]][4, 1],
                        tax_classification[[1]][5, 1], tax_classification[[1]][6, 1], tax_classification[[1]][7, 1], tax_classification[[1]][8, 1],
                        strsplit(tax_classification[[1]][9, 1], " ")[[1]][2])

final_string