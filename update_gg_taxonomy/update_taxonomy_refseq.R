# Install and load required packages

# Bioconductor Packages
source("http://bioconductor.org/biocLite.R")
if (!"Biostrings" %in% installed.packages()) biocLite("Biostrings")
library("Biostrings")
# CRAN packages

if (!"dplyr" %in% installed.packages()) install.packages("dplyr")
library("dplyr")
if (!"taxize" %in% installed.packages()) install.packages("taxize")
library("taxize")
if (!"devtools" %in% installed.packages()) install.packages("devtools")
library("devtools")
if (!"rBLAST" %in% installed.packages()) install_github("mhahsler/rBLAST")
library("rBLAST")





### Function for determining if original greengenes taxonomy is incomplete.
# Incomplete is not assigned at specified taxonomic level.
# Returns boolean TRUE if taxonomy is incomplete.

is_taxonomy_incomplete <- function(taxonomy, level = spcs) {
    str_splt <- strsplit(taxonomy, ";")
    switch(level,
            spcs = {
            # Species level
            return(length(str_splt[[1]]) < 7 || length(strsplit(str_splt[[1]][7], "__")[[1]]) < 2)
            },
            gen = {
            # Genus level
            return(length(str_splt[[1]]) < 6 || length(strsplit(str_splt[[1]][6], "__")[[1]]) < 2)
            },
            fam = {
            # Family
            return(length(str_splt[[1]]) < 5 || length(strsplit(str_splt[[1]][5], "__")[[1]]) < 2)
            })
}





### Function for parsing the NCBI's taxonomy into greengenes format.
# Returns string of taxonomy in greengenes format.
# What is NCBI's format?

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





### Function for blasting sequences against microbial refseq 16s database and returning the taxonomy from NCBI's taxonomy database.
# Orders results by bits and percent of identity. Then chooses the first match with identity percent.
# Returns vector of size 3 including:
# 1) bool if % indentity is at least given number, def. 97 and bits >= 100;
# 2) float of blast result's % identity;
# 3) taxonomy in ncbi's format.

blast_n_get_ncbi_tax <- function(seq, perc_ident = 97, min_bits = 100) {
    # Blasts sequence and return a table of results. Orders it by BITS from greater value.
    blast_result_table <- arrange(predict(microbial_database, seq), desc(Bits))
    # variable to store if we found a sequence with at least percent id value.
    min_percent_found <- FALSE
    # For each row of the results table.
    for (row in 1:nrow(blast_result_table)) {
        # if current result has a percent id of at least given value.
        if (blast_result_table[row, ]$Perc.Ident >= perc_ident) {
            # set to true as we found a seq with at least given percent of id.
            min_percent_found <- TRUE
            # make current row the blast result.
            blast_result <- blast_result_table[row,]
            break
        }
    }
    # if we didnt found a seq with at least given perc. id.
    if (min_percent_found == FALSE) {
        # result with higher bits value is chosen.
        blast_result <- blast_result_table[1,]
    }
    # Return our results vector. 
    return(c(
           ((blast_result["Perc.Ident"] >= perc_ident) & (blast_result["Bits"] >= min_bits)),
           (blast_result["Perc.Ident"]),
           (blast_result["Bits"]),
           # Grab taxonomy from ncbi taxonomy server
           (classification(genbank2uid(id = blast_result["SubjectID"][1, 1]),db = "ncbi"))))
}











# Check if blast is in PATH.

blast_path <- Sys.which("blastn")

blast_path

as.character(blast_path)[1]

length(blast_path)

if (length(blast_path) > 1) {
    print("Blast found at:")
    print(blast_path)
}else {
    "Blast not found! Install Blast"
}



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
#tax_table <- read.table("C:/Users/marce/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/chelo/hidro_agave_diversidad/2_resultados/9_tsv_gg_sk/taxonomy.tsv", sep = "\t", header = TRUE)
tax_table <- read.table("C:/Users/marce/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/chelo/picrust_analyses_13_8/17_picrust_tsv_gg_13_8/taxonomy.tsv", sep = "\t", header = TRUE)
tax_table
# Load sequences to be analyzed. Ids corresponding to taxonomy OTUs.
data_fasta <- readDNAStringSet("C:/Users/marce/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/chelo/picrust_analyses/14_closedRef_forPICRUSt_13_8/dna-sequences.fasta")
# Check fasta sequences
data_fasta

# This is the copy of taxonomy file that we are modifying.
new_tax2 <- tax_table
new_tax2$Taxon <- as.character(new_tax2$Taxon)
new_tax2











############ MAIN ############


update_taxonomy_refseq <- function(current_seq_taxonomy, data_fasta, entrez, level = 97, phyl_group = "bacteria") {

}

for (i in 1:length(data_fasta)) {
    print(sprintf("OTU: %s / %s", i, length(data_fasta)))
    # Get current sequence info from fasta file.
    current_sequence <- data_fasta[i,]
    # Get the OTU id from fasta file.
    seq_id <- current_sequence@ranges@NAMES[1]
    # Get the current taxonomy for current OTU
    current_seq_taxonomy <- as.character((new_tax2 %>% filter(Feature.ID == seq_id) %>% select(Taxon))[1, 1])
    # If current taxonomy is incomplete, update
    if (is_taxonomy_incomplete(current_seq_taxonomy)) {
        # Blast and get new taxonomy from ncbi
        new_ncbi_taxonomy = blast_n_get_ncbi_tax(current_sequence, 97)
        # If ident. perc is above specified.
        if (new_ncbi_taxonomy[[1]]) {
            # Replace taxonomy and ident. perc.
            new_tax2 <- new_tax2 %>% mutate(Confidence = replace(Confidence, which(Feature.ID == seq_id), (new_ncbi_taxonomy[2])))
            new_tax2 <- new_tax2 %>% mutate(Taxon = replace(Taxon, which(Feature.ID == seq_id), parse_ncbi_to_gg(new_ncbi_taxonomy[4])))
        }
    }
    # Wait 0.2 seconds to prevent ncbi's server to stop the process.
    Sys.sleep(0.2)
}

new_tax3 <- apply(new_tax2, 2, as.character)

new_tax3

write.table(new_tax3, file = "./taxonomy_updated_picrust_13_8.tsv", sep = "\t", row.names = FALSE)

# To do:
# Tests for individual functions.
# Turn main into a function.
# Add flexibility for choosing the taxonomic level of analyses.
# if you want to check all