# Install and load required packages

# Bioconductor Packages

if (!"Biostrings" %in% installed.packages()) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
}

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

blast_n_get_ncbi_tax <- function(seq, perc_ident = 97, min_bits = 100, microbial_database) {
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





#### 


update_taxonomy_refseq <- function(taxonomy_table, data_fasta, level = spcs, phyl_group = "bacteria", update_all = FALSE, microbial_database_exists = FALSE) {
    # selecting working percent identity 
    switch(level,
            spcs = {
    # Species level
    percent = 97
    },
            gen = {
    # Genus level
    percent = 95
    },
            fam = {
    # Family
    percent = 90
    })

    # Iterate over taxonomy table
    for (otu_entry in taxonomy_table) {
        # grab taxonomy for each entry. Is a string?!!! ######## CHECK ######
        current_taxonomy <- as.character(taxonomy_table[otu_entry,] %>% select(Taxon))[1, 1])
        # If update_all TRUE or if not If current taxonomy is incomplete.
        if (update_all | is_taxonomy_incomplete(current_seq_taxonomy)) {
            # get otu_id ######## CHECK ######
            current_id <- taxonomy_table[otu_entry,] %>% select(Feature.ID))[1, 1])
            # get otu sequence from fasta file. ######## CHECK ######
            current_sequence <- data_fasta[data_fasta@ranges@NAMES == current_id]
            # get new taxonomy
            new_ncbi_taxonomy <- blast_n_get_ncbi_tax(seq = current_sequence, perc_ident = percent, min_bits = 100, microbial_database)
            # If ident. perc is above specified.
            if (new_ncbi_taxonomy[[1]]) {
                # Replace taxonomy and ident. perc.
                taxonomy_table <- taxonomy_table %>% mutate(Confidence = replace(Confidence, which(Feature.ID == seq_id), (new_ncbi_taxonomy[2])))
                taxonomy_table <- taxonomy_table %>% mutate(Taxon = replace(Taxon, which(Feature.ID == seq_id), parse_ncbi_to_gg(new_ncbi_taxonomy[4])))
            }
        }
        # Wait 0.2 seconds to prevent ncbi's server to explode.
        Sys.sleep(0.2)
    }
    return(taxonomy_table)
}



# To do:
# Tests for individual functions.
# fucntion for downldng microbial database