
fungi_database <- blast(db = "./fungi_analysis/fungi.ITS.fna/fungi.ITS.fna")


fungi_fasta <- readDNAStringSet("./fungi_analysis/dna-sequences.fasta")

fungi_old_tax <- read.table("./data_for_tests/fungi_analysis/taxonomy.tsv", sep = "\t", header = FALSE)


new_tax <- fungi_old_tax
new_tax$V2 <- as.character(new_tax$V2)
head(new_tax)



is_taxonomy_incomplete <- function(taxonomy) {
    str_splt <- strsplit(taxonomy, ";")
    # Species level
    return(length(str_splt[[1]]) < 7 || length(strsplit(str_splt[[1]][7], "__")[[1]]) < 2 || strsplit(str_splt[[1]][7], "__")[[1]][1] == "unidentified")
    # Genus level
    #return(length(str_splt[[1]]) < 6 || length(strsplit(str_splt[[1]][6], "__")[[1]]) < 2)
    # Family
    #return(length(str_splt[[1]]) < 5 || length(strsplit(str_splt[[1]][5], "__")[[1]]) < 2)
}

blast_n_get_ncbi_tax <- function(seq, perc_ident = 97) {
    blast_result <- arrange(predict(fungi_database, seq), desc(Bits), desc(Perc.Ident))[1,]
    return(c((blast_result["Perc.Ident"] >= perc_ident), (blast_result["Perc.Ident"]), (classification(genbank2uid(id = blast_result["SubjectID"][1, 1]), db = "ncbi"))))
}


#fungi_fasta = fungi_fasta[1530:1533,]

head(fungi_fasta)



for (i in 1:length(fungi_fasta)) {
    print(sprintf("OTU: %s / %s", i, length(fungi_fasta)))
    # Get current sequence info from fasta file.
    current_sequence <- fungi_fasta[i,]
    # Get the OTU id from fasta file.
    seq_id <- current_sequence@ranges@NAMES[1]
    # Get the current taxonomy for current OTU
    current_seq_taxonomy <- as.character((new_tax %>% filter(V1 == seq_id) %>% select(V2))[1, 1])
    # If current taxonomy is incomplete, update
    if (is_taxonomy_incomplete(current_seq_taxonomy)) {
        tryCatch({
            # Blast and get new taxonomy from ncbi
            new_ncbi_taxonomy = blast_n_get_ncbi_tax(current_sequence, 97)
            # If ident. perc is above specified.
                if (new_ncbi_taxonomy[[1]]) {
                    # Replace taxonomy and ident. perc.
                    new_tax <- new_tax %>% mutate(V3 = replace(V3, which(V1 == seq_id), (new_ncbi_taxonomy[2])))
                    new_tax <- new_tax %>% mutate(V2 = replace(V2, which(V1 == seq_id), parse_ncbi_to_gg(new_ncbi_taxonomy[3])))
                }
            },
                error = function(e) {
                    #message("Here's the original error message:\n")
                    #message(e)
                },
                warning = function(e) {
                    #message("Here's the original error message:\n")
                    #message(e)
                }
        )
    }
    # Wait 0.2 seconds to prevent ncbi's server to stop the process.
    Sys.sleep(0.1)
}

new_tax3 <- apply(new_tax, 2, as.character)

write.table(new_tax3, file = "./taxonomy_updated__fungi_species.tsv", sep = "\t", row.names = FALSE)
