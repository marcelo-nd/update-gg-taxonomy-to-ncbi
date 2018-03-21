fungi_database <- blast(db = "./fungi_db/fungi_db")

fungi_fasta <- readDNAStringSet("C:/Users/marce/Desktop/hongos/dna-sequences.fasta")

fungi_old_tax <- read.table("C:/Users/marce/Desktop/hongos/taxonomy.tsv", sep = "\t", header = FALSE)

head(fungi_fasta)

fungi_fasta[1,]

head(predict(fungi_database, fungi_fasta[1,]))

fungi_fasta_test <- head(fungi_fasta, 5)

fungi_tax_test <- head(fungi_old_tax, 5)

fungi_tax_test

fungi_new_tax2 <- fungi_tax_test
fungi_new_tax2$V2 <- as.character(fungi_new_tax2$V2)
fungi_new_tax2

data_fasta <- fungi_fasta_test

data_fasta

new_tax2 <- fungi_new_tax2
new_tax2

blast_n_get_ncbi_tax <- function(seq, perc_ident = 97) {
    blast_result <- arrange(predict(fungi_database, seq), desc(Bits), desc(Perc.Ident))[1,]
    print(blast_result)
    return(c((blast_result["Perc.Ident"] >= perc_ident), (blast_result["Perc.Ident"]), (classification(genbank2uid(id = blast_result["SubjectID"][1, 1]), db = "ncbi"))))
}

for (i in 1:length(data_fasta)) {
    print(sprintf("OTU: %s / %s", i, length(data_fasta)))
    # Get current sequence info from fasta file.
    current_sequence <- data_fasta[i,]
    # Get the OTU id from fasta file.
    seq_id <- current_sequence@ranges@NAMES[1]
    # Get the current taxonomy for current OTU
    current_seq_taxonomy <- as.character((new_tax2 %>% filter(V1 == seq_id) %>% select(V2))[1, 1])
    # If current taxonomy is incomplete, update
    if (is_taxonomy_incomplete(current_seq_taxonomy)) {
        # Blast and get new taxonomy from ncbi
        new_ncbi_taxonomy = blast_n_get_ncbi_tax(current_sequence, 90)
        # If ident. perc is above specified.
        if (new_ncbi_taxonomy[[1]]) {
            # Replace taxonomy and ident. perc.
            new_tax2 <- new_tax2 %>% mutate(V3 = replace(V3, which(V1 == seq_id), (new_ncbi_taxonomy[2])))
            new_tax2 <- new_tax2 %>% mutate(V2 = replace(V2, which(V1 == seq_id), parse_ncbi_to_gg(new_ncbi_taxonomy[3])))
        }
    }
    # Wait 0.2 seconds to prevent ncbi's server to stop the process.
    Sys.sleep(0.2)
}

new_tax2