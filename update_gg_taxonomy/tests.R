source("./update_taxonomy_refseq.R")

# TESTS for function is_taxonomy_incomplete
### Function for determining if original greengenes taxonomy is incomplete.
# Incomplete is not assigned at specified taxonomic level.
# Returns boolean TRUE if taxonomy is incomplete.

current_string <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__buchneri"
is_taxonomy_incomplete(current_string, "spcs")
is_taxonomy_incomplete(current_string, "gen")
is_taxonomy_incomplete(current_string, "fam")

current_string2 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__"
is_taxonomy_incomplete(current_string2, "spcs")
is_taxonomy_incomplete(current_string2, "gen")
is_taxonomy_incomplete(current_string2, "fam")

current_string3 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__;s__"
is_taxonomy_incomplete(current_string3, "spcs")
is_taxonomy_incomplete(current_string3, "gen")
is_taxonomy_incomplete(current_string3, "fam")

current_string4 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae"
is_taxonomy_incomplete(current_string4, "spcs")
is_taxonomy_incomplete(current_string4, "gen")
is_taxonomy_incomplete(current_string4, "fam")

current_string5 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;"
current_string5 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__;g__;s__"
current_string5 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales"
is_taxonomy_incomplete(current_string5, "spcs")
is_taxonomy_incomplete(current_string5, "gen")
is_taxonomy_incomplete(current_string5, "fam")





# Reading files for next tests

# fasta file
data_fasta <- readDNAStringSet("C:/Users/marce/OneDrive/[9] Doc/[6] Papers/Paper metagnom/2_diversidad/8_dna_sequences.fasta")
# Reading microbial data base for tests
microbial_database <- blast(db = "./16SMicrobialDB/16SMicrobial")
# set NCBI's entrez api key
Sys.setenv(ENTREZ_KEY = "ed4870836e8f61529227d9176a7c4a994c07")
# Check that key variable is in path.
getkey(service = "entrez")
# Check if blast is in PATH.
Sys.which("blastn")

head(data_fasta)

fasta1 <- data_fasta[1,]
fasta1
fasta2 <- data_fasta[10,]
fasta2
fasta3 <- data_fasta[30,]
fasta3

# Converting fasta file to text.
#as.character(fasta1[[1]])
# Searching fasta file by otu_id
#data_fasta[data_fasta@ranges@NAMES == "4fdb872697ff4712d1408c2a31c881ef"]



# TESTS blast_n_get_ncbi_tax

ref_seq_taxonomy1 <- blast_n_get_ncbi_tax(fasta1, microbial_database = microbial_database)
ref_seq_taxonomy1

ref_seq_taxonomy2 <- blast_n_get_ncbi_tax(fasta2, microbial_database = microbial_database)
ref_seq_taxonomy2

ref_seq_taxonomy3 <- blast_n_get_ncbi_tax(fasta3, microbial_database = microbial_database)
ref_seq_taxonomy3




# TESTS parse_ncbi_to_gg

parse_ncbi_to_gg(ref_seq_taxonomy1[4])
parse_ncbi_to_gg(ref_seq_taxonomy2[4])
parse_ncbi_to_gg(ref_seq_taxonomy3[4])


# TESTS update_taxonomy_refseq

##########################################################################################

data_fasta <- readDNAStringSet("C:/Users/marce/OneDrive/[9] Doc/[6] Papers/Paper metagnom/2_diversidad/8_dna_sequences.fasta")

# set NCBI's entrez api key
Sys.setenv(ENTREZ_KEY = "ed4870836e8f61529227d9176a7c4a994c07")

# Check that key variable is in path.
getkey(service = "entrez")

# Check if blast is in PATH.
Sys.which("blastn")





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

new_tax3 <- apply(new_tax2, 2, as.character)

new_tax3

write.table(new_tax3, file = "./taxonomy_updated_picrust_13_8.tsv", sep = "\t", row.names = FALSE)


