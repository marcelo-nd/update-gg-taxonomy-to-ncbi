# TESTS
#current_string <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__buchneri"
#is_taxonomy_incomplete(current_string, "spcs")
#is_taxonomy_incomplete(current_string, "gen")
#is_taxonomy_incomplete(current_string, "fam")

#current_string2 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__"
#is_taxonomy_incomplete(current_string2, "spcs")
#is_taxonomy_incomplete(current_string2, "gen")
#is_taxonomy_incomplete(current_string2, "fam")

#current_string3 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__;s__"
#is_taxonomy_incomplete(current_string3, "spcs")
#is_taxonomy_incomplete(current_string3, "gen")
#is_taxonomy_incomplete(current_string3, "fam")

#current_string4 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae"
#is_taxonomy_incomplete(current_string4, "spcs")
#is_taxonomy_incomplete(current_string4, "gen")
#is_taxonomy_incomplete(current_string4, "fam")

#current_string5 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;"
#current_string5 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__;g__;s__"
#current_string5 <- "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales"
#is_taxonomy_incomplete(current_string5, "spcs")
#is_taxonomy_incomplete(current_string5, "gen")
#is_taxonomy_incomplete(current_string5, "fam")



# TESTS parse
# faltan

# TESTS blast



##########################################################################################


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


# Test for manipulating fasta file sequences. (TEMP)
data_fasta <- readDNAStringSet("F:/1_Programs/OneDrive/[9] Doc/[6] Papers/Paper metagnom/2_diversidad/8_dna_sequences.fasta")

head(data_fasta)

fasta1 <- data_fasta[1,]

as.character(fasta1[[1]])

data_fasta[data_fasta@ranges@NAMES == "4fdb872697ff4712d1408c2a31c881ef"]