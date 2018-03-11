suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("input", help = "meerkat .variants files to process", nargs = "+")
parser$add_argument("-o", "--output", help = "output filename", default = "meerkat_structural_variants.txt")

args <- parser$parse_args()

# Read in unfiltered transposon insertion data
#setwd("/home/sfrenk/Documents/Projects_sup/braxton_paper/meerkat")
#file_list <-  c("TRT1_2.variants", "EXO1_A.variants")

# Compile df
print("Reading and compiling data...")

vars <- data.frame("type" = character(), "mechanism" = character(), "cluster_id" = character(), "sup_pairs" = character(), "split_reads" = character(), "chr" = character(), "start" = numeric(), "end" = numeric(), "size" = character(), "opt1" = character(), "opt2" = character(), "opt3" = character(), "opt4" = character(), "opt5" = character(), "opt6" = character(), "opt7" = character(), stringsAsFactors = FALSE, sample = character())

classes <- sapply(vars, class)

for (i in args$input){
    data <- read.table(i, stringsAsFactors = FALSE, sep = "\t", fill = TRUE, col.names = colnames(vars), colClasses = classes)
    data$sample <- i
    vars <- rbind(vars, data)
}

vars[is.na(vars)] <- "empty"
vars[vars==""] <- "empty"

################ Genotyping ##################
print("Genotyping...")

# Initialize data frame to store variants
unique_vars <- data.frame("type" = character(), "mechanism" = character(), "cluster_id" = numeric(), "sup_pairs" = numeric(), "split_reads" = numeric(), "chr" = character(), "start" = numeric(), "end" = numeric(), "size" = numeric(), "opt1" = character(), "opt2" = character(), "opt3" = character(), "opt4" = character(), "opt5" = character(), "opt6" = character(), "opt7" = character(), "idx" = numeric(), stringsAsFactors = FALSE)

# Add genotype column for each sample
for (i in args$input){
    unique_vars[i] <- numeric()
}

# For each variant, check add the variant to "unique_vars" if it is not already listed. If it's already in "unique_var", update the sample genotype.
idx <- 1
for (i in 1:nrow(vars)) {
    
    name <- vars[i, "sample"]
    match <- unique_vars %>% filter(type == vars[i, "type"], chr == vars[i, "chr"], abs(start - vars[i, "start"]) < 100, abs(end - vars[i, "end"]) < 100)

    if (nrow(match) == 0) {
        
        # If the variant is new, create a new entry for the genotype table
        new_entry <- vars[i,] %>% select(-sample)
        new_entry$idx = idx
        idx <- idx + 1
        for (j in args$input){
            new_entry[j] <- 0
        }
        new_entry[name] <- 1
        unique_vars <- rbind(unique_vars, new_entry)
        
    } else {
        # If variant is already present, update the genotype for the sample.
        unique_vars[unique_vars$idx == match[1, "idx"], name] <- 1
    }
    
}

unique_vars <- unique_vars[!is.na(unique_vars),]

write.table(unique_vars, args$output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
