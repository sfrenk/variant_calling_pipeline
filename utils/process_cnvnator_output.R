#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

# Compiles multiple CNVnator output files into a single VCF-style file

parser <- ArgumentParser()
parser$add_argument("input", help = "CNVnator output files to process", nargs = "+")
parser$add_argument("-o", "--output", help = "output filename", default = "cnvnator_variants.txt")

args <- parser$parse_args()

# Read in unfiltered transposon insertion data

# Compile df
print("Reading and compiling data...")

vars <- data.frame("sample" = character(), "type" = character(), "chrom" = character(), "start" = numeric(), "end" = numeric(), "size" = numeric(), "normalized_RD" = numeric(), "eval1" = numeric(), "eval2" = numeric(), "eval3" = numeric(), "eval4" = numeric(), "q0" = numeric(), stringsAsFactors = FALSE)

classes <- sapply(vars, class)

for (i in args$input){
    data <- read.table(i, stringsAsFactors = FALSE, sep = "\t", col.names = c("type", "coords", "size", "normalized_RD", "eval1", "eval2", "eval3", "eval4", "q0"), fill = TRUE)
    data$chrom <- sapply(data$coords, function(x) str_match(x, "(.*):")[,2])
    data$start <- sapply(data$coords, function(x) as.numeric(str_match(x, ":([0-9]+)")[,2]))
    data$end<- sapply(data$coords, function(x) as.numeric(str_match(x, "-([0-9]+)")[,2]))
    data$sample <- i
    data <- data[, colnames(vars)]
    vars <- rbind(vars, data)
}


################ Genotyping ##################
print("Genotyping...")

# Add idx column to vars. This will be used to extract unique CNVs at the end
vars$idx <- 0

# Initialize data frame to store variants
unique_vars <- data.frame("sample" = character(), "type" = character(), "chrom" = character(), "start" = numeric(), "end" = numeric(), "size" = numeric(), "normalized_RD" = numeric(), "eval1" = numeric(), "eval2" = numeric(), "eval3" = numeric(), "eval4" = numeric(), "q0" = numeric(), stringsAsFactors = FALSE)


# Add genotype column for each sample
for (i in args$input){
    unique_vars[i] <- numeric()
}

# For each variant, check add the variant to "unique_vars" if it is not already listed. If it's already in "unique_var", update the sample genotype.
idx <- 1
for (i in 1:nrow(vars)) {
    
    name <- vars[i, "sample"]
    
    # Check if variant is already in the unique_vars database
    match <- unique_vars %>% filter(type == vars[i, "type"], chrom == vars[i, "chrom"], abs(start - vars[i, "start"]) < 100, abs(end - vars[i, "end"]) < 100)
        
    if (nrow(match) == 0) {
        
        # If the variant is new, create a new entry for the genotype table
        new_entry <- vars[i,] %>% select(-sample)
        new_entry$idx <- idx
        vars[i, "idx"] <- idx
        idx <- idx + 1
        for (j in args$input){
            new_entry[j] <- 0
        }
        new_entry[name] <- 1
        unique_vars <- rbind(unique_vars, new_entry)
        
    } else {
        # If variant is already present, update the genotype for the sample.
        unique_vars[unique_vars$idx == match[1, "idx"], name] <- 1
        
        # Expand coordinates to cover range of new entry
        unique_vars[unique_vars$idx == match[1, "idx"], "start"] <- min(c(vars[i, "start"]), unique_vars[unique_vars$idx == match[1, "idx"], "start"])
        unique_vars[unique_vars$idx == match[1, "idx"], "end"] <- max(c(vars[i, "end"]), unique_vars[unique_vars$idx == match[1, "idx"], "end"])
       
    }
    
}


print("Writing output...")
unique_vars <- unique_vars[complete.cases(unique_vars),]

write.table(unique_vars, args$output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Get events that are only present in a single sample
unique_vars$nsamp <- rowSums(unique_vars[,args$input])
unique_idx <- unique_vars %>% filter(nsamp == 1) %>% .$idx
unique_vars_unique <- vars %>% filter(idx %in% unique_idx)
write.table(unique_vars_unique, paste0(args$output, ".unique.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
