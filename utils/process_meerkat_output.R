#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

# Compiles multiple Meerkat output files into a single VCF-style file

parser <- ArgumentParser()
parser$add_argument("input", help = "meerkat .variants files to process", nargs = "+")
parser$add_argument("-o", "--output", help = "output filename", default = "meerkat_structural_variants.txt")

args <- parser$parse_args()

# Read in unfiltered transposon insertion data

# Compile df
print("Reading and compiling data...")

vars <- data.frame("type" = character(), "mechanism" = character(), "cluster_id" = character(), "sup_pairs" = character(), "split_reads" = character(), "chr" = character(), "start" = numeric(), "end" = numeric(), "opt1" = character(), "opt2" = character(), "opt3" = character(), "opt4" = character(), "opt5" = character(), "opt6" = character(), "opt7" = character(), "opt8" = character(), stringsAsFactors = FALSE, sample = character())

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
unique_vars <- data.frame("type" = character(), "mechanism" = character(), "cluster_id" = numeric(), "sup_pairs" = numeric(), "split_reads" = numeric(), "chr" = character(), "start" = numeric(), "end" = numeric(), "opt1" = numeric(), "opt2" = character(), "opt3" = character(), "opt4" = character(), "opt5" = character(), "opt6" = character(), "opt7" = character(), "opt8" = character(), "idx" = numeric(), stringsAsFactors = FALSE)

# Add genotype column for each sample
for (i in args$input){
    unique_vars[i] <- numeric()
}

# For each variant, check add the variant to "unique_vars" if it is not already listed. If it's already in "unique_var", update the sample genotype.
idx <- 1
for (i in 1:nrow(vars)) {
    
    name <- vars[i, "sample"]
    
    # Need to process insertions and interchromosomal translocations differently from other events - need to take into account donor as well as event locus
    if (grepl("ins", vars[i, "type"])){
        
        match <- unique_vars %>% filter(type == vars[i, "type"])
        match <- match %>% filter(chr == vars[i, "chr"], abs(start - vars[i, "start"]) < 100, abs(end - vars[i, "end"]) < 100, opt1 == vars[i, "opt2"], abs(as.numeric(opt2) - as.numeric(vars[i, "opt3"])) < 100, abs(as.numeric(opt3) - as.numeric(vars[i, "opt4"])) < 100)
        
        
    } else if (vars[i, "type"] == "transl_inter"){
        # Interchromsomal translations have second chromosome as opt1 instead of size, so the features are shifted forward by one column

        match <- unique_vars %>% filter(type == vars[i, "type"]) 
        match <- match %>% filter(chr == vars[i, "chr"], abs(start - vars[i, "start"]) < 100, abs(end - vars[i, "end"]) < 100, opt1 == vars[i, "opt1"], abs(as.numeric(opt2) - as.numeric(vars[i, "opt2"])) < 100, abs(as.numeric(opt3) - as.numeric(vars[i, "opt3"])) < 100)
        
    } else{
        
        match <- unique_vars %>% filter(type == vars[i, "type"], chr == vars[i, "chr"], abs(start - vars[i, "start"]) < 100, abs(end - vars[i, "end"]) < 100)
        
    }
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


print("Writing output...")
unique_vars <- unique_vars[complete.cases(unique_vars),]

write.table(unique_vars, args$output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
