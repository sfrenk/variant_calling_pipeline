suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("input", help = "jitterbug file to process")
parser$add_argument("-o", "--output", help = "output filename", default = "jitterbug_genotypes.txt")

args <- parser$parse_args()

# Read in transposon insertion data
transposons <- read.table(args$input, sep = "\t", stringsAsFactors = FALSE, col.names = c("chr", "start", "end", "type", "sample"))

# Create empty dataframe to be populated with genotyping data
unique_transposons <- data.frame("chr" = character(), "start" = numeric(), "end" = numeric(), "type" = character(), "idx" = numeric(), stringsAsFactors = FALSE)

for (i in unique(transposons$sample)){
    unique_transposons[i] <- numeric()
}

################ Genotyping ##################
print("Genotyping...")

# For each insertion, add the insertion to "unique_transposons" if it is not already listed. If it's alread there, update the sample genotype.
idx <- 1
for (i in 1:nrow(transposons)) {
    
    name <- transposons[i, "sample"]
        
    match <- unique_transposons %>% filter(type == transposons[i, "type"], chr == transposons[i, "chr"], abs(start - transposons[i, "start"]) < 100, abs(end - transposons[i, "end"]) < 100)
        
    if (nrow(match) == 0) {
        
        # If the insertion is new, create a new entry for the genotype table
        new_entry <- transposons[i,] %>% select(-sample)
        new_entry$idx = idx
        idx <- idx + 1
        for (j in unique(transposons$sample)){
            new_entry[j] <- 0
        }
        new_entry[name] <- 1
        unique_transposons <- rbind(unique_transposons, new_entry)
        
    } else {
        # If variant is already present, update the genotype for the sample.
        unique_transposons[unique_transposons$idx == match[1, "idx"], name] <- 1
    }
    
}


print("Writing output...")
unique_transposons <- unique_transposons[complete.cases(unique_transposons),]

write.table(unique_transposons, args$output, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
