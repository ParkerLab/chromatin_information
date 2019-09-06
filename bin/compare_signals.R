options(stringsAsFactors = F)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-m", "--motif"), type="character", default=NULL, 
                help="motif name", metavar="character"),
    make_option(c("-d", "--dir"), type="character", default=NULL, 
                help="Input directory", metavar="character"),
    make_option(c("-o", "--out"), type="character", default=NULL, 
                help="Output directory", metavar="character"),
    make_option(c("-n", "--ntiles"), type="integer", default=20, 
                help="Number of tiles to split signals", metavar="integer")
)
opt_parser <- OptionParser(option_list=option_list);
args <- parse_args(opt_parser)



write("User parameters:", stdout())
write(sprintf("Input dir: %s", args$dir), stdout())
write(sprintf("Output dir: %s", args$out), stdout())
write(sprintf("File handle: %s\n", args$motif), stdout())


## Helper functions

# Read all data for a given motif
read_motif <- function(handle, ntiles = args$ntiles){
    colnames <- c("chrom", "start", "end", "motif", "pwm", "strand", "atac_signal", "tf_signal")
    
    # Create dummy output in case files are empty and don't crash snakemake
    if(file.size(in_plus) == 0 | file.size(in_plus) == 0){
        write("One or more input files are empty. Creating dummy output files and exiting.", stdout())
        write.table(data.frame(), out_plus)
        write.table(data.frame(), out_minus)
        quit(save = 'no', status = 0)
    }
    
    plus  <- read.table(in_plus, col.names = colnames) %>% 
        mutate(cohesin = "plus")
    minus <- read.table(in_minus, col.names = colnames) %>% 
        mutate(cohesin = "minus")
    
    out <- bind_rows(plus, minus)
    out <- out %>% 
        mutate(info = handle) %>% 
        separate(info, into = c("factor", "motif2", "experiment", "rad21_exp"), sep = "__") %>% 
        select(-motif2) %>% 
        mutate(index = 1:nrow(out)) %>% 
        mutate(chip_quantile = ntile(tf_signal, ntiles),
               atac_quantile = ntile(atac_signal, ntiles),
               pwm_quantile = ntile(pwm, ntiles))
    return(out)
}

# Subsetting function
equal_split <- function(x){
    counts <- min(as.numeric(table(x$cohesin)))
    x_split <- split(x, x$cohesin)
    x_split <- lapply(x_split, sample_n, counts, replace = F)
    x <- bind_rows(x_split)
    return(x)
}

## RUN

indir  <- args$dir
outdir <- args$out
motif  <- args$motif

# # Debug
# indir = "work/intersections_plus_minus/buenrostro_rep1"
# motif = "IKZF2__IKZF2_3__ENCFF489GBB-IKZF2__ENCFF002CHR"
# outdir = "work/comparable_motifs/buenrostro_rep1"

in_plus   <- file.path(indir, sprintf("%s.plus.bed", motif))
in_minus  <- file.path(indir, sprintf("%s.minus.bed", motif))
out_plus  <- file.path(outdir, sprintf("%s.plus.bed", motif))
out_minus <- file.path(outdir, sprintf("%s.minus.bed", motif))

# Load data
write("Reading data...", stdout())
df <- read_motif(motif)

write("Trimming...", stdout())
# Trim first with FIMO
df_trim_subset <- split(df, with(df, list(rad21_exp, experiment, pwm_quantile)))
df_trim_subset <- bind_rows(lapply(df_trim_subset, equal_split))

# Then with ChIP-seq
df_trim_subset <- split(df_trim_subset, with(df_trim_subset, list(rad21_exp, experiment, chip_quantile)))
df_trim_subset <- bind_rows(lapply(df_trim_subset, equal_split))

# And ATAC-seq signal
df_trim_subset <- split(df_trim_subset, with(df_trim_subset, list(rad21_exp, experiment, atac_quantile)))
df_trim_subset <- bind_rows(lapply(df_trim_subset, equal_split))

# Write output
write("Writing output...", stdout())

plus_df <- df_trim_subset %>% 
    filter(cohesin == "plus") %>% 
    select(chrom:cohesin)

minus_df <- df_trim_subset %>% 
    filter(cohesin == "minus") %>% 
    select(chrom:cohesin)

write.table(plus_df, out_plus, row.names = F, col.names = F, quote = F, sep ="\t")
write.table(minus_df, file = out_minus, row.names = F, col.names = F, quote = F, sep = "\t")

write(sprintf("Wrote %s", out_plus), stdout())
write(sprintf("Wrote %s", out_minus), stdout())
write("Done!", stdout())