options(stringsAsFactors = F)
library(data.table)
library(entropy)
library(ggplot2)
library(optparse)
library(RColorBrewer)
library(scales)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))
source('/home/albanus/scripts/multiplot.R')

library(compiler)
enableJIT(3)

option_list = list(
    # Mandatory
    make_option(c("-f", "--file"), type="character", default=NULL, 
                help="Signal file", metavar="character"),
    make_option(c("-o", '--output'), type="character", default=NULL, 
                help="output file handle", metavar="character"),
    make_option(c("-n", '--name'), type="character", default=NULL, 
                help="Motif name", metavar="integer"),
    
    # Optional
    make_option(c("--step"), type = 'integer', default = 3, 
                help="Overlap window step (default = 3)"),
    make_option(c("--maxfrags"), type = 'integer', default = 500000, 
                help="Maximum number of fragments (default = 500,000)"),
    make_option(c("--pseudo_prob"), type = 'numeric', default = 0, 
                help=paste("Add pseudo-counts with binomial prob X",
                           "to sparse data when calculating information content",
                           "(default = 0)")),
    
    # Plotting parameters
    make_option(c("--noplot"), action = 'store_true', default = F,
                help="Don't generate plot"),
    make_option(c("--split"), action = 'store_true', default = F, 
                help="Split plotting output in different files"),
    make_option(c("--randomplot"), action = 'store_true', default = F,
                help="Also generate the randomized plot"),
    make_option(c("--maxpoints"), type = 'integer', default = 250000, 
                help="Maximum points in v-plot (default = 250,000)"),
    make_option(c("--alpha"), type = 'numeric', default = 0.025, 
                help="Plot alpha (default = 0.025)"),
    make_option(c("--size"), type = 'numeric', default = 0.025, 
                help="Plot point size (default = 0.025)"),
    make_option(c("--ylim"), type = 'numeric', default = NA, 
                help="Y-axis limit value for enrichments (default = automatic)"),
    make_option(c("--ylim2"), type = 'numeric', default = 1000, 
                help="Y-axis limit value for V-plots (default = 1,000)"),
    make_option(c("--xlim"), type = 'numeric', default = 500, 
                help="X-axis limit for V-plots"),
    make_option(c("--xrange"), type = 'numeric', default = 500, 
                help="X-axis range to calculate information")
)

opt_parser <- OptionParser(option_list = option_list)
args <- parse_args(opt_parser)

# setwd('/lab/work/albanus/residenceTime_ver3')
# args$file <- 'tmp/vsignal/lab_gm12878/by_0.05/CTCF_known2.vsignal_allRegions.gz'
# args$file <- 'tmp/vsignal/buenrostro_rep1/by_0.05/AHR_1.vsignal_allRegions.gz'
# args$maxfrags <- 50000
# args$output <- 'test_vplot'
# args$name <- 'test_vplot'

input    <- args$file
outfile  <- paste(args$output, '.out', sep = '')
plot_out <- paste(args$output, '.png', sep = '')
title    <- args$name
plot     <- !args$noplot
split    <- args$split
step     <- args$step
xrange   <- args$xrange

# Define which bins should be required for the information content enrichments
lim1 <- 25        # center
lim2 <- c(50, 70) # sides1

# Minimum fragment size to include in the analyses
lowerFragmentSizeLimit <- 41

# Probability of pseudo-counts added to data
pseudo_prob <- args$pseudo_prob

#
## Main functions
#

read_n_trim <- function(input, minSize = lowerFragmentSizeLimit){
    # Read data
    d <- as.data.frame(fread(input))
    if(nrow(d) == 0){
        stop("Parsed vsignal file is empty. Closing.")
    }
    
    # Fetch number of reads in pChART regions before filtering
    pchart_reads_center <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim1)
    pchart_reads_side1  <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) >= lim2[1] & 
                               abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim2[2])
    pchart_reads_raw <<- pchart_reads_center + pchart_reads_side1
    
    # Trim to our region of interest
    d <- d[d$FragmentSize >= minSize,]
    maxSize  <<- max(d$FragmentSize)
    maxCoord <- max(d$FragmentMidpointDistanceToFeatureMidpoint)
    minCoord <- min(d$FragmentMidpointDistanceToFeatureMidpoint)
    coord_span <<- length(minCoord:maxCoord)
    
    # Fetch number of reads in pChART regions after filtering
    pchart_reads_center <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim1)
    pchart_reads_side1  <- sum(abs(d$FragmentMidpointDistanceToFeatureMidpoint) >= lim2[1] & 
                                   abs(d$FragmentMidpointDistanceToFeatureMidpoint) <= lim2[2])
    pchart_reads_trimmed <<- pchart_reads_center + pchart_reads_side1
    
    return(d)
}

# Randomize data
randomize <- function(d){
    fragdists <- d$FragmentMidpointDistanceToFeatureMidpoint
    fragdists <- sample(fragdists, replace = F)
    d$FragmentMidpointDistanceToFeatureMidpoint <- fragdists
    return(d)
}

preprocess <- function(d){
    # Transform fragment sizes and coordinates in factors to avoid 0-count sizes not being included
    minSize = min(d$FragmentSize)
    maxSize = max(d$FragmentSize)
    coords <- d$FragmentMidpointDistanceToFeatureMidpoint
    minCoord <- min(coords)
    maxCoord <- max(coords)
    d$FragmentSize <- factor(d$FragmentSize, levels = minSize:maxSize)
    d$FragmentMidpointDistanceToFeatureMidpoint <- factor(coords, levels = minCoord:maxCoord)
    
    
    # Make random set
    d2 <- list(all = d,
               random = randomize(d))
    
    # Make counts matrix
    for(i in 1:length(d2)){
        d2[[i]] <- table(d2[[i]]$FragmentMidpointDistanceToFeatureMidpoint, d2[[i]]$FragmentSize)
        d2[[i]] <- as.data.frame.matrix(d2[[i]])
    }
    
    # Output
    return(d2)
}

compress_matrix <- function(dat){
    # Compress matrix in 5(y-axis) and 11(x-axis) bp windows to boost signal
    # comp_mat: rows=10bp position coordinate bins; cols=1bp fragment sizes
    # with 3 bp overlap (default) on the position bins only
    backstep <- 10 - step
    nrow_overlap <<- floor((nrow(dat) - backstep + 1) / step)
    comp_mat <- matrix(nrow = nrow_overlap, ncol = ncol(dat))
    comp_mat[1,] <- as.numeric(apply(dat[1:11,], 2, sum))
    end_previous <- 11; i <- 1
    positions <- 5
    while(end_previous < coord_span){
        i <- i + 1
        start <- end_previous - backstep
        end <- start + 10
        if(end >= coord_span){
            end <- coord_span
        }
        end_previous <- end
        # print(paste(i, start, end, sep = "-"))
        comp_mat[i,] <- as.numeric(apply(dat[start:end,], 2, sum))
        positions[i] <- start + 4
    }
    positions <- positions - xrange
    rownames(comp_mat) <- positions
    # print(paste(i, start, end, sep = "-"))
    # print(nrow_overlap)
    rm(start, end, end_previous, i)
    
    # Add pseudo counts? prob=0 means a vector of 0's is added (no effect)
    comp_mat <- comp_mat + rbinom(length(comp_mat), 1, pseudo_prob)
    
    return(comp_mat)
}

# Normalized entropy function
normentropy <- function(x){
    normentro <- (maxentro - entropy(x)) / maxentro
    return(normentro)
}

## Subsample input fragment dataframe into required number of fragments
subsample <- function(frag_df, nfrags){
    if(!is.null(frag_df$FeatureNumber)){
        if(nrow(frag_df) > nfrags){
            counts <- as.data.frame(table(frag_df$FeatureNumber))
            counts$Var1 <- as.numeric(as.character(counts$Var1))
            counts <- counts[sample(1:nrow(counts), replace = F),]
            counts$cumulative <- cumsum(counts$Freq)
            counts$keep <- counts$cumulative <= nfrags
            to_keep <- counts$Var1[counts$keep]
        } else{
            to_keep <- unique(frag_df$FeatureNumber)
        }
        subsampled_data <- filter(frag_df, FeatureNumber %in% to_keep)
    } else {
        if(nrow(frag_df) > nfrags){
            to_keep <- sample(nrow(frag_df), nfrags, replace = F)
            subsampled_data <- frag_df[to_keep,]
        } else{
            subsampled_data <- frag_df
        }
    }
    return(subsampled_data)
}

#
## Run analyses
#

## Step 1: Read and pre-process data

raw_data   <- read_n_trim(input)
usedReads  <- nrow(raw_data)

# Define maximum entropy
size_range <- range(raw_data$FragmentSize)
maxentro  <- entropy(rep(1, length(size_range[1]:size_range[2])))

subsampled_data <- subsample(raw_data, args$maxfrags)

processed_data     <- preprocess(raw_data)
processed_data_sub <- preprocess(subsampled_data)

## Step 2: Make compressed matrices

top_motifs         <- lapply(processed_data[c('all')], compress_matrix)
top_motifs_sub     <- lapply(processed_data_sub[c('all')], compress_matrix)
shuffled_data      <- compress_matrix(processed_data$random)
shuffled_data_sub  <- compress_matrix(processed_data_sub$random)

## Step 3: Calculate entropies

entro_pipeline <- function(top_significant){
    
    my_shuffled_data <- shuffled_data
    
    t_significant_entro <- list()
    shuffled_entro <- list()
    norm_entro <- list()
    
    # Generate the vector of information contents/enrichment
    t_significant_entro$full <- apply(top_significant, 1, normentropy)
    shuffled_entro$full <- apply(my_shuffled_data, 1, normentropy)
    norm_entro$full <- log2(t_significant_entro$full / shuffled_entro$full)
    
    # Make information content/enrichment data frame
    posinfo <- data.frame(motif          = title,
                          signif_motifs  = t_significant_entro$full,
                          shuffled       = shuffled_entro$full,
                          normalized     = norm_entro$full,
                          pos = as.numeric(rownames(top_significant)))
    
    return(posinfo)
}

for(i in c('all')){
    out_variable  <- paste('entro_data_', i, sep = '')
    out_variable2 <- paste('entro_data_sub_', i, sep = '')
    assign(out_variable, entro_pipeline(top_motifs[[i]]))
    assign(out_variable2, entro_pipeline(top_motifs_sub[[i]]))
}; rm(i, out_variable)

posinfo <- cbind(entro_data_all[,1:3],
                 all_normalized     = entro_data_all$normalized,
                 position           = entro_data_all$pos)
posinfo <- as.data.frame(posinfo)

posinfo_sub <- cbind(entro_data_sub_all[,1:3],
                 all_normalized     = entro_data_sub_all$normalized,
                 position           = entro_data_sub_all$pos)
posinfo_sub <- as.data.frame(posinfo_sub)

# Output information content enrichments
write.table(posinfo,  file = gzfile(paste(args$output, '.posinfo.gz',sep='')), sep = "\t",
            row.names = F, quote = F)

write.table(posinfo_sub,  file = gzfile(paste(args$output, '.posinfo_sub.gz',sep='')), sep = "\t",
            row.names = F, quote = F)

## Step 4: Calculate (non-normalized) pChART score

# Define which bins should be required for each track
vic_cols <- cbind(abs(posinfo$position) <= lim1,                                        # all_center
                  abs(posinfo$position) >= lim2[1] & abs(posinfo$position) <= lim2[2])  # all_sides1

# Fetch bins from positional information matrix and process signal
fetch_data <- function(posinfo, cols_to_get){
    subset_dat <- posinfo$all_normalized[cols_to_get]
    subset_dat[is.na(subset_dat)] <- 0 # Change NA's to 0
    subset_dat[subset_dat < 0] <- 0    # Set negative values to 0
    total_info <- mean(subset_dat)     # Returns the mean of the signal in the region
    return(total_info)
}
all_center <- fetch_data(posinfo, vic_cols[,1])
all_sides  <- fetch_data(posinfo, vic_cols[,2])
pchart     <- all_center + all_sides

all_center_sub <- fetch_data(posinfo_sub, vic_cols[,1])
all_sides_sub  <- fetch_data(posinfo_sub, vic_cols[,2])
pchart_sub <- all_center_sub + all_sides_sub

# Make pChART output
output <- data.frame(center = all_center,
                     side1  = all_sides,
                     pchart = pchart,
                     pchart_sub = pchart_sub,
                     nreads = usedReads,
                     nreads_sub = nrow(subsampled_data),
                     pchart_reads_raw = pchart_reads_raw,
                     pchart_reads_trimmed = pchart_reads_trimmed)

if(!is.null(raw_data$FeatureNumber)){
    output$nmotifs     <- max(raw_data$FeatureNumber)
    output$nmotifs_sub <- length(unique(subsampled_data$FeatureNumber))
}

write.table(output, file = outfile, sep = '\t', row.names = F, quote = F)


## Step 5: Plot information content per region

if(plot){
    enableJIT(0)
    # Subsample data for plotting
    d_to_plot <- subsample(subsampled_data, args$maxpoints)
    n_kept    <- nrow(d_to_plot)
    
    # ggplot-friendly data frames
    posinfo_melt <- entro_data_sub_all %>% 
        gather('type', 'value', 2:4) %>% 
        mutate(type = factor(type, 
                             levels = c("signif_motifs","shuffled","normalized"),
                             labels = c('High signal regions', 'Random', 'Enrichment'),
                             ordered = T)) %>% 
        mutate(category = ifelse(type %in% c('High signal regions', 'Random'), 
                                 'V-plot information\ncontent', 'Chromatin\ninformation (log2)'))
    
    # Add dummy rows for manual axis limits (temporary)
    if(!is.null(args$ylim)){
        dummy <- head(posinfo_melt[posinfo_melt$category %in% "Chromatin\ninformation (log2)",], 2)
        dummy$pos <- 510
        dummy$value <- c(args$ylim, args$ylim * -1)
        posinfo_melt <- rbind(posinfo_melt, dummy)
    }
    
    # Make plots
    p1 <- ggplot(posinfo_melt, aes(x = pos , y = value, color = type)) + 
        geom_line(na.rm = T) + 
        geom_hline(aes(yintercept = 0, linetype = category), posinfo_melt) + 
        labs(x = 'ATAC-seq fragment midpoint position\nrelative to motif center (bp)', 
             y = NULL, color = NULL) +
        facet_wrap(~ category, ncol = 1, strip.position = 'left', scales = 'free_y') +
        scale_x_continuous(breaks = seq(-xrange, xrange, 250), expand = c(0.0125, 0), 
                           limits = c(-args$xlim - 5, args$xlim + 5)) +
        scale_color_manual(breaks = c('High signal regions', 'Random'),
                           values = c("#E41A1C","#377EB8", 'black')) +
        scale_linetype_manual(breaks = c('Enrichment', 'Normalized information'), values = c(2, 0)) +
        guides(linetype = F, color = F) +
        theme_bw(base_size = 14) + 
        theme(
            panel.grid.minor = element_blank(),
            strip.background = element_blank()
        )
    # plot(p1)
    
    p2 <- ggplot(d_to_plot, aes(x     = FragmentMidpointDistanceToFeatureMidpoint, 
                                y     = FragmentSize, 
                                color = FragmentPositionRelativeToFeatureMidpoint)) + 
        geom_point(size = args$size, alpha = args$alpha) + 
        scale_x_continuous(expand = c(0.0125, 0)) +
        scale_y_continuous(expand = c(0.0125, 0)) +
        coord_cartesian(ylim = c(lowerFragmentSizeLimit, args$ylim2),
                        xlim = c(-args$xlim, args$xlim)) +
        scale_colour_manual(values = c("blue", "red", "black"), guide = F) +
        labs(title = title, y = '', x = '', 
             subtitle = paste0(n_kept, ' / ', nrow(raw_data), ' points plotted (', args$maxpoints,
                               ' requested)')) + 
        theme_bw(base_size = 10) + 
        theme(
            strip.background = element_rect(fill="gray90", colour=FALSE),
            panel.border  = element_rect(colour="gray90"),
            panel.spacing = unit(1, "lines"),
            plot.title    = element_text(vjust = 1, hjust = 0.5),
            axis.title.y  = element_text(vjust = 1, margin = margin(l = 8))
        )
    # plot(p2)
    
    layout <- matrix(c(1,1, 1,1, 1,1, 1,1,
                       2,2, 2,2, 2,2), ncol = 2, byrow = T)
    
    png(file = plot_out, width = 6, height = 8, units = 'in', res = 120)
    multiplot(layout = layout, p2, p1)
    dev.off()
    
    if(split == T){
        ggsave(p1, filename = paste(args$output, '.posinfo.pdf', sep = ''), device = cairo_pdf,
               width = 4, height = 3)
        
        p2 <- ggplot(d_to_plot, aes(
          x = FragmentMidpointDistanceToFeatureMidpoint, 
          y = FragmentSize, 
          color = FragmentPositionRelativeToFeatureMidpoint)
        ) + 
          geom_point(size = args$size, alpha = args$alpha) + 
          scale_color_manual(values = c("blue", "red", "black"), guide = F) +
          labs(y = '', x = '') + 
          scale_x_continuous(expand = c(0,0)) +
          scale_y_continuous(expand = c(0,0)) +
          coord_cartesian(ylim = c(lowerFragmentSizeLimit, args$ylim2),
                          xlim = c(-args$xlim, args$xlim)) +
          theme_bw(base_size = 12) + 
          theme(panel.grid.minor = element_blank(),
                axis.ticks = element_blank(),
                panel.border = element_blank(),
                axis.text = element_blank())
        ggsave(p2, filename = paste(args$output, '.vplot.png', sep = ''), width = 5, height = 3.7,
               dpi = 150)
        
    }
}

# Plot randomized data
if(args$randomplot){
    d_shuf_to_plot <- randomize(d_to_plot)
    
    d_shuf_to_plot[,3] <- as.numeric(as.character(d_shuf_to_plot[,3]))
    d_shuf_to_plot[,4] <- as.numeric(as.character(d_shuf_to_plot[,4]))
    
    p_rand <- ggplot(d_shuf_to_plot, aes(x = FragmentMidpointDistanceToFeatureMidpoint, 
                                         y = FragmentSize)) + 
        geom_point(size = args$size, alpha = args$alpha) + 
        scale_x_continuous(expand = c(0.0125, 0)) +
        scale_y_continuous(expand = c(0.0125, 0)) +
        labs(title = title, y = 'Fragment size (bp)', x = 'Fragment distance to motif center (bp)') + 
        theme_bw(base_size = 10) + 
        theme(
            strip.background = element_rect(fill="gray90", colour=FALSE),
            panel.border  = element_rect(colour="gray90"),
            panel.spacing = unit(1, "lines"),
            plot.title    = element_text(vjust=1, hjust = 0.5),
            axis.title.y  = element_text(vjust=1)
        )
    # plot(p_rand)
    ggsave(p_rand, filename = paste0(args$output, '.random.png'), width = 6, height = 5)
}