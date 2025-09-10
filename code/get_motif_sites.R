library(TOP)
library(data.table)
library(optparse)
options(scipen=999)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--tf",type="character")
parser <- add_option(parser,"--motif",type="character")
parser <- add_option(parser,"--motif_dir",type="character")
parser <- add_option(parser,"--threshP",type="character",default="1e-5")
parser <- add_option(parser,"--threshPWM",type="character",default="10")
parser <- add_option(parser,"--flank",type="integer",default=100)
parser <- add_option(parser,"--outdir",type="character")
out    <- parse_args(parser)
tf_name     <- out$tf
motif_name  <- out$motif
motif_dir  <- out$motif_dir
thresh_pValue  <- out$threshP
thresh_pwmscore <- out$threshPWM
flank          <- out$flank
output_dir     <- out$outdir
# print(out)
rm(parser,out)

# Example:
# tf_name = 'CTCF'
# motif_name = 'MA0139.1'
# motif_dir = '/project2/xinhe/kevinluo/footprint_clustering/data/motifs/JASPAR2024_CORE_non-redundant_pfms_meme'
# thresh_pValue = 1e-5
# flank = 100
# output_dir = '/project2/xinhe/kevinluo/footprint_clustering/processed_data/hg38'

sequence_file = '/project2/xinhe/kevinluo/footprint_clustering/data/ref_genome/hg38.fa'
blacklist_file = '/project2/xinhe/kevinluo/footprint_clustering/data/ENCODE_blacklist/blacklist.hg38.bed.gz'
outname = paste(tf_name, motif_name, thresh_pValue, sep = '_')

chromsize_file = '/project2/xinhe/kevinluo/footprint_clustering/data/ref_genome/hg38.chrom.sizes'
if (!file.exists(chromsize_file)) {
  index_fa(sequence_file, chromsize_file=chromsize_file)
}

# Find motif matches
cat("Find motif matches ...\n")
fimo_motif_matches(motif_file=file.path(motif_dir, paste0(motif_name, ".meme")),
                   sequence_file=sequence_file,
                   thresh_pValue=as.numeric(thresh_pValue),
                   outname=file.path(output_dir, paste0(motif_name, '_', thresh_pValue, '.fimo.txt')),
                   fimo_path='fimo')

# Get candidate sites
cat("Get candidate sites ...\n")
sites <- process_candidate_sites(fimo_file=file.path(output_dir, paste0(motif_name, '_', thresh_pValue, '.fimo.txt')),
                                 flank=flank,
                                 thresh_pValue=as.numeric(thresh_pValue),
                                 thresh_pwmscore=as.numeric(thresh_pwmscore),
                                 blacklist_file=blacklist_file,
                                 chr_order = paste0('chr', c(1:22)))
cat(nrow(sites), "candidate sites.\n")
saveRDS(sites, file = file.path(output_dir, paste0(outname, '.candidate.sites.rds')))

data.table::fwrite(sites[,1:6],
                   file.path(output_dir, paste0(outname, '.candidate.sites.bed')), sep = '\t',
                   col.names = FALSE)
