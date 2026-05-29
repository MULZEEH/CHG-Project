# install.packages("DNAcopy")

library(DNAcopy)

# getting args from command line
# @params: copynumber_file: path to the input file with copy number data
# @params: output_folder: path to the folder where the output files will be saved

# USAGE
#
# Usage: Rscript CBS.R <copynumber_file> <output_folder>
#

args <- commandArgs(trailingOnly = TRUE)

cn=read.table(args[1],header=TRUE,fill=TRUE,stringsAsFactor=FALSE)
output_folder <- args[2]

# folder = "~/Documents/HumanGenomics/05_SomaticCopyNumberCalling/Data/"
# cn <- read.table(file.path(folder,"SCNA.copynumber.called"),header=T)

# pdf(file.path(output_folder,"SegPlot.pdf"))

# plot(cn$raw_ratio,pch=".",ylim=c(-2.5,2.5))
# plot(cn$adjusted_log_ratio,pch=".",ylim=c(-2.5,2.5))
CNA.object <-CNA(genomdat = cn$adjusted_log_ratio, 
                 chrom = cn$chrom,
                 maploc = cn$chr_start, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, min.width=2,
                undo.splits="sdundo", #undoes splits that are not at least this many SDs apart.
                undo.SD=3,verbose=1)

# plot(segs,plot.type="w")

# dev.off()

segs2 = segs$output
write.table(segs2, file=file.path(output_folder,"SCNA.copynumber.called.seg"), row.names=F, col.names=T, quote=F, sep="\t")