# lzhang
# Shogo says A22i file should have a 4th column of the renumbered ID
# otherwise see this error when using OPTION readAInverse [file]
# Reading A22-inverse from file: "../pregs_gbs_a22/Gi" # this is a binary file
# Wrong file or format
# header in the file:GINVERSE  ; expected header:A22INVERSE

library(dplyr)

###### read
# fake a22i
a22i <- read.table("blupf90/vce/pregs_gbs_a22/Gi_ascii", header=F, 
                   col.names=c("i", "j", "gii"))

xrefid <- read.table("output/gbs_a22_clean_XrefID", header = F,
                     col.names = c("id", "id_orig"))
xrefid_hap <- read.table("output/haplotype_clean_XrefID", header = F,
                     col.names = c("id", "id_orig"))

all(xrefid==xrefid_hap) # sanity check

xrefid$index <- 1:nrow(xrefid)
