##############################
## Convert 1kg SVs to bed file with useful annotations
## may need to setwd()
##############################
library(vcfR)
##############################
kg = read.vcfR('ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf') # change path as needed
meta = kg@fix
meta = as.data.frame(meta, stringsAsFactors = F)
rm(kg)
meta['#chr'] = paste("chr", meta$CHROM, sep="")
meta$start = as.integer(meta$POS)
meta$deletion = ifelse(grepl("<CN[01]>", meta$ALT), T, F)
meta$duplication = ifelse(grepl("<CN[3456789]>", meta$ALT), T, F)
meta$mei = grepl("<INS:ME:.*$", meta$ALT)
meta$inversion = grepl("<INV>", meta$ALT)
meta$end = gsub("(^.*\\;END=)(\\d*)(\\;.*$)", "\\2", meta$INFO)

meta$svtype = gsub("(^.*SVTYPE=)(.*)(\\;.*)", "\\2", meta$INFO)
while (any(grepl("\\;", meta$svtype))){
  meta$svtype = gsub("(^.*)(\\;.*$)", "\\1", meta$svtype)
}

## stack
dels = meta[meta$deletion,]
dels$SV_type = "DELETION"
dups = meta[meta$duplication,]
dups$SV_type = "DUPLICATION"
invs = meta[meta$inversion,]
invs$SV_type = "INVERSION"
meis = meta[meta$mei,]
meis$SV_type = "MEI"
meis$end = meis$start + 1
stack = rbind.data.frame(dels, dups, invs, meis, stringsAsFactors = F)
stack.sub = stack[,c("#chr", "start", "end", "ID" ,"SV_type")]


for(sv in unique(stack.sub$SV_type)){
  write.table(stack.sub[stack.sub$SV_type==sv,], paste("1kg_", sv, ".bed", sep=""), quote=F, sep="\t", row.names=F)
}

