## preprocess gnomad SVs
# split gnomAD SV bed file into separate bed files by SV type

gnomad_bed = $1
gunzip $gnomad_bed

awk -F"\t" '$5 == "DEL" { print $1"\t"$2"\t"$3"\t"$4 }' gnomad_v2_sv.sites.bed > gnomad_deletions.bed

awk -F"\t" '$5 == "INV" { print $1"\t"$2"\t"$3"\t"$4 }' gnomad_v2_sv.sites.bed > gnomad_inversions.bed

awk -F"\t" '$5 == "DUP" { print $1"\t"$2"\t"$3"\t"$4 }' gnomad_v2_sv.sites.bed > gnomad_duplications.bed

awk -F"\t" '$9 ~ /melt/  { print $1"\t"$2"\t"$3"\t"$4 }' gnomad_v2_sv.sites.bed > gnomad_meis.bed

sed -i 's/^/chr/g' gnomad_deletions.bed
sed -i 's/^/chr/g' gnomad_inversions.bed
sed -i 's/^/chr/g' gnomad_duplications.bed
sed -i 's/^/chr/g' gnomad_meis.bed








