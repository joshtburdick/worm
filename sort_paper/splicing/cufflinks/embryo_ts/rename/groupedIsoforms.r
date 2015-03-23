# Groups isoforms to be consistent with the nomenclature used
# in the paper.

source("git/utils.r")

r = read.tsv("git/sort_paper/splicing/cufflinks/embryo_ts/cufflinksSummary.tsv")
colnames(r)[1] = "Time"
colnames(r)[2] = "Isoform_name"

# mapping from UCSC RefSeq ID to isoform name used in the paper
f = NULL
f[ c("c") ] = "LIN-3A"
f[ c("a", "b", "e") ] = "LIN-3B"
f[ c("d") ] = "LIN-3C"
f[ c("f.1", "f.2") ] = "LIN-3D"
names(f) = paste0("F36H1.4", names(f))

r$"Gene" = f[ r$Isoform_name ]

# sum these, and convert to a table
r1 = t(tapply(r$FPKM, list(Gene=r$Gene, Time=r$Time), sum))




write.tsv(r1, file="git/sort_paper/splicing/cufflinks/embryo_ts/rename/groupedIsoforms.tsv")


