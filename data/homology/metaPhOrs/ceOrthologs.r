# Gets Ce orthologs from metaPhOrs database.

source("git/data/name_convert.r")

data.dir = "~/data/ftp/phylomedb.org/metaphors/release-201405/"

# get Ce orthology
ortho.ce = read.table(gzfile(paste0(data.dir, "orthologs/6239.txt.gz")),
  sep="\t", comment.char="", header=TRUE, as.is=TRUE)
colnames(ortho.ce)[1] = "taxid1"

# get mapping from metaPhOrs database to external gene names
# (FIXME: need to convert transcript names to gene names)
p = pipe(paste0("gunzip -c ", data.dir, "ext2meta.txt.gz | fgrep CAEEL"), "r")
ids = read.table(p, sep="\t", comment.char="", header=FALSE, as.is=TRUE)
colnames(ids) = c("extid", "dbid", "version", "protid")

ids$extid = rename.gene.name.vector(ids$extid)



