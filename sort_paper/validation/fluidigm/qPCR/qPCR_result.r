# Analysis of the qPCR results.

source("git/utils.r")

data.dir = "data/expression/qPCR/20150428 daf-19 hlh-6/"

r1 = read.tsv(paste0(data.dir, "results1.tsv"))
r2 = read.tsv(paste0(data.dir, "results2.tsv"))

# omit some controls
r2 = r2[ r2$"Sample Name" != "NORT" , ]
r2 = r2[ r2$"Sample Name" != "Water" , ]

# slight name tweaking
r1$"Sample Name" = sub("daf1a", "daf19", r1$"Sample Name")
r1$"Sample Name" = sub("hih6", "hlh6", r1$"Sample Name")
r2$"Sample Name" = sub("hih6", "hlh6", r2$"Sample Name")
r1$"Target Name" = sub("hih6", "hlh6", r1$"Target Name")
r2$"Target Name" = sub("Phat5", "phat5", r2$"Target Name")

# replace any rows in r1 that are present in r2
r1$sampleAndTarget = paste(r1$"Sample Name", r1$"Target Name")
r2$sampleAndTarget = paste(r2$"Sample Name", r2$"Target Name")
r2 = r2[ order(match(r2$sampleAndTarget, r1$sampleAndTarget)) , ]

# XXX relies on there being the same number of rows here
i = which(r1$sampleAndTarget %in% r2$sampleAndTarget)
# stopifnot(r1[i,"sampleAndTarget"] == r2[,"sampleAndTarget"])
# r1[i,] = r2


