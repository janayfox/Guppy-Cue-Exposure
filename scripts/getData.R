library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

meth.unite <- readMethylDB("./guppyWGBS_dev_DB/methylBase_meth_unit.txt.bgz")

meth.data <- getData(meth.unite)

saveRDS(meth.data, file = "../guppyWGBS/meth_data.rds")
