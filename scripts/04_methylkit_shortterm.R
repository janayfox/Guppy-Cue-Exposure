library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

file.list.shortterm = list("../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC10F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC10M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC11F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC11M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC13F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC13M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC14M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC15F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC15M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC16F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC16M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC1F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC1M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC2F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC2M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC3F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC3M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC4F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC4M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC5F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC5M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC6F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC6M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC7F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC7M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC8F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC8M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC9F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC9M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C10F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C10M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C11F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C11M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C12F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C12M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C13F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C13M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C14F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C14M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C15F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C15M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C1F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C1M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C2F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C2M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C3F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C3M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C4F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C4M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C5F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C5M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C6F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C6M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C7F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C7M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C8F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C8M.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C9F.CpG_report.txt.gz",
                     "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C9M.CpG_report.txt.gz")


myobj_shortterm=methRead(file.list.shortterm,
                   sample.id=list("ST2AC10F","ST2AC10M","ST2AC11F","ST2AC11M","ST2AC13F",
                                  "ST2AC13M","ST2AC14F","ST2AC14M","ST2AC15F","ST2AC15M",
                                  "ST2AC16F","ST2AC16M","ST2AC1F","ST2AC1M","ST2AC2F",
                                  "ST2AC2M","ST2AC3F","ST2AC3M","ST2AC4F","ST2AC4M",
                                  "ST2AC5F","ST2AC5M","ST2AC6F","ST2AC6M","ST2AC7F",
                                  "ST2AC7M","ST2AC8F","ST2AC8M","ST2AC9F","ST2AC9M",
                                  "ST2C10F","ST2C10M","ST2C11F","ST2C11M","ST2C12F",
                                  "ST2C12M","ST2C13F","ST2C13M","ST2C14F","ST2C14M",
                                  "ST2C15F","ST2C15M","ST2C1F","ST2C1M","ST2C2F",
                                  "ST2C2M","ST2C3F","ST2C3M","ST2C4F","ST2C4M",
                                  "ST2C5F","ST2C5M","ST2C6F","ST2C6M","ST2C7F",
                                  "ST2C7M","ST2C8F","ST2C8M","ST2C9F", "ST2C9M"),
                   assembly="guppyWGBS_shortterm",
                   pipeline="bismarkCytosineReport",
                   treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                   context="CpG",
                   mincov = 10,
                   dbtype = "tabix",
                   dbdir = "guppyWGBS_shortterm_DB"
)

saveRDS(myobj_shortterm, file = "./guppyWGBS/myobj_shortterm.rds")

meth=unite(myobj_shortterm, destrand=TRUE)
myDiff=calculateDiffMeth(meth ,mc.cores=10)
diffMeth <- diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
saveRDS(meth, file = "../guppyWGBS/meth_shortterm.rds")
saveRDS(myDiff, file = "../guppyWGBS/myDiff_shortterm.rds")
saveRDS(diffMeth, file = "../guppyWGBS/diffMeth_shortterm.rds")


#meth.min_8=unite(myobj_shortterm,min.per.group=8L, destrand=TRUE)
#myDiff_min_8=calculateDiffMeth(meth.min_8 ,mc.cores=10)
#diffMeth_8 <- diffMethPerChr(myDiff_min_8,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
#saveRDS(meth.min_8, file = "../guppyWGBS/meth.min_8_shortterm.rds")
#saveRDS(myDiff_min_8, file = "../guppyWGBSmyDiff_min_8_shortterm.rds")
#saveRDS(diffMeth_8, file = "../guppyWGBS/diffMeth_8_shortterm.rds")

#meth.min_1=unite(myobj_shortterm,min.per.group=1L, destrand=TRUE)
#myDiff_min_1=calculateDiffMeth(meth.min_1 ,mc.cores=10)
#diffMeth_1 <- diffMethPerChr(myDiff_min_1,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
#saveRDS(meth.min_1, file = "../guppyWGBS/meth.min_1_shortterm.rds")
#saveRDS(myDiff_min_1, file = "../guppyWGBS/myDiff_min_1_shortterm.rds")
#saveRDS(diffMeth_1, file = "../guppyWGBS/diffMeth_1_shortterm.rds")

