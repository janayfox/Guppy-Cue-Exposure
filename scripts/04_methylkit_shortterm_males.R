library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

file.list.shortterm = list("../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC10M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC11M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC13M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC14M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC15M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC16M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC1M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC2M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC3M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC4M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC5M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC6M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC7M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC8M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC9M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C10M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C11M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C12M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C13M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C14M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C15M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C1M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C2M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C3M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C4M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C5M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C6M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C7M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C8M.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C9M.CpG_report.txt.gz")


myobj_shortterm=methRead(file.list.shortterm,
                         sample.id=list("ST2AC10M","ST2AC11M",
                                        "ST2AC13M","ST2AC14M","ST2AC15M",
                                        "ST2AC16M","ST2AC1M",
                                        "ST2AC2M","ST2AC3M","ST2AC4M",
                                        "ST2AC5M","ST2AC6M",
                                        "ST2AC7M","ST2AC8M","ST2AC9M",
                                        "ST2C10M","ST2C11M",
                                        "ST2C12M","ST2C13M","ST2C14M",
                                        "ST2C15M","ST2C1M",
                                        "ST2C2M","ST2C3M","ST2C4M",
                                        "ST2C5M","ST2C6M",
                                        "ST2C7M","ST2C8M","ST2C9M"),
                         assembly="guppyWGBS_shortterm_males",
                         pipeline="bismarkCytosineReport",
                         treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                         context="CpG",
                         mincov = 10,
                         dbtype = "tabix",
                         dbdir = "guppyWGBS_shortterm_DB_males"
)

#filter out sites in the 99.9th percentile of coverage (PCR bias)
filt.myobj_shortterm=filterByCoverage(myobj_shortterm,lo.count=NULL,lo.perc=NULL,
                                      hi.count=NULL, hi.perc=99.9)

#normalize by median coverage
norm.filt.myobj_shortterm=normalizeCoverage(filt.myobj_shortterm, method="median")

saveRDS(norm.filt.myobj_shortterm, file = "../guppyWGBS/myobj_norm_filt_shortterm_males.rds")

meth=unite(norm.filt.myobj_shortterm, min.per.group=2L, destrand=TRUE)
myDiff=calculateDiffMeth(meth,mc.cores=2)
diffMeth <- getMethylDiff(myDiff, difference = 15, qvalue = 0.0125)
diffMethChr <- diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=25)

saveRDS(meth, file = "../guppyWGBS/meth_shortterm_males.rds")
saveRDS(myDiff, file = "../guppyWGBS/myDiff_shortterm_males.rds")
saveRDS(diffMeth, file = "../guppyWGBS/diffMeth_shortterm_males.rds")
saveRDS(diffMethChr, file = "../guppyWGBS/diffMethChr_shortterm_males.rds")