library("S4Vectors", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("IRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("GenomicRanges", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

library("methylKit", lib.loc="/home/janayfox/R/x86_64-pc-linux-gnu-library/4.2")

file.list.shortterm = list("../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC10F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC11F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC13F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC14F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC15F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC16F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC1F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC2F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC3F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC4F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC5F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC6F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC7F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC8F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2AC9F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C10F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C11F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C12F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C13F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C14F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C15F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C1F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C2F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C3F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C4F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C5F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C6F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C7F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C8F.CpG_report.txt.gz",
                           "../guppyWGBS/results/bismark_methylation_calls/stranded_CpG_report/ST2C9F.CpG_report.txt.gz")


myobj_shortterm=methRead(file.list.shortterm,
                         sample.id=list("ST2AC10F","ST2AC11F","ST2AC13F",
                                        "ST2AC14F","ST2AC15F",
                                        "ST2AC16F","ST2AC1F","ST2AC2F",
                                        "ST2AC3F","ST2AC4F",
                                        "ST2AC5F","ST2AC6F","ST2AC7F",
                                        "ST2AC8F","ST2AC9F",
                                        "ST2C10F","ST2C11F","ST2C12F",
                                        "ST2C13F","ST2C14F",
                                        "ST2C15F","ST2C1F","ST2C2F",
                                        "ST2C3F","ST2C4F",
                                        "ST2C5F","ST2C6F","ST2C7F",
                                        "ST2C8F","ST2C9F"),
                         assembly="guppyWGBS_shortterm_females",
                         pipeline="bismarkCytosineReport",
                         treatment=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                         context="CpG",
                         mincov = 10,
                         dbtype = "tabix",
                         dbdir = "guppyWGBS_shortterm_females_DB"
)

#filter out sites in the 99.9th percentile of coverage (PCR bias)
filt.myobj_shortterm=filterByCoverage(myobj_shortterm,lo.count=NULL,lo.perc=NULL,
                                hi.count=NULL, hi.perc=99.9)

#normalize by median coverage
norm.filt.myobj_shortterm=normalizeCoverage(filt.myobj_shortterm, method="median")

saveRDS(norm.filt.myobj_shortterm, file = "../guppyWGBS/myobj_norm_filt_shortterm_females.rds")

meth=unite(norm.filt.myobj_shortterm, min.per.group=2L, destrand=TRUE)
myDiff=calculateDiffMeth(meth ,mc.cores=2)
diffMeth <- getMethylDiff(myDiff, difference = 15, qvalue = 0.0125)
diffMethChr <- diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.0125, meth.cutoff=25)

saveRDS(meth, file = "../guppyWGBS/meth_shortterm_females.rds")
saveRDS(myDiff, file = "../guppyWGBS/myDiff_shortterm_females.rds")
saveRDS(diffMeth, file = "../guppyWGBS/diffMeth_shortterm_females.rds")
saveRDS(diffMethChr, file = "../guppyWGBS/diffMethChr_shortterm_females.rds")
