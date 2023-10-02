#!/bin/bash
# download UCSC scripts
wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/gtfToGenePred
wget http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/genePredToBed

#make them executable
chmod +x ./gtfToGenePred ./genePredToBed 

#convert Gtf to genePred
./gtfToGenePred guppy.gtf.gz guppy.genePred

#convert genPred to bed12
./genePredToBed guppy.genePred guppy.bed12

# sort bed12
sort -k1,1 -k2,2n guppy.bed12 > guppy.sorted.bed