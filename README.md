#Data and code associated with papers: 
**Developmental Paper:** Developmental behavioural plasticity and DNA methylation patterns in response to predation stress in Trinidadian guppies

**Short-term Paper:** Rapid neural DNA methylation responses to predation stress in Trinidadian guppies

#DOIs: TBD

#Experiment Summary: 

**Developmental Paper:** Juvenile guppies were exposed to either alarm cue or control cue throughout development three times a week for eight weeks. After a 22-week period without exposure to cue guppies went through behavioural assays (open field and shoaling). Immediately following behavioural assays, we dissected out brains for whole genome bisulfite sequencing. 

**Short-term Paper:** Adult guppies were exposed to either alarm cue or control cue one time. Behaviour was measured for 5 minute before and after cue exposure. At set time points following cue exposure (0.5h, 1h, 4h, 24h, 72h), brains were dissected out for whole genome bisulfite sequencing.

#File structure: Raw and clean behavioural and weight data files are provided but can also be downloaded along with metadata from Dryad (Accession # TBD). Sequencing data are stored on SRA (Accession Numbers TBD). Analysis scripts are provided in separate folders for each project (dev_exp and shortTerm_exp) and then each analysis type (01_behavioural_analysis and 02_methylation_analysis). Within each folder, each script should be run in order of numbering. 

All variables and abbreviations defined in metadata. Missing data = NA. 

All required packages and functions are listed at the top of each script. 

#Setup 
1. Download project repository from Github or Dryad (When published: https://doi.org/10.5061/dryad.dbrv15fch)
2. Download data folders:
 - Behavioural data: from Github or Dryad (When published: https://doi.org/10.5061/dryad.dbrv15fch)
 - Sequencing data: from SRA (Accession #TBD)

 #Code Runnnig
 All code for analysis is provided in the github repository or in scripts.zip file on Dryadd. Each folder contains all scripts for a specific type of analysis (01_behavioural_analysis and 02_methylation_analysis) and the files within the folder should be run in order of numbering. Files with a and b lettering represent code paired with its slurm submission file. Note that any code with this pairing was written to be run on Compute Canada HPC system and may need to be adapted to your needs if you do not have access. All required software is commented at top of each script. Uncomment and run install lines if required. Analysis was originally run on R v4.3.2. The nf-core/methylseq pipeline v1.6.1 was used to process sequencing reads. Other details on versions of software used can be found in manuscript.