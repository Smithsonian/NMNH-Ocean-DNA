# SI-Ocean-DNA
Smithsonian NMNH Ocean DNA - data management and analysis

## Data Management Guide
```mermaid
graph TD;

  GenoHub[**GenoHub**
Demultiplexed and compressed sequence reads in FASTQ format. Files should end in “.fastq.gz” or “.fq.gz”. Is there any consistent naming scheme?]
  Rename("`Rename FASTQ files following Best Practices Guide`")
  Analyses(Run quality/adapter trimming, mitogenome assembly, etc)
  Scratch[**Hydra Scratch**
/scratch/???/USER_ID
40 TB. Not backed up. Might set up automatic file purging to keep space clean.]
  Store[**Hydra Store**
/store/???/USER_ID
40 TB. Not backed up. For non-active projects or large raw data files. Drive system is slower, can't be used for active analysis]
  XDrive[**X Drive**
MAC: \\SI-OCIO-QNAS2\NMNH-DEPARTMENT 
WINDOWS: SMB://SI-OCIO-QNAS2/NMNH-DEPARTMENT
Unlimited storage. Incrementally backed up daily, fully backed up weekly.]

  GenoHub--download raw FASTQ files-->Scratch
  subgraph " "
    Scratch-->Rename
    Rename-->Analyses
    Rename--copy renamed raw FASTQs-->Store
    Analyses--copy trimmed reads and final results-->Store
  end
  Store--Dan runs monthly backup-->XDrive

```
