# SI-Ocean-DNA
Smithsonian NMNH Ocean DNA - data management and analysis

## Sequence Data File Names -- Best Practices Guide

Required fields:
- Voucher/Catalog ID
- Taxonomic ID (Family-Genus-Species?)
- anything else?

NO SPACES in file name

Use delimiter (underscore? period? dash?) to seperate fields in the file name. Delimiter must NOT be used within the fields.

## Data Management Guide
```mermaid
graph TD;

  GenoHub[**GenoHub**
Demultiplexed and compressed sequence reads in FASTQ format. Files should end in “.fastq.gz” or “.fq.gz”. Is there any consistent naming scheme?]
  Rename("`Rename FASTQ files following Best Practices Guide`")
  Analyses(Run quality/adapter trimming, mitogenome assembly, etc)
  Scratch[(**Hydra Scratch**
/scratch/???/USER_ID
40 TB. Not backed up. Might set up automatic file purging to keep space clean.)]
  Store[(**Hydra Store**
/store/???/USER_ID
40 TB. Not backed up. For non-active projects or large raw data files. Drive system is slower, can't be used for active analysis)]
  XDrive[(**X Drive**

MAC -- &bsol;&bsol;SI-OCIO-QNAS2&bsol;NMNH-DEPARTMENT
WINDOWS -- SMB://SI-OCIO-QNAS2/NMNH-DEPARTMENT
Unlimited storage? Incrementally backed up daily, fully backed up weekly.)]

  Move1[/download raw FASTQ files/]
  Move2[/copy renamed raw FASTQs/]
  Move3[/copy clean reads and final results/]
  Move4[/Dan runs monthly backup/]

  GenoHub-->Move1
  Move1-->Scratch
  subgraph " "
    Scratch-->Rename
    Rename-->Analyses
    Rename-->Move2
    Move2-->Store
    Analyses-->Move3
    Move3-->Store
  end
  Store-->Move4
  Move4-->XDrive

  classDef process stroke:black,color:white,fill:#159BD7,stroke-dasharray: 5 5
  classDef storage stroke:black,color:white,fill:#159BD7
  classDef ccr stroke:black,color:white,fill:#159BD7

  class Rename,Analyses,Move1,Move2,Move3,Move4 process
  class GenoHub,Scratch,Store,XDrive storage
```
