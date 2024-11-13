# SI-Ocean-DNA
Smithsonian NMNH Ocean DNA - data management and analysis

```mermaid
graph TD;

  GenoHub[**GenoHub**
Raw, demultiplexed sequence reads in FASTQ format. Files shoule end in “.fastq.gz” or “.fq.gz”. Is there any consistent naming scheme?]
  Rename("`Rename FASTQ files following Best Practices Guide`")
  Scratch[**Scratch**]
  Store[**Store**]
  XDrive[**X Drive**]

  GenoHub--download FASTQ files-->Rename
  subgraph Hydra
    Rename-->Scratch
    Rename-->Store
    Scratch-->Store
  end
  Store-->XDrive

```
