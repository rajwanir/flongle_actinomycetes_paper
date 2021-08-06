# Scripts used in preparation of flongle actinomycetes paper.




Main text figures:

| Figure no.       |Description    |File          |
| ------------- | ------------- |------------- |
| 1             |Effect of sequencing coverage  |src/plot_assembly_quality.R  |
| 2 | Sequencing run details  |src/sequencing_qc.R  |
| 3 | Assembly stats  |src/plot_assembly_qc_sample.R |
| 4 | CDPS  |src/QTOF/postdereplication_prism.R |
| 5 | lassopeptide  |src/QTOF/ga6002_lassopeptide.R |
| 6 | MiBiG BGCs BLAST  |src/QTOF/postdereplication_mibig.R |

Supplementary figures:
  
  
  
| Figure no.       |Description    |File          |
| ------------- | ------------- |------------- |
| S1             |Genome sizes JGI  |src/avg_genome_sizes.R.R  |
| S2 | Coverage and contig edge  |src/plot_assembly_quality.R  |
| S3 | Read length and assm quality |src/plot_assembly_qc_readlen.R |
| S4 | Read length and contig edge |src/plot_assembly_qc_readlen.R |
| S6 | taxonomy and BGC |src/16s.R |
| S7 | GNPS spectra |src/QTOF/postdereplication_gnps.R |
| S10 | aminocoumarin BGC |src/plot_cluster_aminocoumarin.R |
  
  
  
Parameters for primary analyses:
  
  
  
| Analysis                 |File                |
| ------------------------ | -------------------|
| Base calling             |src/guppy.sh        |
| Demultiplexing           |src/qcat.sh |
| Assembly           |src/canu.sh |
| Assembly polishing          |src/racon.sh and medaka.sh |
| AntiSMASH           |src/antismash.sh |
