# UTR_library
Project: Designing a comprehensive 3' UTR library of the human genome
Authors:
- Kim Insigne
- Brenda Ji 

The objective of this project is to design a 3' UTR library to test for novel functional regulatory elements that may provide insight into the quantitative effects of genetic variation in the human genome. Since we want to see how each sub-seqeunce of each UTR works, we will use a tiling system to create our library that consists of 150 bp probes with 50 bp overlaps. In addition, we will use Taylor Ward's massively parallel reporter assay to test for the quantitative effects, so the library size will be limited to 100,000 probes. This requires a major reduction in the size of the library, since we start with ~100,000 UTRs alone. To narrow down the library, we decided to take conservation scores into consideration, along with polyadenylation signal sequences, poly(A) sites, and miRNA binding sites, which are all known regulatory elements that affect the gene expression. 

The documented Beaker version of this project is also available here: ***https://pub.beakernotebook.com/publications/9db300e2-6362-11e6-9fcd-8fb7ce4441b8*** (Also available in the UTR_project folder: ur_library_documented.bkr)
This Beaker notebook contains the main code for creating the library. 

This README doc will help with navigation through the UTR_project folder. 

==========================================================================================
Files:
annotated_3UTR.tsv
- I'm not entirely sure how this is different from above haha, but it is used in the main code of this project 
	Format:
	chr,start,end,strand,gene_id,transcript_id
	

annotated_3UTR_v2.csv
- Same as above, but without duplicates or UTRs that are a single base 
	Format:
	chr,start,end,strand,gene_id,transcript_id


annotated_3utrs_distinct_final.csv
- All the distinct 3' UTRs without any overlap with original start and end interval objects as strings parsed	
	Format:
	chr,gene_id,strand,start,end,length
	

annotated_3utrs_distinct.csv
- All the distinct 3' UTRs without any overlap with just the string version of the interval objects
	Format: 
	chr,gene_id,regions,strand

conserved_8mer_3UTR_motifs.txt
- Just the conserved 8mer 3' UTR motifs from the paper/supplemental info:
- http://www.nature.com/nature/journal/v434/n7031/full/nature03441.html 
- ^see supplemental information ~29 for 3' UTR related motifs:
- http://www.nature.com/nature/journal/v434/n7031/extref/nature03441-s1.pdf
- The authors separated the motifs into clusters, which have been divided into txt files in the folder, 3UTR_motif_clusters


cxcl7.fasta
- Stuffer sequence for probes that are less than 150 bps


hg19.apadb_v2_UTR3.txt
- from bed file: http://tools.genxpro.net/apadb/download/ (v2! and then grepped for just 3'UTR PA sites)
	Format:
	Chromosome,PA site start,PA site end,PA site identifier (The hugo name of the corresponding gene and the number of the PA site (the site at the outermost 3' end corresponds to 1)),number of supporting MACE reads,mapping strand,region in the gene (3'UTR, exon, intron…),median in the cluster,mode in the cluster,lost miRNA binding sites if the PA site isn't used (if there are any)


hg38_LO.apadb_v2_3UTR.csv
- Contains the same info as above, but with the PA site start and end lifted over to hg38 from hg19 
	Format:
	chr,PA_start,PA_end,PA_id,MACE_reads,strand,gene_region,cluster_median, cluster_mode, lost_miRNA_binding_sites


longest_UTR_per_gene_PA_sites.csv
- PA sites from APADB merged with the longest_UTR_per_gene. Done by checking all the PA_start locations to see if they are within the range of a given UTR's start and end positions. 
	Format: 
	chr,strand,gene_id,start,end,length,PA_start,PA_end


longest_UTR_per_gene_PAS_range.csv
- PAS locations in the longest_UTR_per_gene based on the 12 putative hexameter variants. At the poster session someone from Hillary's lab mentioned that it might be a good idea to look for a motif at the before the PAS and a GU motif after the PAS to make sure that it's an actual PAS. This file contains just the locations of the PASs that have already been binned into different length ranges for the easiest way of plotting. 
	Format:
	range+50,locs_pas (range.50 in R)


longest_UTR_per_gene.csv
- Based off of annotated_3utrs_distinct_final.csv. The longest (canonical) UTR transcript per gene
	Format:
	chr,strand,gene_id,start,end,length


longest_UTR_tiled_annotated.csv
- "Annotated" version of longest_UTR_tiled_probes.csv with the name parsed and the sequence omitted. 
	Format:
	chr,gene_id,utr_length,probe_start,probe_end


longest_UTR_tiled_probes.csv
- After tiling and stuffing the UTR sequences from longest_UTR_per_gene.csv, this is the cvs form of the dictionary of the probe name and sequence. 
	Format:
	name, sequence


==========================================================================================

Python Scripts:
get_PA_sites.py
get_PAS.py

cons_scores:
get_cons.py (+ get_cons.sh to run on server)

==========================================================================================

R markdowns:
PAS_hists.Rmd
polyA_site_hists.Rmd
utr_freq_length_hists.Rmd

cons_scores:
longest_utrs_cons_graphs.Rmd

==========================================================================================

Beaker files:
utr_library_documented.bkr
- Main code for 3' UTR library project 
UTRNotes.bkr 
- No code, just my notes on background stuff/PAS/polyA/miRNA/etc. May or may not be helpful?

==========================================================================================

Folders:
1. 3UTR_motif_clusters
- 8mer motifs of the clusters as text files
- also contains miRNA unrelated 3' UTR motifs 


2. cons_scores
The ????.phastCons100way.wigFix.gz files are from: 
http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons100way/hg19.100way.phastCons/
I think??—so I'm not sure if they have to be lifted over to hg38?


all_chrs_avg_cons_scores.csv
- Average conservation scores for every 15 bps of the longest UTR transcript per gene in the string form of a list
	Format:
	chr,strand,gene_id,start,end,length,avg_cons_score


all_chrs_total_avg_cons_scores.csv
- The total average conservation score of the longest UTR transcript per gene
	Format:
	chr,strand,gene_id,start,end,length,avg_cons_score


all_probes_cons_scores.csv
- The average conservation score every 15 bps over probes of 150 bps in length. I feel like there could be a bug in the code? Because I think technically there should be 10 averages, but there are only 9?
	Format:
	"chr","gene_id","utr_length","probe_start","probe_end","0","15","30","45","60","75","90","105","120","total_average"


END
==========================================================================================
