# Gene transcription responses of *Daphnia galeata* to the prolonged presence of natural predators
## Introduction
*Daphnia* is a genus of branchiopod crustaceans belonging to the sub-order Cladocera (order Diplotraca) and comprises the so-called water fleas. Daphniids exhibit a folded (‘bivalve’) carapace that encloses the entire trunk (but not the cephalon as it does in other diplostracans) and serves mainly as a brood chamber (Brusca et al., 2016). Water fleas are distributed worldwide in freshwater basins and ponds where they play a key role as bioindicators of environmental health and quality (Miyakawa et al., 2010 and references therein). In addition, *Daphnia pulex* is the first crustacean sensu stricto to have had its genome sequenced and, together with other congeneric species, now counts on a personal genome web-database, the [wFleaBase](http://wfleabase.org/) (Colbourne et al., 2005).
Besides genomics and ecotoxicology, water fleas have become also good models to study developmental processes and their interactions with environment, in the modern context of the Eco-Evo-Devo (ecological evolutionary developmental biology). As a matter of fact, *Daphnia* spp. show the capacity to form defensive morphological structures (such as the helmets of *Daphnia cucullata* and *Daphnia longispina* and the neck teeth of *Daphnia pulex*; Tams et al., 2020) in response to the detection of non-lethal chemical cues (named kairomoes) from predators. These morphological defenses in *Daphnia*, also known as predator-induced polyphenisms (Miyakawa et al., 2010), act by lowering the capture success of predators and exhibit a transgen-erational effect too (Agrawal et al., 1999).

In this work, a differential expression analysis of one clone of *Daphnia galeata* were performed in order to address the genetic response and basis of this water flea species to the exposure of fish kairomones. *Daphnia galeata* shows in fact life history trait (such as the age at first reproduction, the somatic growth rate and the body length) variations in populations grown up in the presence of fish kairomones, but not severe morphological changes (Tams et al., 2018, 2020).

## Materials and Methods
### Data collection
Total RNA were extracted from pools of 20 egg-bearing adults of *Daphnia galeata* and sequencing was thereafter performed for 12 samples, assuring three replicates per genotype (M6 and M9) per condition of growth (control environment and kairomone-rich environment; Table 1).
| RUN ACCESSION NUMBER |	GENOTYPE |	CONDITION OF GROWTH |	ID |
|:--------------------:|:---------:|:--------------------:|:--:|
|ERR2929116     	     |M9         |Control               |M9_C|
|ERR2929117*           |M6         |Control	              |M6_C|
|ERR2929118            |M9         |Control	              |M9_C|
|ERR2929119*	         |M6         |Control	              |M6_C|
|ERR2929120*	         |M6         |Control	              |M6_C|
|ERR2929121	           |M9	       |Control	              |M9_C|
|ERR2929122*	         |M6	       |Kairomone	            |M6_F|
|ERR2929123	           |M9	       |Kairomone	            |M9_F|
|ERR2929124*	         |M6	       |Kairomone	            |M6_F|
|ERR2929125	           |M9	       |Kairomone	            |M9_F|
|ERR2929126*	         |M6	       |Kairomone	            |M6_F|
|ERR2929127	           |M9	       |Kairomone	            |M9_F|

Raw RNA-sequencing reads of *Daphnia galeata* (BioProject: [PRJEB29887](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB29887/)) were downloaded from NCBI SRA. **Prefetch v2.8.0** and **FASTQ-dump v2.8.0** from a Bioconda (Grüning et al., 2018) installation of the NCBI SRA-toolkit were used to retrieve SRA and FASTQ files, respectively. FASTQ files were given a header syntax suitable for the subsequent Trinity assembly (<code>--defline-seq '@$sn[_$rn]/$ri'</code>) and were split to keep paired reads separated (<code>--split-files</code>). **vdb-validate v2.8.0** was used to analyze the downloaded SRA data for corruption and other problems, while **FastQC v0.11.7** was used to check read quality. Only 6 runs out of 12 were eventually downloaded: 3 (ERR2929117, ERR2929119 and ERR2929120) were obtained from a *Daphnia galeata* clone bred in absence of fish kairomones (control); conversely, 3 (ERR2929122, ERR2929124, ERR2929126) were obtained from a *Daphnia galeata* clone bred in presence of ide (*Leucuscus idus*) kairomones. The whole process of downloading and quality checks was auto-mated as following:

~~~
#download sras and fastqs, then check read quality
  for i in ERR29291{17,19,20,22,24,26}; do
  #create a directory for every run
	mkdir "$i"
	cd "$i"
  
	#download sra files and validate them
	prefetch "$i"
	vdb-validate "$i".sra
  
	#get fastqs from sra and check their quality
	fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files "$i".sra
	fastqc *fastq -o . -f fastq
  
	#move to parental directory
	cd ..
done

#move read quality checks into a dedicate folder
mkdir reads_quality_check
mv */*html reads_quality_check/
~~~

### Read trimming and filtering
A Bioconda installation of **Trimmomatic v0.39** (Bolger et al., 2014) was used to filter low quality reads and remove Illumina adapters. The first 3 leading and the 5 trailing low-quality bases were removed and the quality of other bases were checked using a sliding window of length 4 and a required average quality of 15. The minimum length of reads to keep was set to 75, considering a total read length of 101 (LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75). The process was automated through all the downloaded FASTQs as following:

~~~
#trimmomatic through every fastq file. NB: the current directory contains a folder for every run named after its accession number
for i in ERR29291*; do

	#create a directory to store trimmed fastqs for every run
	mkdir "$i"_trimmed
	
	trimmomatic PE -threads 5 -phred33 "$i"/"$i"_1.fastq "$i"/"$i"_2.fastq "$i"_trimmed/"$i"_1_trim.fastq "$i"_trimmed/unpaired_1 "$i"_trimmed/"$i"_2_trim.fastq "$i"_trimmed/unpaired_2 ILLUMINACLIP:/usr/local/anaconda3/share/trimmomatic-0.39-2/adapters/ TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 2> "$i"_trimmed/stats_trimmomatic

	#remove upaired reads
	rm "$i"_trimmed/unpaired*

done

#remove all the directories with raw reads
rm ERR29291{17,19,20,22,24,26}
~~~

### Reference transcriptome de-novo assembly
In order to obtain a reference transcriptome for subsequent analyses, all reads from the various runs were merged into a single FASTQ file; read pairs were kept separated. The obtained FASTQ files had a size of 46GB each, which correspond to more than 150M reads each. Thus, FASTQ files were subsampled in order to obtain two files of 10GB each (corresponding to more than 33M reads), which make the assembler better perform. A de-novo assembly was then produced using a Bioconda installation of **Trinity v2.1.1** (Grabherr et al., 2011), with a maximum accessible memory set to 30GB (~1GB RAM per 1M Illumina reads, as Trinity user manual suggests; accessed 06/06/2021). The following command lines were launch in this phase:

~~~
#merge fastqs in one file; NB: the current directory contains a folder for every trimmed run
cat */*1_trim.fastq > all_1_trim.fastq
cat */*2_trim.fastq > all_2_trim.fastq

#count the number of reads for every merged fastq
grep -c "@" all_*

#total number of reads each: 154592585 (46G)
#reads to keep each: 154592585/4.6 = 33607084

#read subsampling
bash random_subsampling_PE.sh all_1_trim.fastq all_2_trim.fastq 33607084

#trinity de-novo assembly
Trinity --seqType fq --left all_1_trim.fastq --right all_2_trim.fastq --CPU 10 \   --max_memory 30G
~~~

The <code>random_subsampling_PE.sh</code> bash script can be accessed [here](scripts/random_subsampling_PE.sh).

### Assembly completeness and quality assessment
Trinity output FASTA was onelinerized and then manipulated in order to obtain headers in the form of <code>>Dgal_trinityID</code>. The following command lines were run:

~~~
#rename trinity output headers to obtain ">Dgal_trinityID"; NB: dgal_ref.fna is the Trinity FASTA output
sed -E 's/ len.+$//; s/_/./g; s/^>/>Dgal_/' dgal_ref.fna > dgal_ref_rn.fna

#onelinerize FASTA
cat dgal_ref.rn.fna | awk '/^>/ {printf("\n%s\n",$0);next;} {printf("%s",$0);} END {printf("\n");}' | tail -n +2 > dgal_ref.rn.oneline.fna
~~~

Trinity de-novo assembly quality ad completeness was checked using **BUSCO v5** (Seppey et al., 2019) on the web application **gVolante2** (Nishimura et al., 2019). ‘arthropoda_odb10’ was used as the reference ortholog set and the sequence type was set to ‘coding/transcribed (nucleotide)’.
After this first assessment, the redundancy of the reference transcriptome assembly was reduced using a Bioconda installation of **CD-HIT v4.8.1** (Fu et al., 2012), with a sequence identity threshold of 0.9 and a tolerance for redundancy of 1.

~~~
#reduce redundancy
cd-hit-est -i dgal_ref_rn_oneline.fna -o dgal_ref_rn_oneline_rd.fna -T 12 -t 1 -c 0.9
~~~

### Read mapping
Aligning sequenced reads to a reference genome/transcriptome is the first step in many comparative genomics pipelines, including pipelines for variant calling, isoform quantitation and differential gene expression (Langmead & Salzberg, 2012). Here, downloaded reads were mapped to the Trinity reference transcriptome using a Bioconda installation of **Bowtie2 v2.4.2** (Langmead & Salzberg, 2012). Reference transcriptome were indexed using the Bowtie built-in tool (<code>bowtie2-build</code>) and reads were then mapped disabling the discordant alignments.

~~~
#index reference fasta
bowtie2-build dgal_ref_rn_oneline_rd.fna dgal_ref_rn_oneline_rd.fna

for i in ERR29291{17,19,20,22,24,26}; do
	#map reads to the reference	
	bowtie2 -x dgal_ref_rn_oneline_rd.fna -1 "$i"_1_trim.fastq - "$i"_2_trim.fastq  -S "$i"_mapped_bt2.sam --no-discordant -p 12
done
~~~

The resulting Sequence Alignment/Map (SAM) files were soon after converted into the less disk-space demanding Binary Alignment/Map (BAM) format using a Bioconda installation of **Mtools v1.3**Li et al., 2009). BAM files were then sorted and filtered to keep only (i) proper paired mapping reads (<code>-f 0x2</code>) and (ii) reads with a minimum mapping quality of 30 (<code>-q 30</code>); conversely, the filtering excluded (iii) secondary alignments (<code>-F 0x200</code>). Raw count statistics were retrieved for both raw mapping reads and filtered mapping reads. The process was automated as following:

~~~
#convert sam to bam, then sort, index and filter bam files. NB: the current directory tree is as follow:
#mapping/
#├── 01_dgal_ref_rn_oneline_rd.fna
#├── 02_raw_mapping
#│	└── {raw bam files}
#├── 03_sorted_mapping
#│	├── 01_rawcount_stats
#│	│    └── {raw count statistics}
#│	└── {sorted and indexed bam files}
#├── 04_sorted_filtered_mapping
#│	├── 01_rawcount_stats
#│	│    └── {raw count statistics}
#│	└── {sorted and indexed and filtered bam files}
#└── {sam files}

for i in ERR29291{17,19,20,22,24,26}; do

	#convert sam to bam
	samtools view -b "$i"_mapped_bt2.sam > 02_raw_mapping/"$i"_mapped_bt2.bam

	#sort bam
	samtools sort -@ 10 02_raw_mapping/"$i"_mapped_bt2.bam > 03_sorted_mapping/"$i"_mapped_bt2_sorted.bam
	
	#index sorted bam
	samtools index 03_sorted_mapping/"$i"_mapped_bt2_sorted.bam
	
	#get raw count statistics of sorted bam
	samtools idxstats 03_sorted_mapping/"$i"_mapped_bt2_sorted.bam > 03_sorted_mapping/01_rawcount_stats/"$i"_rawmapping_stats.tsv
	
	#filter sorted bam
	samtools view -h -@ 10 -f 0x2 -F 256 -q 30 -b 03_sorted_mapping/"$i"_mapped_bt2_sorted.bam > 04_sorted_filtered_mapping/"$i"_mapped_bt2_sorted_filtered.bam
	
	#index sorted and filtered bam
	samtools index 04_sorted_filtered_mapping/"$i"_mapped_bt2_sorted_filtered.bam
	
	#get raw count statistics of sorted and filtered bam
	samtools idxstats 04_sorted_filtered_mapping/"$i"_mapped_bt2_sorted_filtered.bam > 04_sorted_filtered_mapping/01_rawcount_stats/"$i"_rawmapping_stats.tsv
done

#remove sam files
rm *sam
~~~

Raw count statistics for both unfiltered and filtered BAMs were merged into three different tab-separeted-value files, one for mapping unfiltered reads, one for unmapping unfiltered reads and one for filtered reads. The Unix <code>join</code> command was use to this purpose.

~~~
#the code hereafter refers to filtered read-count files; it was then ap-plied to the other raw count statistics as well
#
#remove the last line of raw count file (which reports the total amount of unmapped reads) and then sort it
for i in ERR29291*; do sed '$d' "$i" | sort > "$i"_sorted; done

#join files two by two
join -j 1 -o 1.1,1.3,2.3 ERR2929117_rawmapping_stats_filtered.tsv_sorted ERR2929119_rawmapping_stats_filtered.tsv_sorted > 1719_joined.tmp
join -j 1 -o 1.1,1.2,1.3,2.3 1719_joined.tmp ERR2929120_rawmapping_stats_filtered.tsv_sorted > 171920_joined.tmp
join -j 1 -o 1.1,1.2,1.3,1.4,2.3 171920_joined.tmp ERR2929122_rawmapping_stats_filtered.tsv_sorted > 17192022_joined.tmp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.3 17192022_joined.tmp ERR2929124_rawmapping_stats_filtered.tsv_sorted > 1719202224_joined.tmp
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.3 1719202224_joined.tmp ERR2929126_rawmapping_stats_filtered.tsv_sorted > all_rawcounts_filtered.txt

rm *tmp
~~~

### Transcriptome GO annotation
The obtained reference transcriptome was scanned in order to get open reading frames (ORFs), which were subsequently used to get a gene ontology (GO) annotation. ORFs encoding at least 100 amino acids were first extracted using a Bioconda installation of **[TransDecod-er.LongOrfs v5.5.0](https://github.com/TransDecoder/TransDecoder/wiki)** and then annotated using a Bioconda installation of **Diamond v2.0.8.146** (Buchfink et al., 2015) against the UniProt database (accessed on 07/06/2021). The e-value threshold for the Diamond search was set to 10–5 (<code>--evalue 1e-05</code>) and just the best hit for every query sequence was kept (<code>--max-target-seqs 1</code>). The likely coding regions of ORFs were then predicted using <code>TransDecoder.Predict</code> and integrating the Diamond homology inference.

~~~
#get ORFs from reference transcriptome
TransDecoder.LongOrfs -t 01_dgal_ref_rn_oneline_rd.fna

#annotate ORFs with diamond
diamond blastp --query 01_dgal_ref_rn_oneline_rd.fna.transdecoder_dir/longest_orfs.pep --db /var/local/uniprot/uniprot_sprot.fasta-norep_per_diamond.dmnd --evalue 1e-05 --max-target-seqs 1 --threads 10 --outfmt 6 --out 02_dgal_diamond_uniprot_annotation.tsv

#predict coding regions
TransDecoder.Predict -t 01_dgal_ref_rn_oneline_rd.fna --retain_blastp_hits 02_dgal_diamond_uniprot_annotation.tsv

#onelinerize predicted protein fasta file
cat 01_dgal_ref_rn_oneline_rd.fna.transdecoder.pep | awk '/^>/ {printf("\n%s\n",$0);next;} {printf("%s",$0);}  END {printf("\n");}' | tail -n +2 > 01_dgal_ref_rn_oneline_rd.fna.transdecoder_oneline.pep

#customize fasta entries to keep only geneID
sed -Ei '/^>/ s/ .+$//' 01_dgal_ref_rn_oneline_rd.fna.transdecoder_oneline.pep
~~~

The obtained predicted proteome was then used to perform a GO annotation on the web application [PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/) with the Argot scoring function (Törönen et al., 2018).

### Differential expression analysis
The differential expression analysis was performed with the R package **NOIseq v2.36.0** (Tarazona et al., 2015). The dataset was first filtered using the built-in function <code>filtered.data</code> to keep only reads with a count-per-million (CPM) greater than 10 and with a cutoff for the variation coefficient of 100. False discovery rate (FDR) was used as the method for the multiple test correction. Data quality was checked both before and after filtering using the built-in function
<code>explo.plot</code> to plot saturation values and raw count distributions among samples.
Data were then normalized using the between-sample Trimmed Mean of M-values (TMM) method implemented in the built-in function <code>tmm</code>. No length correction was applied.
The differential expression analysis itself was afterward performed using the built-in function <code>noiseqbio</code>; null values were replaced with 0.1. Differentially expressed features were then re-covered using the built-in function <code>degenes</code> with a q-value of 0.95; all differentially expressed features were returned. Differential expression results were then visualized using the built-in function <code>DE.plot</code>.

The R code used to perform the differential expression analysis can be found [here](scripts/diffExpr.analysis.R).

Genes found to be differentially expressed were then analysed to check whether they show an enrichment in some functional annotation. The R package **topGO v2.44** (Alexa & Rahnenfuhrer, 2021) was therefore run using the PANNZER2 functional annotation of the refer-ence proteome and the NOIseq results. FDRs for every differently expressed gene were used as gene scores to perform the Kolmogorov-Smirnov test statistics. Furthermore, each GO category was tested independently. The final set of best enriched GO terms was selected for every GO major class on the basis of FDR < 0.05 by visually checking the results. Eventually, the results of the GO enrichment analyses were visualized using the web application **[REVIGO](http://revigo.irb.hr/)** (Supek et al., 2011) against the whole UniProt database with an allowed similarity of 0.7.

The R code used to perform the Go enrichment analysis can be found [here](scripts/GOenrichment.R).

## Results and discussion
### Read and assembly general statistics
Raw read quality was generally high across all runs and the trimming process keep about the 90% of the initial raw reads (figure **[here](trimming_statistics/stats_trimmomatic_summ.txt)**). The raw reference transcriptome shows a BUSCO score of 59.9% (in term of single copy orthologs) but after the redundancy removal, this value grows up to 78.9%. The overall BUSCO completeness score (in term of both single-copy and duplicated orthologs) is close to 98% (figure **[here](trimming_statistics/stats_trimmomatic_summ.txt)**). The transcriptome showed a GC content of ~42% (which is consistent with the GC content of other daphniids \[e.g., Lee et al., 2019; Ye et al., 2017]) and an N50 value of ~2k nt. The filtered tran-scriptome used as reference for the subsequent analyses consists of 46,613 nucleotide sequences. From this, 25,262 coding sequences were retrieved by the TransDecoder analysis.

The number of reads that successfully mapped onto the reference transcriptome were 343,219,541, while those that don’t were 5,237,449 (figure **[here](mapping_statistics/mapping.png)**. After filtering, the total amount of mapping reads dropped to 269,526,372.

The number of transcripts showing a CPM value greater than 10 and a variation coefficient between the two sample growth conditions (control environment and kairomone-rich environment) greater than 100 were 7,690 out of 46,613. Sequencing depth proved to be very high even without any filtering operation, since it was greater than 78% for every samples. As a matter of fact, after filtering, sequencing depth reached 100% for every sample. Interestingly, three samples (ERR2929117, ERR2929122, ERR2929124) show an higher number of detected features than the other three (ERR2929119, ERR2929120, ERR2929126), however with no links with the condition of growth.

### Differentially expressed genes

The NOISeq analysis retrieved 1,034 differentially expressed features between the control and the kairomone-exposed water fleas, most of which resulted from an up-regulation of that feature in the kairomone-exposed samples in comparison to the control.

The topGO analysis eventually identified a significant enrichment in 63 GO categories, 35 of which belongs to the biological process major class, 10 to the cellular component and 18 to the molecular function. Some differentially expressed features were related to the response to toxic substances and to detoxification processes. This finding may suggest that kairomones are molecules that trigger Daphnia galeata metabolic pathways related to the removal of xenobiotics, as they are detected by water fleas as foreign compounds. In addition, the analyses found an enrichment in features related to (i) the response to lipids, (ii) to proteolysis and (iii) to carbohydrate binding processes, suggesting that *Daphnia galeata* may be able to sense a variety of fish-released molecules. A significant enrichment was also found in processes related to blood vessel morphogenesis and vasculature development. Despite these categories are strongly related to vertebrate animals and their circulatory system, the genes involved may be likely acting on some other *Daphnia galeata* morphogenetic processes which control somatic growth rate and the general developmental pattern. Not by chance, *Daphnia galeata* shows a decrease of the (i) age at first reproduction, (ii) somatic growth rate and (iii) body length when in presence of fish kairomones (Tams et al., 2018). Accordingly, an enrichment in the (i) structural molecule activity, (ii) structural constituents of cuticle and components of membranes has been found as well.

Kairomones and the effect they produce on preys is a sparkling example of the role the environment have on certain biological processes, such as population structure, feeding and repro-ductive modes (Yousuf et al., 2020). In the case here investigated, fish kairomones seems to act mostly on *Daphnia galeata* developmental processes as well as on many detoxification pathways. The actual molecular and genetic responses needs however further investigations in order to characterize in a deeper way the identity of differentially expressed transcripts and the putative post-transcriptional regulation processes. Furthermore, a survey with a larger number of different *Daphnia galeata* clones is needed in order to unveil whether differences in the responses to kairomones can be detected. Such differences, if found, may contribute to the understanding of the role of developmental plasticity in the field of Eco-Evo-Devo.

## References
* Agrawal, A. A., Laforsch, C., & Tollrian, R. (1999). Transgenerational induction of defences in animals and plants. Nature, 401(6748), 60–63. https://doi.org/10.1038/43425
* Alexa, A., & Rahnenfuhrer, J. (2021). topGO: Enrichment Analysis for Gene Ontology (R package version 2.44.0). https://doi.org/10.18129/B9.bioc.topGO
* Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
* Brusca, R. C., Moore, W., & Shuster, S. M. (2016). Invertebrates. Sinauer Associates Inc.
* Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59–60. https://doi.org/10.1038/nmeth.3176
* Colbourne, J. K., Singan, V. R., & Gilbert, D. G. (2005). wFleaBase: The Daphnia genome database. BMC Bioinformatics, 6, 1–5. https://doi.org/10.1186/1471-2105-6-45
* Fu, L., Niu, B., Zhu, Z., Wu, S., & Li, W. (2012). CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics, 28(23), 3150–3152. https://doi.org/10.1093/bioinformatics/bts565
* Grabherr, M. G., Haas, B. J., Yassour, M., Levin, J. Z., Thompson, D. A., Amit, I., Adiconis, X., Fan, L., Raychowdhury, R., Zeng, Q., Chen, Z., Mauceli, E., Hacohen, N., Gnirke, A., Rhind, N., di Palma, F., Birren, B. W., Nusbaum, C., Lindblad-Toh, K., … Regev, A. (2011). Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nature Biotechnology, 29(7), 644–652. https://doi.org/10.1038/nbt.1883
* Grüning, B., Dale, R., Sjödin, A., Chapman, B. A., Rowe, J., Tomkins-Tinch, C. H., Valieris, R., & Köster, J. (2018). Bioconda: sustainable and comprehensive software distribution for the life sciences. Nature Methods, 15(7), 475–476. https://doi.org/10.1038/s41592-018-0046-7
* Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357–359. https://doi.org/10.1038/nmeth.1923
* Lee, B. Y., Choi, B. S., Kim, M. S., Park, J. C., Jeong, C. B., Han, J., & Lee, J. S. (2019). The genome of the freshwater water flea Daphnia magna: A potential use for freshwater molecular ecotoxicology. Aquatic Toxicology, 210(January), 69–84. https://doi.org/10.1016/j.aquatox.2019.02.009
* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
* Miyakawa, H., Imai, M., Sugimoto, N., Ishikawa, Y., Ishikawa, A., Ishigaki, H., Okada, Y., Miyazaki, S., Koshikawa, S., Cornette, R., & Miura, T. (2010). Gene up-regulation in response to predator kairomones in the water flea, Daphnia pulex. BMC Developmental Biology, 10(1), 45. https://doi.org/10.1186/1471-213X-10-45
* Nishimura, O., Hara, Y., & Kuraku, S. (2019). Evaluating Genome Assemblies and Gene Models Using gVolante. In Methods in Molecular Biology (Vol. 1962, pp. 247–256). https://doi.org/10.1007/978-1-4939-9173-0_15
* Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUSCO: Assessing Genome Assembly and Annotation Completeness. In Notes on the Greek Text of Genesis (Vol. 1962, pp. 227–245). SBL Press. https://doi.org/10.1007/978-1-4939-9173-0_14
* Supek, F., Bošnjak, M., Škunca, N., & Šmuc, T. (2011). Revigo summarizes and visualizes long lists of gene ontology terms. PLoS ONE, 6(7). https://doi.org/10.1371/journal.pone.0021800
* Tams, V., Lüneburg, J., Seddar, L., Detampel, J.-P., & Cordellier, M. (2018). Intraspecific phenotypic variation in life history traits of Daphnia galeata populations in response to fish kairomones. PeerJ, 6(10), e5746. https://doi.org/10.7717/peerj.5746
* Tams, V., Nickel, J. H., Ehring, A., & Cordellier, M. (2020). Insights into the genetic basis of predator‐induced response in Daphnia galeata. Ecology and Evolution, 10(23), 13095–13108. https://doi.org/10.1002/ece3.6899
* Tarazona, S., Furió-Tarí, P., Turrà, D., Pietro, A. Di, Nueda, M. J., Ferrer, A., & Conesa, A. (2015). Data quality aware analysis of differential expression in RNA-seq with NOISeq R/Bioc package. Nucleic Acids Research, 43(21), gkv711. https://doi.org/10.1093/nar/gkv711
* Törönen, P., Medlar, A., & Holm, L. (2018). PANNZER2: a rapid functional annotation web server. Nucleic Acids Research, 46(W1), W84–W88. https://doi.org/10.1093/nar/gky350
* Ye, Z., Xu, S., Spitze, K., Asselman, J., Jiang, X., Ackerman, M. S., Lopez, J., Harker, B., Raborn, R. T., Thomas, W. K., Ramsdell, J., Pfrender, M. E., & Lynch, M. (2017). A New Reference Genome Assembly for the Microcrustacean Daphnia pulex. G3 Genes|Genomes|Genetics, 7(5), 1405–1416. https://doi.org/10.1534/g3.116.038638
* Yousuf, D. J., Phulia, V., Bhat, I. A., Nazir, M. I., Rasool, I., Hb, R.-A., & Ds, B. (2020). Kairomones: Interspecific Chemical Signalling System in Aquatic Ecosystems. World Journal of Aquaculture Research & Development, 2(January), 1007.
