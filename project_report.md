# Gene transcription responses of *Daphnia galeata* to the prolonged presence of natural predators
## Introduction
*Daphnia* is a genus of branchiopod crustaceans belonging to the sub-order Cladocera (order Diplotraca) and comprises the so-called water fleas. Daphniids exhibit a folded (‘bivalve’) carapace that encloses the entire trunk (but not the cephalon as it does in other diplostracans) and serves mainly as a brood chamber (Brusca et al., 2016). Water fleas are distributed worldwide in freshwater basins and ponds where they play a key role as bioindicators of environmental health and quality (Miyakawa et al., 2010 and references therein). In addition, *Daphnia pulex* is the first crustacean sensu stricto to have had its genome sequenced and, together with other congeneric species, now counts on a personal genome web-database, the [wFleaBase](http://wfleabase.org/@blank) (Colbourne et al., 2005).
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

Raw RNA-sequencing reads of *Daphnia galeata* (BioProject: [PRJEB29887](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB29887/@blank)) were downloaded from NCBI SRA. **Prefetch v2.8.0** and **FASTQ-dump v2.8.0** from a Bioconda (Grüning et al., 2018) installation of the NCBI SRA-toolkit were used to retrieve SRA and FASTQ files, respectively. FASTQ files were given a header syntax suitable for the subsequent Trinity assembly (<code>--defline-seq '@$sn[_$rn]/$ri'</code>) and were split to keep paired reads separated (<code>--split-files</code>). **vdb-validate v2.8.0** was used to analyze the downloaded SRA data for corruption and other problems, while **FastQC v0.11.7** was used to check read quality. Only 6 runs out of 12 were eventually downloaded: 3 (ERR2929117, ERR2929119 and ERR2929120) were obtained from a *Daphnia galeata* clone bred in absence of fish kairomones (control); conversely, 3 (ERR2929122, ERR2929124, ERR2929126) were obtained from a *Daphnia galeata* clone bred in presence of ide (*Leucuscus idus*) kairomones. The whole process of downloading and quality checks was auto-mated as following:

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

The <code>random_subsampling_PE.sh</code> bash script can be accessed [here](random_subsampling_PE.sh/@blank).
