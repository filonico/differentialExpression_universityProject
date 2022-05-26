#!/bin/bash

#random subsampling paired end reads
#$1 = left fastq
#$2 = right fastq
#$3 = number of reads to keep

paste $1 $2 | awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' | shuf | head -n $3 | sed 's/\t\t/\n/g' | awk -F "\t" '{print $1 > "'$1'_subsampl.fastq"; print $2 > "'$2'_subsampl.fastq"}'
