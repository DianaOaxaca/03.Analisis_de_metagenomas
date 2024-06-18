#!/usr/bin/bash
#DianaOaxaca
#Subseq random reads with seqtk

FASTQ=$(ls data/*.fastq)

for FILE in $FASTQ; do
        OUT=$(echo $FILE | sed 's/data\///g' | sed 's/.fastq//')
        CMDline='seqtk sample -s100 '$FILE' 10000000 > data/'$OUT'_10M.fastq'
        CMDrun=$(seqtk sample -s100 $FILE 10000000 > data/$OUT'_10M.fastq')
        echo -e $CMDline"\n"$CMDrun
done
