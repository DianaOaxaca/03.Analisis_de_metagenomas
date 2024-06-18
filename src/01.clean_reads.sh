#!/bin/bash
#DianaOaxaca 110923

#Remover lecturas de baja calidad, eliminar adaptadores, eliminar secuencias con Ns y remover lecturas menores a 75 pb

mkdir -p results/{02.trimgalore}

FASTQ=$(ls data/raw/*.fastq | sed 's/_.*.$//' | sed 's/data\///' |sed 's/raw\///' | sort -u)
declare -a FASTQs=("$FASTQ")

for FILE in ${FASTQs[@]}; do
        echo "Empieza la limpieza con la muestra: " $FILE
        TRIMline='trim_galore --fastqc --illumina -j 45 --paired data/raw/'$FILE'_1_10M.fastq  data/raw/'$FILE'_2_10M.fastq -o results/02.trimgalore/'$FILE'_trimgalore'
        echo "La linea de comandos es: " $TRIMline
        $TRIMline
done
