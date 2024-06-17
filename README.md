# Análisis de datos metagenómicos

Práctica de metagenómica para el tópico **Hackeando las comunidades microbianas** impartido en el Instituto de Ciencias del Mar y Limnología (**ICMyL**) de la **UNAM**. El material de esta práctica inicialmente se preparó para los Talleres Internacionales de Bioinformática (**TIB2022**) organizados por la Sociedad Mexicana de Bioinformática (SMB) y la Comunidad de Desarrollo de Software Bioinformático (**CDSB**). El material se ha ido modificando con el tiempo y sigue en desarrollo. 

Autoras: Mirna Rosas Vázquez Landa [@MirnaVazquez](https://github.com/mirnavazquez) y Diana Hernández Oaxaca [@DianaOaxaca](https://github.com/DianaOaxaca)



#### **Datos de la fermentación del pulque**

**Sinopsis**: En este trabajo se secuenciaron librerias de Illumina MiSeq (2 X 149 pb) de muestras a diferentes tiempos de fermentación: Aguamiel, 0, 3, 6 y 12 horas. En general observaron que la abundancia de los géneros cambió durante la fermentación y se asoció con una disminución de sacarosa, aumento de etanol y ácido láctico (Medidos por HPLC). Predijeron, entre otros, la biosíntesis de folato y vitamina B. Uno de los MAG obtenidos correspondió a Saccharomyces cereviseae relacionado filogenéticamente con S. cerevisiae aislado de sake y bioetanol y otro correspondió a Zymomonas mobilis que propusieron como un nuevo linaje.

Trabajaremos con un subset de los datos originales.

**Los datos crudos**: PRJNA603591
https://www.ebi.ac.uk/ena/browser/view/PRJNA603591

**El articulo**: Genomic profiling of bacterial and fungal communities and their predictive functionality during pulque fermentation by whole-genome shotgun sequencing
https://www.nature.com/articles/s41598-020-71864-4

**Referencia:** Chacón-Vargas, K., Torres, J., Giles-Gómez, M. et al. Genomic profiling of bacterial and fungal communities and their predictive functionality during pulque fermentation by whole-genome shotgun sequencing. Sci Rep 10, 15115 (2020). https://doi.org/10.1038/s41598-020-71864-4

### 0. Espacio de trabajo

```bash
# Creamos los directorios en donde empezaremos a trabajar
mkdir -p Análisis_de_metagenomas/{data/{raw,clean},results}
```

```bash
#Accedemos al directorio principal
cd Análisis_de_metagenomas
```

### 1. Análisis de calidad y filtrado de las lecturas

Una herramienta común para verificar la calidad de las lecturas es **[FastQC](https://github.com/s-andrews/FastQC),** genera reportes que nos permiten darnos una idea de como están las lecturas, los posibles problemas y que futuros análisis son necesarios.

Para obtener los reportes y el filtrado de calidad ejecutamos la siguiente linea:

```bash
mkdir -p results/{01.fastqc,02.trimgalore}
```

```bash
cd data/raw
ln -s ../../curso_metagenomas/pulque/ .
```

Activamos el ambiente patra **fastqc**

```bash
conda activate fastqc
```

Ejecutamos  fastqc

```bash
fastqc data/raw/pulque/*.fastq -o results/01.fastqc/
```

Desactivamos el ambiente

```bash
conda deactivate
```

Ahora activamos el ambiente donde se encuentra **multiqc**

```bash
conda activate metagenomas
```

```bash

multiqc results/01.fastqc/*.zip -o results/01.fastqc/multiqc
```

La herramienta **[TrimGalore](https://github.com/FelixKrueger/TrimGalore)** nos permite eliminar lecturas de baja calidad, adaptadores, etc. Y con  [**MultiQC**](https://github.com/ewels/MultiQC), podemos ver las calidades del conjunto de lecturas. Existen otros programas para limpiar las lecturas como **[Trimmomatic](https://github.com/usadellab/Trimmomatic)**.

```bash
#Ejecutamos trim_galore para filtrar
trim_galore --fastqc -j 45 --paired data/raw/pulquet0_1_10M.fastq data/raw/pulquet0_2_10M.fastq -o results/02.trimgalore/pulquet0_trimgalore
```

Para evaluar los resultados de calidad tras el filtrado, podemos ejecutar **multiqc**.

Creamos el directorio donde alojaremos las entradas de multiqc y movemos estos archivos al nuevo directorio.

```bash
mkdir -p results/02.trimgalore/zips
mv results/02.trimgalore/*/*.zip results/02.trimgalore/zips/
```

Ahora si ejecutamos multiqc que esta en en el ambiente de qiime2, por eso lo vamos a activar. OJO, multiqc puede estar alojado en un ambiente independiente o instalarse fuera de ambientes, en el servidor esta dentro de qiime pero no es una regla.

```bash
#Activamos el ambiente
#conda activate /botete/diana/.conda/envs/qiime2-2022.11
#Ejecutamos multiqc
multiqc results/02.trimgalore/zips/*.zip -o results/02.trimgalore/zips/multiqc 
```

Movemos nuestros datos limpios al subdirectorio `data/clean`

```bash
mv results/02.trimgalore/*/*.fq data/clean/
```

Y descativamos el ambiente

```bash
conda deactivate
```

### 2. Ensamble metagenómico

En este emplo ensamblaremos todas las muestras que dan lugar a la fermentación del pulque, es decir, haremos un **coensamble**. Sin embargo, también podríamos correr un esamble por muestra de manera independiente. 

Ensamblaremos con **[MetaSPAdes](https://github.com/ablab/spades)** y **[Megahit](https://github.com/voutcn/megahit)**, después compararemos ambos coensambles basándonos en longitudes de los contigs con QUAST y el porcentaje de lecturas que dan lugar al coensamble, esto con bbmap de **[bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)**. 

Para acceder a metaspades y megahit activamos el ambiente metagenomics

```bash
conda activate metagenomics
```

#### **[Megahit](https://github.com/voutcn/megahit)**

Ejecutamos megahit ...

```bash
megahit -t 34 --k-list 21,33,55,77,99,127 --min-contig-len 1000 -1 data/clean/fermentation_1.fastq -2 data/clean/fermentation_2.fastq -o results/03.megahit 
```

#### **[MetaSPAdes](https://github.com/ablab/spades)** 

Para SPAdes necesitamos crear el directorio en el que se alojarán los resultados

```bash
mkdir -p results/03.metaspades
```

Y lo ejecutamos ....

```bash
spades.py --meta -k 21,33,55,77,99,127 -t 30 -1 data/clean/fermentation_1.fastq -2 data/clean/fermentation_2.fastq -o results/03.metaspades 
```

Desactivamos el ambiente

```bash
conda deactivate
```

En spades no pudimos delimitar el tamaño mínimo de contig, por lo tanto tenemos que eliminar los de menor tamaño.

Veamos cuanto miden los contigs de nuestro ensamble.

```bash
infoseq -auto -nocolumns -delimiter '\t' -only -noheading -name -length results/03.metaspades/contigs.fasta | sort -k2n | cut -f2 | sort | uniq -c | sort -k2n | less
```

Vamos a filtrarlo ... Para el filtrado necesitamos entrar al directorio de resultados de spades y ejecutar el script desde ahí.

```bash
python filter_contig_length_sp.py 1000 1 contigs.fasta
```

Comparemos el número de contigs que había antes y después de filtrar.

```bash
grep -c '>' contigs.fasta 
173273
grep -c '>' contigs.fltr.fasta 
18673
```

### 3. Comparación de los ensambles

Creamos el directorio para los mapeos

```bash
mkdir -p results/04.depth/
```

####  **[bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)**

Mapeamos las lecturas al ensamble de metaspades con **bbmap**

```bash
#Hacer ligas simbolicas y correr dentro del directorio 04.depth
nohup bbmap.sh ref=results/03.metaspades/contigsftr.fasta in=data/clean/fermentation_1.fastq in2=data/clean/fermentation_2.fastq out=results/04.depth/fermentation.metaspades.sam   scafstats=fermentation_metaspades.scafstats &
```

Bien, veamos cuantas lecturas mapearon al ensamble

```bash
#Número de lecturas iniciales
95588026
#Obtenemos las que mapearon
cut -f8 results/04.depth/fermentation_metaspades.scafstats | grep -v assignedReads | awk '{sum += $1;}END{print sum;}'
91148282
#Sacamos el porcentaje
95.35 %
```

Hagamos lo mismo para el de megahit

```bash
nohup bbmap.sh ref=results/03.megahit/final.contigs.fa in=data/clean/fermentation_1.fastq in2=data/clean/fermentation_2.fastq out=results/04.depth/fermentation.megahit.sam maxindel=80 scafstats=results/04.depth/fermentation_megahit.scafstats &
```

Cuántas mapearon?

```bash
95588026
cut -f8 results/04.depth/fermentation_megahit.scafstats | grep -v assignedReads | awk '{sum += $1;}END{print sum;}'
91218466
95.428
```

Ok, falta observar las longitudes de los contigs y el tamaño del ensamble

#### [QUAST](https://github.com/ablab/quast)

Activemos el ambiente metagenomics para ejecutar **QUAST**

```bash
conda activate metagenomics
```

Ahora ejecutemos **QUAST** ...

```bash
quast.py results/03.metaspades/contigsftr.fasta results/03.megahit/final.contigs.fa -o results/05.quast
```

Desactivemos el ambiente

```bash
conda deactivate
```

Y veamos la salida de QUAST para determinar la mejor opción.

### 4. Binning

Ahora que sabemos que ensamble usar para el binning, empezemos.

Primero vamos a borrar el archivo .sam de megahit porque es muy muy pesado y no lo necesitamos.

```bash
rm results/04.depth/fermentation.megahit.sam
```

Ahora si, para poder agrupar en *bins* el coensamble, los *binneadores* requieren de un archivo con las coverturas de cada contig. Afortunadamente ya habíamos mapeado las lecturas, asi que ya tenemos un archivo de mapeo, nos falta convertirlo a formato bam y ordenarlo. Para ello ejecutamos **[samtools](http://www.htslib.org/doc/samtools.html)**

 

```bash
samtools view -bShu results/04.depth/fermentation.metaspades.sam | samtools sort -@ 80 -o results/04.depth/fermentation_metaspades_sorted.bam
samtools index results/04.depth/fermentation_metaspades_sorted.bam
```

Y borramos el `.sam`

```bash
rm results/04.depth/fermentation.metaspades.sam
```

Ahora vamos a obtener la información en el formato que necesitan los programas que harán el *binning*. Para ello usamos `jgi_summarize_bam_contig_depths` .

Activamos el ambiente de metabat

```bash
conda activate metabat
```

Ahora si, obtengamos la información del mapeo

```bash
jgi_summarize_bam_contig_depths --outputDepth results/04.depth/fermentation_metaspades_depth.txt results/04.depth/fermentation_metaspades_sorted.bam 
```

#### [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)

Ejecutemos **metabat2**

```bash
metabat2 -i results/03.metaspades/contigsftr.fasta -a results/04.depth/fermentation_metaspades_depth.txt -o results/06.metabat/metabat_bin -t 50 -m 1500
# Y desactivemos el ambiente
conda deactivate
```

#### [Maxbin](https://sourceforge.net/p/maxbin/code/ci/master/tree/)

Ahora ejecutemos **Maxbin**

Activa el ambiente

```bash
conda activate maxbin
```

Creamos el directorio de estos resultados

```bash
mkdir -p results/07.maxbin/
```

Ahora si, corremos Maxbin

```bash
run_MaxBin.pl -thread 48 -min_contig_length 1500 -contig results/03.metaspades/contigsftr.fasta -out results/07.maxbin/maxbin  -abund results/04.depth/fermentation_metaspades_depth.txt
```

Y desactivamos el ambiente

```bash
conda deactivate
```

Ya casi....

#### [Vamb](https://github.com/RasmussenLab/vamb)

Corramos **Vamb**

```bash
/botete/diana/.local/bin/vamb
vamb --fasta results/03.metaspades/contigsftr.fasta --jgi results/04.depth/fermentation_metaspades_depth.txt --minfasta 500000 --outdir results/08.vamb
```

### 5. Refinamiento

#### [DAS_Tool](https://github.com/cmks/DAS_Tool)

Tenemos los resultados de tres *binneadores*, muchos de estos *bins* estarán repetidos, ejecutaremos DAS_Tool para desreplicar estos bins, el flujo para DAS_Tool es el siguiente:

Primero activamos el ambiente.

```bash
conda activate das_tool
```

Creamos un directorio de trabajo para los resultados

```bash
mkdir -p results/09.dastool
```

Y creamos archivos tabulares que son legibles para DAS_Tool

```bash
#Para metabat
Fasta_to_Contig2Bin.sh -i results/06.metabat/ -e fa > results/09.dastool/fermentation_metabat.dastool.tsv
#Para maxbin
Fasta_to_Contig2Bin.sh -i results/07.maxbin/ -e fasta > results/09.dastool/fermentation_maxbin.dastool.tsv
#Para vamb
Fasta_to_Contig2Bin.sh -i results/08.vamb/bins/ -e fna > results/09.dastool/fermentation_vamb.dastool.tsv
```

Ahora si, ejecutemos DAS_Tool:

```bash
DAS_Tool -i results/09.dastool/fermentation_maxbin.dastool.tsv,results/09.dastool/fermentation_metabat.dastool.tsv,results/09.dastool/fermentation_vamb.dastool.tsv -l maxbin,metabat,vamb -c results/03.metaspades/contigsftr.fasta -o results/09.dastool/fermentation_bins -t 12 --write_bins
```

#### [CheckM](https://ecogenomics.github.io/CheckM/)

Ya desreplicamos los *bins* que obtuvimos, pero nos falta evaluar la calidad de estos, para ello ejecutaremos **CheckM**

```bash
checkm lineage_wf results/09.dastool/fermentation_bins_DASTool_bins/ results/10.checkm/ -x fa -t 40  -f results/10.checkm/checkm_dastool_bins.txt
```



Recordemos que hay más herramientas para desreplicar como **dRep** y también para refinar los bins que obtengamos, como **RefineM**. Si te es posible, detente a observar tus datos, analiza los resultados, prueba lo que puedas y toma las mejores decisiones.



## Ejercicio!!!

Bien, ya hemos visto el flujo de trabajo, revisamos algunas salidas y las tomas de decisiones. Ahora te toca a tí.

¿Cómo serían los resultados si en lugar de un coensamble se analizan las muestras individuales? Averiguémoslo entre todXs.

En equipos, ejecuten el flujo de análisis a partir de los ensamble y su comparación hasta el binning. Tomen sólo las lecturas de una muestra, anotaremos los resultados en una tabla compartida y al final los discutiremos.

1. Elige un equipo y que muestra analizarán

2. Crea los directorios de trabajo, el directorio principal será el del nombre d ela muestra que analizarán

3. Crea una liga simbólica de los datos limpios correspondientes a la muestra que eligieron y ponlos en tu directorio `data` los datos se encuentran en:

   ```bash
   /botete/diana/Hackeando_las_comunidades_microbianas_v1/03.Analisis_de_metagenomas/data/clean/
   ```

4. Para remover los contigs cortos puedes copiar el script `filter_contig_length_sp.py` al directorio de resultados de **metaspades**, su ubicación es:

   ```bash
   /botete/diana/Hackeando_las_comunidades_microbianas_v1/03.Analisis_de_metagenomas/src/filter_contig_length_sp.py
   ```

   

