# Análisis de datos metagenómicos

Práctica para el taller de metagenómica de la Red Mexicana de Bioinformática **RMB** organizado en el **CIBNOR Junio 2024**

El material de esta práctica inicialmente se preparó para los Talleres Internacionales de Bioinformática (**TIB2022**) organizados por la Sociedad Mexicana de Bioinformática (SMB) y la Comunidad de Desarrollo de Software Bioinformático (**CDSB**). Se ha usado en talleres intersemestrales del Instituto de Ciencias del Mar y Limnología (**ICMyL**) de la **UNAM** y el tópico de posgrado **Hackeando las comunidades microbianas**. El material se ha ido modificando con el tiempo y sigue en desarrollo. 

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
ln -s ../../../curso_metagenomas/pulque/ .
```

Activamos el ambiente patra **fastqc**

```bash
conda activate fastqc
```

Ejecutamos  fastqc

```bash
cd ~/Análisis_metagenomas/
```

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
Y corremos multiqc

```bash

multiqc results/01.fastqc/*.zip -o results/01.fastqc/multiqc
```
Descarga y visualiza el archivo multiqc_report.html

En una terminal nueva sin conectarse al servidor, es decir que sea la terminal de tu computadora

```bash
scp USUARIO@IP:/home/USUARIO/Análisis_de_metagenomas/results/01.fastqc/multiqc/multiqc_report.html .
```
Te pedirá tu contrase;a del servidor


La herramienta **[TrimGalore](https://github.com/FelixKrueger/TrimGalore)** nos permite eliminar lecturas de baja calidad, adaptadores, etc. Y con  [**MultiQC**](https://github.com/ewels/MultiQC), podemos ver las calidades del conjunto de lecturas. Existen otros programas para limpiar las lecturas como **[Trimmomatic](https://github.com/usadellab/Trimmomatic)**.

Vamos a limpiar los datos de **Aguamiel**

```bash
trim_galore --illumina --fastqc -j 15 --paired data/raw/pulque/Pulque-AM_SRR10997050_1_10M.fastq data/raw/pulque/Pulque-AM_SRR10997050_2_10M.fastq -o results/02.trimgalore/pulqueAM_trimgalore
```

**NOTA** Limpiamos los de aguamiel, pero faltan otros tiempos. Podemos hacer un ciclo for para correr todos o correrlos por separado.

Para evaluar los resultados de calidad tras el filtrado, podemos ejecutar **multiqc** otra vez pero sobre los datos ya limpios


Ejecutamos multiqc

```bash
multiqc results/02.trimgalore/pulqueAM_trimgalore/*.zip -o results/02.trimgalore/multiqc 
```

Movemos nuestros datos limpios al subdirectorio `data/clean`


```bash
mv results/02.trimgalore/*/*.fq data/clean/
```

### 2. Ensamble metagenómico

En este emplo ensamblaremos todas las muestras que dan lugar a la fermentación del pulque, es decir, haremos un **coensamble**. Sin embargo, también podríamos correr un esamble por muestra de manera independiente. 

Ensamblaremos con **[MetaSPAdes](https://github.com/ablab/spades)** y **[Megahit](https://github.com/voutcn/megahit)**, después compararemos ambos coensambles basándonos en longitudes de los contigs con QUAST y el porcentaje de lecturas que dan lugar al coensamble, esto con bbmap de **[bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)**. 

Para este coensamble concatenamos todas las librerias limpias con `cat` y lo nombramos fermentation **OJO** aqui ya tenemos todas las librerias limpias y copiadas en el directorio `data/clean`

```bash
cat results/02.trimgalore/*/*_1.fq > data/clean/fermentation_1.fastq
cat results/02.trimgalore/*/*_2.fq > data/clean/fermentation_2.fastq
```

#### **[Megahit](https://github.com/voutcn/megahit)**

Ejecutamos megahit ...

**Coensamble**
```bash
megahit -t 12 --k-list 21,33,55,77,99,127 --min-contig-len 1000 -1 data/clean/fermentation_1.fastq -2 data/clean/fermentation_2.fastq -o results/03.megahit 
```
**Solo aguamiel**

```bash
nohup megahit -t 12 --k-list 21,33,55,77,99,127 --min-contig-len 1000 -1 data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq -2 data/clean/Pulque-AM_SRR10997050_2_10M_val_2_val_2.fq -o results/03.megahit_AM &
```

#### **[MetaSPAdes](https://github.com/ablab/spades)** 

Para SPAdes necesitamos crear el directorio en el que se alojarán los resultados

**coensamble**
```bash
mkdir -p results/03.metaspades
```

Y lo ejecutamos ....
**coensamble**

```bash
 spades.py --meta -k 21,33,55,77,99,127 -t 12 -1 data/clean/fermentation_1.fastq -2 data/clean/fermentation_2.fastq -o results/03.metaspades &
```
**solo aguamiel**

```bash
nohup spades.py --meta -k 21,33,55,77,99,127 -t 12 -1 data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq -2 data/clean/Pulque-AM_SRR10997050_2_10M_val_2_val_2.fq -o results/03.metaspades_AM &
```


En spades no pudimos delimitar el tamaño mínimo de contig, por lo tanto tenemos que eliminar los de menor tamaño.

Veamos cuanto miden los contigs de nuestro ensamble.

```bash
infoseq -auto -nocolumns -delimiter '\t' -only -noheading -name -length results/03.metaspades_AM/contigs.fasta | sort -k2n | cut -f2 | sort | uniq -c | sort -k2n | less
```

Vamos a filtrarlo ... Para el filtrado necesitamos entrar al directorio de resultados de spades y ejecutar el script desde ahí.

```bash
cd results/03.metaspades_AM/
```

Copia y pega el contenido del script `src/filter_contig_length_sp.py`

```bash
python filter_contig_length_sp.py 1000 1 contigs.fasta
```

Comparemos el número de contigs que había antes y después de filtrar.

```bash
grep -c '>' contigs.fasta 

grep -c '>' contigs.fltr.fasta 

```

### 3. Comparación de los ensambles

Creamos el directorio para los mapeos

```bash
mkdir -p results/04.depth/
```

####  **[bbtools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/)**

Mapeamos las lecturas al ensamble de metaspades con **bbmap**

```bash
nohup bbmap.sh ref=results/03.metaspades/contigsftr.fasta in=data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq in2=data/clean/Pulque-AM_SRR10997050_2_10M_val_2_val_2.fq out=results/04.depth/Pulque-AM.metaspades.sam   scafstats=Pulque-AM_metaspades.scafstats > nohup_bbmap_metaspades.out &
```

Bien, veamos cuantas lecturas mapearon al ensamble

Número de lecturas iniciales

```bash
wc -l data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq
```
**NOTA** este número entre 4 para sacar el numero acercado de lecturas 

Obtenemos las que mapearon

```bash
cut -f8 results/04.depth/Pulque-AM_metaspades.scafstats | grep -v assignedReads | awk '{sum += $1;}END{print sum;}'
```
#Sacamos el porcentaje


Hagamos lo mismo para el de megahit

```bash
nohup bbmap.sh ref=results/03.megahit/final.contigs.fa in=data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq in2=data/clean/Pulque-AM_SRR10997050_2_10M_val_2_val_2.fq out=results/04.depth/Pulque-AM.megahit.sam maxindel=80 scafstats=results/04.depth/Pulque-AM_megahit.scafstats > nohup_bbmap_megahit.out &
```

Cuántas mapearon?

```bash
cut -f8 results/04.depth/Pulque-AM_megahit.scafstats | grep -v assignedReads | awk '{sum += $1;}END{print sum;}'
```


Ok, falta observar las longitudes de los contigs y el tamaño del ensamble

#### [QUAST](https://github.com/ablab/quast)

Activemos el ambiente metagenomas para ejecutar **QUAST**

```bash
conda activate metagenomas
```

Ahora ejecutemos **QUAST** ...

```bash
quast.py results/03.metaspades/contigsftr.fasta results/03.megahit/final.contigs.fa -o results/05.quast
```


Y veamos la salida de QUAST para determinar la mejor opción.

### 4. Binning

Ahora que sabemos que ensamble usar para el binning, empezemos.

Primero vamos a borrar el archivo .sam de megahit porque es muy muy pesado y no lo necesitamos.

```bash
rm results/04.depth/Pulque-AM.megahit.sam
```

Ahora si, para poder agrupar en *bins* el coensamble, los *binneadores* requieren de un archivo con las coverturas de cada contig. Afortunadamente ya habíamos mapeado las lecturas, asi que ya tenemos un archivo de mapeo, nos falta convertirlo a formato bam y ordenarlo. Para ello ejecutamos **[samtools](http://www.htslib.org/doc/samtools.html)**

 

```bash
samtools view -bShu results/04.depth/Pulque-AM.metaspades.sam | samtools sort -@ 80 -o results/04.depth/Pulque-AM_metaspades_sorted.bam
samtools index results/04.depth/Pulque-AM_metaspades_sorted.bam
```

Y borramos el `.sam`

```bash
rm results/04.depth/Pulque-AM.metaspades.sam
```

Ahora vamos a obtener la información en el formato que necesitan los programas que harán el *binning*. Para ello usamos `jgi_summarize_bam_contig_depths` .


Ahora si, obtengamos la información del mapeo

```bash
jgi_summarize_bam_contig_depths --outputDepth results/04.depth/Pulque-AM_metaspades_depth.txt results/04.depth/Pulque-AM_metaspades_sorted.bam 
```

#### [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)

Ejecutemos **metabat2**

```bash
metabat2 -i results/03.metaspades/contigsftr.fasta -a results/04.depth/Pulque-AM_metaspades_depth.txt -o results/06.metabat/metabat_bin -t 50 -m 1500
```

#### [Maxbin](https://sourceforge.net/p/maxbin/code/ci/master/tree/)

Ahora ejecutemos **Maxbin**


Creamos el directorio de estos resultados

```bash
mkdir -p results/07.maxbin/
```

Ahora si, corremos Maxbin

```bash
run_MaxBin.pl -thread 48 -min_contig_length 1500 -contig results/03.metaspades/contigsftr.fasta -out results/07.maxbin/maxbin  -abund results/04.depth/Pulque-AM_metaspades_depth.txt
```

Y desactivamos el ambiente

```bash
conda deactivate
```

Ya casi....

#### [Vamb](https://github.com/RasmussenLab/vamb)

Corramos **Vamb**

Activamos el ambiente de vamb

```bash
conda activate vamb_env
```

```bash
vamb --fasta results/03.metaspades/contigsftr.fasta --jgi results/04.depth/Pulque-AM_metaspades_depth.txt --minfasta 500000 --outdir results/08.vamb
```

Descativamos el ambiente

```bash
conda deactivate
```

### 5. Refinamiento

#### [DAS_Tool](https://github.com/cmks/DAS_Tool)

Tenemos los resultados de tres *binneadores*, muchos de estos *bins* estarán repetidos, ejecutaremos DAS_Tool para desreplicar estos bins, el flujo para DAS_Tool es el siguiente:

Primero activamos el ambiente.

```bash
conda activate das_tool_env
```

Creamos un directorio de trabajo para los resultados

```bash
mkdir -p results/09.dastool
```

Y creamos archivos tabulares que son legibles para DAS_Tool

```bash
#Para metabat
Fasta_to_Contig2Bin.sh -i results/06.metabat/ -e fa > results/09.dastool/Pulque-AM_metabat.dastool.tsv
#Para maxbin
Fasta_to_Contig2Bin.sh -i results/07.maxbin/ -e fasta > results/09.dastool/Pulque-AM_maxbin.dastool.tsv
#Para vamb
Fasta_to_Contig2Bin.sh -i results/08.vamb/bins/ -e fna > results/09.dastool/Pulque-AM_vamb.dastool.tsv
```

Ahora si, ejecutemos DAS_Tool:

```bash
DAS_Tool -i results/09.dastool/Pulque-AM_maxbin.dastool.tsv,results/09.dastool/Pulque-AM_metabat.dastool.tsv,results/09.dastool/Pulque-AM_vamb.dastool.tsv -l maxbin,metabat,vamb -c results/03.metaspades/contigsftr.fasta -o results/09.dastool/Pulque-AM_bins -t 12 --write_bins
```

Desactivamos el ambiente

```bash
conda deactivate
```

#### [CheckM](https://ecogenomics.github.io/CheckM/)

Ya desreplicamos los *bins* que obtuvimos, pero nos falta evaluar la calidad de estos, para ello ejecutaremos **CheckM**

Activamos el ambiente

```bash
conda activate metagenomas
```

```bash
checkm lineage_wf results/09.dastool/Pulque-AM_bins_DASTool_bins/ results/10.checkm/ -x fa -t 40  -f results/10.checkm/checkm_dastool_bins.txt
```



Recordemos que hay más herramientas para desreplicar como **dRep** y también para refinar los bins que obtengamos, como **RefineM**. Si te es posible, detente a observar tus datos, analiza los resultados, prueba lo que puedas y toma las mejores decisiones.



## Ejercicio!!!

Bien, ya hemos visto el flujo de trabajo, revisamos algunas salidas y las tomas de decisiones. Ahora te toca a tí.

¿Cómo serían los resultados si en lugar de un coensamble se analizan las muestras individuales? Averiguémoslo entre todXs.

En equipos, ejecuten el flujo de análisis a partir de los ensamble y su comparación hasta el binning. Tomen sólo las lecturas de una muestra, anotaremos los resultados en una tabla compartida y al final los discutiremos.

1. Elige un equipo y que muestra analizarán

2. Crea los directorios de trabajo, el directorio principal será el del nombre d ela muestra que analizarán

3. Crea una liga simbólica de los datos limpios correspondientes a la muestra que eligieron y ponlos en tu directorio `data`

   

