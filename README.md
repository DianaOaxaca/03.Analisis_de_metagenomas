**Link para ver este contenido en html**

[taller_metagenomas_MAGs](https://dianaoaxaca.github.io/03.Analisis_de_metagenomas/Taller_metagenomas_MAGs.html)


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
mkdir -p Analisis_de_metagenomas/{data/{raw,clean},results}
```

```bash
#Accedemos al directorio principal
cd Analisis_de_metagenomas
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
cd ~/Analisis_metagenomas/
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


**aguamiel**

```bash
nohup megahit -t 12 --k-list 21,33,55,77,99,127 --min-contig-len 1000 -1 data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq -2 data/clean/Pulque-AM_SRR10997050_2_10M_val_2_val_2.fq -o results/03.megahit_AM &
```

#### **[MetaSPAdes](https://github.com/ablab/spades)** 

Para SPAdes necesitamos crear el directorio en el que se alojarán los resultados

```bash
mkdir -p results/03.metaspades
```

**aguamiel**

Activar el ambiente

```bash
conda activate spades_env
```

```bash
nohup spades.py --meta -k 21,33,55,77,99,127 -t 12 -1 data/clean/Pulque-AM_SRR10997050_1_10M_val_1_val_1.fq -2 data/clean/Pulque-AM_SRR10997050_2_10M_val_2_val_2.fq -o results/03.metaspades_AM &
```

Desactivamos el ambiente

```bash
conda deactivate
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

### 5. Refinamiento y Selección

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

Analicemos la salida...

Ahora debemos filtrar los bins que cumplan con criterios de calidad

**Paso 1: remover las lineas que no nos sirven**

```bash
sed -e '1,3d' results/10.checkm/checkm_dastool_bins.txt | sed -e '6d' > results/10.checkm/CheckM-DAS_Tool_bins_mod.txt
```
**Paso 2: seleccionarlos**

Con R

```R
library(tidyverse)
# CheckM -------------------------------------------------------------------####
checkm<-read.table("CheckM-DAS_Tool_bins_mod.txt", sep = "", header = F, na.strings ="", stringsAsFactors= F)
# Extracting good quality bins Megahit ------------------------------------####
colnames(checkm)<-c("Bin_Id", "Marker", "lineage", "Number_of_genomes", 
                         "Number_of_markers", "Number_of_marker_sets", 
                         "0", "1", "2", "3", "4", "5", "Completeness", 
                         "Contamination", "Strain_heterogeneity")  

good_bins<-checkm %>%
  select(Bin_Id, Marker, Completeness, Contamination) %>%
  filter(Completeness >= 50.00 & Contamination <= 10.00)

bins<-good_bins$Bin_Id

write.table(bins, "lista_medium_bins", quote = F, row.names = F, col.names = F)
```

**Paso 3: Copiarlos**

```bash
mkdir -p results/11.Bins
```

```bash
sed 's#bin#cp results/09.dastool/Pulque-AM_bins/bin#g' results/10.checkm/lista_medium_bins | sed 's#$#.fa results/11.Bins#g' > results/11.Bins/copy_bins.sh
```

```bash
bash results/11.Bins/copy_bins.sh
```

**SELECCIÓN DE LOS BINS**
**OPCIÓN CON BASH**


Crea un script con vim o nano que se llame `filter_bins.sh` y copia el siguiente contenido:

```bash
#!/bin/bash

# Archivos y directorios
input_file="results/10.checkm/checkm_dastool_bins.txt"
source_dir="results/09.dastool/Pulque-AM_bins"
mkdir -p "results/11.Bins"
destination_dir="results/11.Bins"


# Leer el archivo línea por línea
while IFS= read -r line; do
  # Saltar la línea de encabezado y las líneas de separación
  if [[ "$line" == *"Completeness"* ]] || [[ "$line" == *"----"* ]]; then
    continue
  fi

  # Extraer los valores necesarios
  bin_id=$(echo "$line" | awk '{print $1}')
  completeness=$(echo "$line" | awk '{print $(NF-2)}')
  contamination=$(echo "$line" | awk '{print $(NF-1)}')

  # Verificar los criterios
  if (( $(echo "$completeness > 50" | bc -l) )) && (( $(echo "$contamination < 10" | bc -l) )); then
    # Formar el nombre del archivo con la extensión correcta
    file_name="${bin_id}.fa.gz"
    
    # Verificar si el archivo existe en el directorio fuente
    if [[ -f "$source_dir/$file_name" ]]; then
      # Copiar el archivo al directorio de destino
      cp "$source_dir/$file_name" "$destination_dir/"
    else
      echo "Archivo no encontrado: $file_name"
    fi
  fi
done < "$input_file"
```
 Ahora ejecutalo desde el directorio principal asi:

```bash
bash filter_bins.sh
```

Taran!!!

Ahora si, exploremos....

```bash
grep '>' results/11.Bins/*.fa
```

### 6. Taxonomía

En esta sección, exploraremos la taxonomía de los bins utilizando [GTDB-tk](https://github.com/Ecogenomics/GTDBTk).

```bash
mkdir -p results/12.GTDBtk
```
Activamos el ambiente

```bash
conda activate 	gtdbtk_env
```

Ahora si, a correrlo!!

```bash
gtdbtk classify_wf --genome_dir results/11.Bins --out_dir results/11.GTDBtk --cpus 15 -x fa
```

Después de ejecutar GTDB-tk, continuaremos en R para visualizar los datos.

Análisis en R

```R
library(tidyverse)

GTDBK <- read.table("results/11.GTDBtk/gtdbtk.bac120.summary.tsv", 
  sep = "\t", header = TRUE, na.strings = "", stringsAsFactors = FALSE) %>%
  as_tibble()

#El archivo contiene información sobre la clasificación taxonómica de los bins.
#Continuamos limpiando y transformando los datos:

pulque_gtdbtk <- GTDBK %>%
  select(user_genome, classification) %>%
  separate(classification, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  rename(Bin_name = user_genome) %>%
  unite(Bin_name_2, c("Bin_name", "Phylum"), remove = FALSE) %>%
  select(Bin_name, Domain, Phylum, Class, Order, Family, Genus, Species)

#Guardamos los datos en un archivo de metadatos:
write.table(pulque_gtdbtk, file = "results/11.GTDBTK/Metadatos.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
```

Visualización de Datos
Creamos un gráfico de barras interactivo que muestra la distribución taxonómica de los bins:

```R
GTDBtk <- pulque_gtdbtk %>%
  count(Domain, Phylum) %>%
  rename(Number_of_MAGs = n) %>%
  ggplot(aes(x = Domain, y = Number_of_MAGs, fill = Phylum)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal()
```

Si prefieres una visualización interactiva:

```R
library(plotly)
GTDBtk_p_fig <- ggplotly(GTDBtk)
```

### 7. Inferencia Metabólica
En esta sección, exploraremos la inferencia metabólica de los bins de nuestro análisis.

**Preparación de Carpetas**
Comenzaremos creando dos carpetas llamadas “Genoma” y “Proteoma” para organizar nuestros archivos:

```bash
mkdir results/11.Bins/Genoma
mkdir results/11.Bins/Proteoma
```

Luego, moveremos todos los archivos con extensión “.fa” a la carpeta “Genoma”:

```bash
mv results/11.Bins/*.fa results/11.Bins/Genoma
```

**Predicción de Proteínas**

Utilizaremos Prodigal para predecir las proteínas a partir de los genomas:

```bash
cd results/11.Bins/Genoma/
for i in $(ls *.fa); do prodigal -i $i -o results/11.Bins/Proteoma/$i.txt -a results/11.Bins/Proteoma/$i.faa ; done
cd
```

Podemos revisar las proteínas predichas con:

```bash
grep ">" results/11.Bins/Proteoma/*.faa
```

**Anotación de Proteínas con KEGG**
Utilizaremos KEGG, una base de datos que proporciona información sobre genes, proteínas, rutas metabólicas y más.

Ahora, utilizaremos [kofam_scan](https://github.com/takaram/kofam_scan) para anotar las proteínas:

Activamos el ambiente

```bash
conda activate kofamscan_env
```

```bash
mkdir -p results/13.KOFAM
cd results/11.Bins/Proteoma/
for i in *.faa; do /RUTA_DE_KOFAM/exec_annotation -o results/13.KOFAM/$i.txt $i --report-unannotated --cpu 14 -p /RUTA_databases/kofamscan_dbs/profiles -k /RUTA_databases/kofamscan_dbs/ko_list; done
cd
```

### 8. Exploración Metabólica con Rbims
Vamos a explorar el metabolismo utilizando [RbiMs](https://github.com/mirnavazquez/RbiMs). Primero, instalamos las bibliotecas necesarias en R:

```R
install.packages("devtools")
library(devtools)
install_github("mirnavazquez/RbiMs")
library(rbims)
library(tidyverse)
A continuación, leemos los resultados de KEGG y los mapeamos con la base de datos de KEGG:

pulque_mapp <- read_ko("08.Kofamscan/02.KO_results/") %>%
    mapping_ko()
Nos centraremos en las vías metabólicas relacionadas con la obtención de energía:

Overview <- c("Central Metabolism", "Carbon Fixation", "Nitrogen Metabolism", "Sulfur Metabolism", "Fermentation", "Methane Metabolism")
Energy_metabolisms_pulque <- pulque_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(rbims_pathway, Overview) 
Visualizamos los datos con un gráfico de burbujas:

plot_bubble(tibble_ko = Energy_metabolisms_pulque,
            x_axis = Bin_name, 
            y_axis = Pathway_cycle,
            analysis = "KEGG",
            calc = "Percentage",
            range_size = c(1, 10),
            y_labs = FALSE,
            x_labs = FALSE)  
Añadiremos metadatos, como la taxonomía:

Metadatos <- read_delim("11.GTDBTK/Metadatos.txt", delim = "\t")
Y generaremos un gráfico de burbujas con metadatos:

plot_bubble(tibble_ko = Energy_metabolisms_pulque,
            x_axis = Bin_name, 
            y_axis = Pathway_cycle,
            analysis = "KEGG",
            data_experiment = Metadatos,
            calc = "Percentage",
            color_character = Class,
            range_size = c(1, 10),
            y_labs = FALSE,
            x_labs = FALSE) 
Exploración de una Vía Específica
Podemos explorar una sola vía, como el “Secretion system,” y crear un mapa de calor para visualizar los genes relacionados con esta vía:

Secretion_system_pulque <- pulque_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(Cycle, "Secretion system")
Y, finalmente, generamos un mapa de calor:

plot_heatmap(tibble_ko = Secretion_system_pulque, 
             y_axis = Genes,
             analysis = "KEGG",
             calc = "Binary")
También podemos agregar metadatos para obtener una visión más completa:

plot_heatmap(tibble_ko = Secretion_system_pulque, 
             y_axis = Genes,
             data_experiment = Metadatos,
             order_x = Phylum,
             analysis = "KEGG",
             calc = "Binary")
plot_heatmap(tibble_ko = Secretion_system_pulque, 
             y_axis = Genes,
             data_experiment = Metadatos,
             order_y = Pathway_cycle,
             order_x = Phylum,
             analysis = "KEGG",
             calc = "Binary")
```

Estas visualizaciones te ayudarán a explorar y comprender mejor los datos metabólicos de tus bins.

****Detente a observar tus datos, 
- revisa la ayuda de cada programa,
- elige los parámetros que mejor se adapten a tus datos,
- analiza los resultados,
- haz varias pruebas,
- visita foros de ayuda,
- toma las mejores decisiones
- Disfruta el proceso****



   

