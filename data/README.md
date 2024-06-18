# Aqui describo los pasos para obtener los datos

## 01. Obtener el filereport del ENA para descargar los datos
https://www.ebi.ac.uk/ena/browser/view/PRJNA603591

```bash
data/filereport_ena_PRJNA603591.txt
```

## 02. Descargar los datos
Dentro de data ejecutar el script download_pulque.sh que toma el contenido de filereport_ena_PRJNA603591.txt y hace un ftp para la descarga, 
ademas verifica la integridad del archivo con md5sum y agrega el nombre de la muestra correspondiente

```bash
cd data/
bash download_pulque.sh
```
Descomprimirlos

```bash
pigz -d *.gz
```

## 03. Tomar un subseq de 10 M de lecturas de los datos para correr la practica rapido

```bash
cd ..

bash src/00.subseq.sh
```

Listo el set de datos está en `data` el contenido es el siguiente:


```bash
tree data/
data/
├── Pulque-AM_SRR10997050_1.fastq
├── Pulque-AM_SRR10997050_1_10M.fastq
├── Pulque-AM_SRR10997050_2.fastq
├── Pulque-AM_SRR10997050_2_10M.fastq
├── Pulque-PQ_SRR10997049_1.fastq
├── Pulque-PQ_SRR10997049_1_10M.fastq
├── Pulque-PQ_SRR10997049_2.fastq
├── Pulque-PQ_SRR10997049_2_10M.fastq
├── Pulque-T0_SRR10997048_1.fastq
├── Pulque-T0_SRR10997048_1_10M.fastq
├── Pulque-T0_SRR10997048_2.fastq
├── Pulque-T0_SRR10997048_2_10M.fastq
├── Pulque-T3_SRR10997047_1.fastq
├── Pulque-T3_SRR10997047_1_10M.fastq
├── Pulque-T3_SRR10997047_2.fastq
├── Pulque-T3_SRR10997047_2_10M.fastq
├── Pulque-T6_SRR10997046_1.fastq
├── Pulque-T6_SRR10997046_1_10M.fastq
├── Pulque-T6_SRR10997046_2.fastq
├── Pulque-T6_SRR10997046_2_10M.fastq
├── download_pulque.sh
└── filereport_ena_PRJNA603591.txt

0 directories, 22 files
```
