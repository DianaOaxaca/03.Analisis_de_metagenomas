# Aqui describo todos los pasos del flujo de analisis

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

Listo el set de datos está en data
