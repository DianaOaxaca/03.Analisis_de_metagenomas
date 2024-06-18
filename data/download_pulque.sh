#!/bin/bash

# Archivo de entrada
input_file="filereport_ena_PRJNA603591.txt"

# Filtrar las líneas que contienen "Pulque" y procesar cada una
grep Pulque-T0 "$input_file" | cut -f4,5,6 | while IFS=$'\t' read -r md5 ftp sample; do
  # Separar los md5sum y los ftp en arrays
  IFS=';' read -r -a md5_array <<< "$md5"
  IFS=';' read -r -a ftp_array <<< "$ftp"

  # Comprobar que ambos arrays tienen la misma longitud
  if [ ${#md5_array[@]} -eq ${#ftp_array[@]} ]; then
    for i in "${!md5_array[@]}"; do
      # Descargar el archivo
      wget -O "${sample}_$(basename ${ftp_array[i]})" "${ftp_array[i]}"

      # Comprobar el md5sum
      echo "${md5_array[i]}  ${sample}_$(basename ${ftp_array[i]})" | md5sum -c -
    done
  else
    echo "Error: Número de md5sum y ftp no coinciden para $sample"
  fi
done
