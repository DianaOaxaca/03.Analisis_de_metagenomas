###### Metabolismo
#Instalar paquetes y cargar librerias
#install.packages("devtools")
library(devtools)
#install_github("mirnavazquez/RbiMs", force=TRUE)
library(rbims)
library(tidyverse)
#leer los resultados de KEEG y mapearlos con el resto de la base de datos      ####
#de KEEG
pulque_mapp<-read_ko("kegg/") %>%
  mapping_ko()
# vamos a enfocarnos en los metabolismos encargados de la obtención de energía.
Overview<-c("Central Metabolism", "Carbon Fixation",
            "Nitrogen Metabolism", "Sulfur Metabolism", "Fermentation",
            "Methane Metabolism")
Energy_metabolisms_pulque<-pulque_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(rbims_pathway, Overview)

# visualizar los datos.
plot_bubble(tibble_ko = Energy_metabolisms_pulque,
            x_axis = Bin_name,
            y_axis = Pathway_cycle,
            analysis="KEGG",
            calc="Percentage",
            range_size = c(1,10),
            y_labs=FALSE,
            x_labs=FALSE)
#metadatos
Metadatos<-read_delim("Metadatos.txt", delim="\t")

plot_bubble(tibble_ko = Energy_metabolisms_pulque,
            x_axis = Bin_name,
            y_axis = Pathway_cycle,
            analysis="KEGG",
            data_experiment = Metadatos,
            calc="Percentage",
            color_character = Class,
            range_size = c(1,10),
            y_labs=FALSE,
            x_labs=FALSE)

#una sola via
Secretion_system_pulque<-pulque_mapp %>%
  drop_na(Cycle) %>%
  get_subset_pathway(Cycle, "Secretion system")

#heatmap
plot_heatmap(tibble_ko=Secretion_system_pulque,
             y_axis=Genes,
             analysis = "KEGG",
             calc="Binary")
#Agregar metadatos
plot_heatmap(tibble_ko=Secretion_system_pulque,
             y_axis=Genes,
             data_experiment = Metadatos,
             order_x = Phylum,
             analysis = "KEGG",
             calc="Binary")


plot_heatmap(tibble_ko=Secretion_system_pulque,
             y_axis=Genes,
             data_experiment = Metadatos,
             order_y = Pathway_cycle,
             order_x = Phylum,
             analysis = "KEGG",
             calc="Binary")

