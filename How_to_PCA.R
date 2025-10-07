#How to PCA 
#Script to inspect and pre-process RNA-seq data counts 
#https://github.com/paulinapglz99

#I used the reference to build a PCA https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf

#Instalación de paquetes necesarios --- ---
# Solo debe ejecutarse la primera vez si no tiene los paquetes instalados.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Instalar los paquetes GEOquery y Biobase
BiocManager::install(c("GEOquery", "Biobase"), force = TRUE)

pacman::p_load("dplyr",         #Para manejar los datos
               "ggplot2",       #Para graficar
               "stringr",       #Para manejar los nombres de genes
               "gridExtra",     #Para hacer un panel 
               "GEOquery",      #Para descargar los datos
               "Biobase",       #Para descargar los datos
               "vroom", 
               "ggfortify",
               "tibble")         #Para descargar los datos

#Descarga del conjunto de datos de expresion de GEO (GSE119834) (Pediatric Glioblastoma)

expresion <- vroom::vroom("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE119834&format=file&file=GSE119834_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
dim(expresion)

#Descarga el objeto GEO para obtener la metadata
gse <- getGEO("GSE119834", GSEMatrix = TRUE, AnnotGPL = FALSE)

metadata <- pData(gse[[1]])

dim(metadata)

#Filtrar metadata

metadata <- metadata[metadata$geo_accession %in% colnames(expresion), ]
dim(metadata)

#En el caso de GSE119834, los controles son 'Neural Stem Cells (NSCs)'.

table(metadata$source_name_ch1)

#Prepare data

matrix <- expresion %>% 
  column_to_rownames("GeneID") %>% 
  as.matrix() %>% 
  t()  # transpose the matrix so that rows = samples and columns = variables

#Look at the first 10 rows and first 5 columns of the matrix
matrix[1:10, 1:10]

#Ahora la Magia!

#PCA prcomp

pca <- prcomp(matrix, retx = TRUE, center = TRUE, scale. = FALSE) 

#PCA to table

pca_df <- pca$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')

pca_df[1:10, 1:10]

#Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

variance_table <- data.frame(
  PC = 1:length(pca$sdev),
  Variance_Percentage = pca$sdev^2 / sum(pca$sdev^2) * 100,
  cumulative_percentage = cumsum(pca$sdev^2 / sum(pca$sdev^2) * 100))

head(variance_table)

#Elbow (Scree plot)

ggplot(variance_table, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = "navyblue",  position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'filtered data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  scale_x_discrete() +
  theme_minimal()

#Ahora hacemos la seleccion de PCs, la verdadera "reduccion de dimensionalidad"

#Table with the PCs explaining the 95% of data

variance_table_95 <- variance_table %>% 
  filter(cumulative_percentage >= 95)

#Plotting PCs explaining the 95% of data

variance_95 <- ggplot(variance_table_95, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = 'navyblue', position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'for the PCs explaining the 95% of data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  theme_minimal()

#Vis
variance_95

#Plot PC1 and PC2

PC1_PC2 <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = paste("PC1 vs PC2"),
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()

#Vis
PC1_PC2

#Plot PC1 and PC3

PC1_PC3 <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC3) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA de datos de Glioblastoma",
       subtitle = "PC1 vs PC3",
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)")) +
theme_minimal()

#Vis
PC1_PC3

#Ahora podemos colorearlo

PCA_color <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = metadata$`cell type:ch1`) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot",
       subtitle = paste("PC1 vs PC2"),
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()

#Vis
PCA_color

#Un plus a la visualizacion, colocar elipses para ver los 

PCA_elipse <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = metadata$`cell type:ch1`) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  #Agrega esta capa para dibujar las elipses
  stat_ellipse(geom = "polygon", aes(fill = metadata$`cell type:ch1`), alpha = 0.2, show.legend = FALSE) +
  labs(title = "PCA",
       subtitle = paste("PC1 vs PC2"),
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)"),
       colour = NULL) +
  theme_minimal()

#Vis
PCA_elipse

#Guardar plots

ggsave( filename = "PCA_elipse.pdf", 
       plot = PCA_elipse, 
      height = 20, 
      width = 20, 
      device = "pdf",  
      units = "cm"
        )

#FIN