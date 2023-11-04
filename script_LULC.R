

########MAPA ZGODNOSCI

setwd("E:/Projekty/projekty_22_23/LULC_stat/siatka/przyciete/ostateczne/koncowe")
library(raster)
library(terra)
library(kableExtra)
library(tmap)
library(sf)
library(ggplot2)
library(tibble)
library(tidyr)
library(plotly)

lc_names <- c("bdot", "lucas", "esa", "esri", "sent", "clc", "urban")

lc_list <- list()
lc_classes <- data.frame(matrix(nrow = 6, ncol = 7))

# Wczytanie danych pokrycia terenu
for (i in 1:length(lc_names)) {
  file_name <- paste0(lc_names[i], ".tif")
  rast <- rast(file_name)
  lc_list[[i]] <- rast
  lc_classes[, i] <- as.data.frame(table(as.data.frame(rast)))[, 2]
  assign(lc_names[i], rast)
}

names(lc_list) <- lc_names
names(lc_classes) <- lc_names

## odejmowanie i reklasyfikacja
lc_list <- lc_list[-1]
lc_names <- lc_names[-1]

for (i in seq_along(lc_list)){
  roznica <- bdot - lc_list[[i]]
  
  m <- c(-0.1, 0.1, 1,
         -6.1, -0.9, 0,
         0.9, 6.1, 0)
  matrix <- matrix(m, ncol = 3, byrow = TRUE)
  r_reclass <- classify(roznica, matrix, right = NA)
  assign(paste0(lc_names[i], "_diff"), r_reclass)
}

result <- clc_diff + esa_diff + esri_diff + lucas_diff + sent_diff + urban_diff

#writeRaster(result, "mapa_zgodnosci_2.tif", datatype = "INT1U")
plot(result, main = "OVaccuracy")















######## MACIERZE PRZEJSC ########


setwd("E:/Projekty/projekty_22_23/LULC_stat/dane")
poznan <- vect("poznan.gpkg")

lc_list <- c(bdot, lucas, esa, esri, sent, clc, urban)
lc_names <- c("bdot", "lucas", "esa", "esri", "sent", "clc", "urban")



#### POZNAN ####

#przycięcie warstw do granicy Poznania
for (i in 1:7) {
  cropped <- crop(lc_list[[i]], poznan, mask = TRUE)
  assign(paste0(lc_names[[i]], "_cropped"), cropped)
}

lc_cropped <- c(bdot_cropped, lucas_cropped, esa_cropped, esri_cropped, sent_cropped, clc_cropped, urban_cropped)

#utworzenie wektorów do obliczenia macierzy przejsc dla poznania
for (i in 1:7) {
  rast <- raster(lc_cropped[[i]])
  factor <- rasterToPoints(rast)
  factor <- as.data.frame(factor)
  factor <- as.factor(factor[, 3])
  assign(paste0(lc_names[[i]], "_crop_factor"), factor)
}

conf_matrix_poznan <- data.frame(matrix(nrow = 6, ncol = 2))
names(conf_matrix_poznan) <- c("accuracy_poznan", "kappa_poznan")
rownames(conf_matrix_poznan) <- lc_names[-1]
lc_factors <- list(lucas_crop_factor, esa_crop_factor, esri_crop_factor, sent_crop_factor, clc_crop_factor, urban_crop_factor)

PA_matrix_poznan <- data.frame(matrix(nrow = 6, ncol = 18))
rownames(PA_matrix_poznan) <- c("Obszary podmokłe", "Obszary wodne", "Zabudowa", "Lasy", "Obszary trawiaste i krzewiaste", "Pola uprawne")
colnames(PA_matrix_poznan)[c(1, 4, 7, 10, 13, 16)] <- paste0(lc_names[-1], "_PA")
colnames(PA_matrix_poznan)[c(2, 5, 8, 11, 14, 17)] <- paste0(lc_names[-1], "_UA")
colnames(PA_matrix_poznan)[c(3, 6, 9, 12, 15, 18)] <- paste0(lc_names[-1], "_F1")

kappa <- 0
### macierz przejść dla urban atlas
UA_matrix <- function(){
  urban_atlas_matrix <<- matrix(c(0, 0, 0, 0, 0, 0, 
                                  125,  67313,    129,   1020,    936,      2,
                                  279,   5636, 769551,  25311, 296850,  15208,
                                  2255,   4353,  31242, 399689,  66256,   9469,
                                  8857,   5045,  33300,  80616, 372198,  96743,
                                  3634,   1148,   2588,  13064,  66193, 237721),
                                nrow = 6, byrow = TRUE,
                                dimnames = list(c("1", "2", "3", "4", "5", "6"),
                                                c("1", "2", "3", "4", "5", "6")))
  
  n <- sum(urban_atlas_matrix)
  OA <<- sum(diag(urban_atlas_matrix)) / n
  
  rowsums <- apply(urban_atlas_matrix, 1, sum)
  p <- rowsums / n
  colsums <- apply(urban_atlas_matrix, 2, sum)
  q <- colsums / n
  expAccuracy <- sum(p*q)
  kappa <<- (OA - expAccuracy) / (1 - expAccuracy)
  
  conf_matrix_poznan[6, 1] <<- OA
  conf_matrix_poznan[6, 2] <<- kappa
  
}
UA_matrix()

#uzupełnienie macierzy przejsc
for (i in 1:6) {
  matrix <- table(lc_factors[[i]], bdot_crop_factor)
  if (i == 6) {
    matrix <- urban_atlas_matrix
  }
  
  ##accuracy
  n <- sum(matrix)
  diag <- diag(matrix)
  OA <- sum(diag) / n
  if (i != 6) {
    conf_matrix_poznan[i, 1] <- OA
  }
  
  ##kappa
  rowsums <- apply(matrix, 1, sum)
  p <- rowsums / n
  colsums <- apply(matrix, 2, sum)
  q <- colsums / n
  expAccuracy <- sum(p*q)
  kappa <- (OA - expAccuracy) / (1 - expAccuracy)
  if (i != 6) {
    conf_matrix_poznan[i, 2] <- kappa
  }
  
  ##producer accuracy
  PA <- diag / colsums
  col_idx <- (i - 1) * 3 + 1 
  PA_matrix_poznan[, col_idx] <- round(PA, 2)
  
  # user accuracy
  UA <- diag / rowsums
  col_idx <- (i - 1) * 3 + 2  
  PA_matrix_poznan[, col_idx] <- round(UA, 2)
  
  ## F1 score
  f1 = 2 * PA * UA / (PA + UA) 
  col_idx <- (i - 1) * 3 + 3  
  PA_matrix_poznan[, col_idx] <- round(f1, 2)
}




#### AGLOMERACJA ####
lc_list <- c(bdot, lucas, esa, esri, sent, clc, urban)
#utworzenie wektorów do obliczenia macierzy przejsc dla aglomeracji
for (i in 1:7) {
  rast <- raster(lc_list[[i]])
  factor <- rasterToPoints(rast)
  factor <- as.data.frame(factor)
  factor <- as.factor(factor[, 3])
  assign(paste0(lc_names[[i]], "_factor"), factor)
}

lc_factors <- list(lucas_factor, esa_factor, esri_factor, sent_factor, clc_factor, urban_factor)

PA_matrix_aglo <- data.frame(matrix(nrow = 6, ncol = 18))
rownames(PA_matrix_aglo) <- c("Obszary podmokłe", "Obszary wodne", "Zabudowa", "Lasy", "Obszary trawiaste i krzewiaste", "Pola uprawne")
colnames(PA_matrix_aglo)[c(1, 4, 7, 10, 13, 16)] <- paste0(lc_names[-1], "_PA")
colnames(PA_matrix_aglo)[c(2, 5, 8, 11, 14, 17)] <- paste0(lc_names[-1], "_UA")
colnames(PA_matrix_aglo)[c(3, 6, 9, 12, 15, 18)] <- paste0(lc_names[-1], "_F1")

conf_matrix_aglo <- data.frame(matrix(nrow = 6, ncol = 2))
names(conf_matrix_aglo) <- c("accuracy_aglo", "kappa_aglo")
rownames(conf_matrix_aglo) <- c("LUCAS", "ESA", "Esri", "Sentinel-2", "CLC", "Urban Atlas")

#uzupełnienie macierzy przejsc
for (i in 1:6) {
  matrix <- table(lc_factors[[i]], bdot_factor)
  
  ##accuracy
  n <- sum(matrix)
  diag <- diag(matrix)
  OA <- sum(diag) / n
  conf_matrix_aglo[i, 1] <- round(OA, 3)
  
  
  ##kappa
  rowsums <- apply(matrix, 1, sum)
  p <- rowsums / n
  colsums <- apply(matrix, 2, sum)
  q <- colsums / n
  expAccuracy <- sum(p*q)
  kappa <- (OA - expAccuracy) / (1 - expAccuracy)
  conf_matrix_aglo[i, 2] <- round(kappa, 3)
  
  
  ##producer accuracy
  PA <- diag / colsums
  col_idx <- (i - 1) * 3 + 1 
  PA_matrix_aglo[, col_idx] <- round(PA, 2)
  
  # user accuracy
  UA <- diag / rowsums
  col_idx <- (i - 1) * 3 + 2  
  PA_matrix_aglo[, col_idx] <- round(UA, 2)
  
  ## F1 score
  f1 = 2 * PA * UA / (PA + UA) 
  col_idx <- (i - 1) * 3 + 3  
  PA_matrix_aglo[, col_idx] <- round(f1, 2)
}

conf_matrix_aglo$accuracy_poznan <- round(conf_matrix_poznan$accuracy_poznan, 3)
conf_matrix_aglo$kappa_poznan <- round(conf_matrix_poznan$kappa_poznan, 3)


conf_matrix_aglo %>%
  kbl() %>%
  kable_minimal()

PA_matrix_aglo %>%
  kbl() %>%
  kable_minimal()




library(DT)
#datatable(PA_matrix_aglo, options = list(pageLength = 6,
#                             scrollX = TRUE, searching = FALSE))

datatable(conf_matrix_aglo, options = list(pageLength = 6,
                                         scrollX = TRUE, searching = FALSE,
                                         scrollY = FALSE,
                                         fontSize = "150%",
                                         dom = 't',
                                         #autoWidth = TRUE,
                                         columnDefs = list(list(width = '100px', fontSize = "150%", targets = "_all")))) %>% formatStyle(names(conf_matrix_aglo),
                              background = styleColorBar(range(0, 1), '#8cc0a5'),
                              color = "#000000",
                              fontSize = "160%",
                              backgroundSize = '98% 88%',
                              backgroundRepeat = 'no-repeat',
                              backgroundPosition = 'center') 
names(conf_matrix_aglo) <- c("AGLOMERACJA - dokładność", "AGLOMERACJA - kappa", "POZNAŃ - dokładność", "POZNAŃ - kappa")


names(PA_matrix_aglo) <- c("LUCAS - precyzja", "LUCAS - czułość", "LUCAS - F1", "ESA - precyzja", "ESA - czułość", "ESA - F1", "Esri - precyzja", "Esri - czułość", "Esri - F1", 
                           "Sentinel-2 - precyzja", "Sentinel-2 - czułość", "Sentinel-2 - F1", "CLC - precyzja", "CLC - czułość", "CLC - F1", "Urban Atlas - precyzja", "Urban Atlas - czułość", "Urban Atlas - F1")
names(PA_matrix_poznan) <- c("LUCAS - precyzja", "LUCAS - czułość", "LUCAS - F1", "ESA - precyzja", "ESA - czułość", "ESA - F1", "Esri - precyzja", "Esri - czułość", "Esri - F1", 
                           "Sentinel-2 - precyzja", "Sentinel-2 - czułość", "Sentinel-2 - F1", "CLC - precyzja", "CLC - czułość", "CLC - F1", "Urban Atlas - precyzja", "Urban Atlas - czułość", "Urban Atlas - F1")


#### STYL DLA WYKRESÓW - artykuł####
styl <- function(){
  style <<- theme(
    text = element_text(family = "Calibri"),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "#ffffff", color = NA),
    panel.background = element_rect(fill = "#ffffff"), 
    axis.title = element_text(size = 40,
                              color = "#222222"), 
    axis.title.x = element_text(vjust= -1),
    axis.title.y = element_text(vjust= 2.2),
    plot.title = element_text(size = 16,
                              color = "#222222",
                              vjust = 2,
                              hjust = 0.5), 
    legend.title = element_text(size = 42,  color = "#222222"),
    legend.text = element_text(size = 40, color = "#222222"),
    axis.text = element_text(size = 40, 
                             color = "#222222"),
    axis.text.x = element_text(angle = 0, vjust = 0.95, hjust=0.5),
    legend.spacing.y = unit(0.2, 'cm'),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin = margin(0, 0, 20, 15),
    panel.grid = element_line(color = "#222222"))
  
}
styl()

#### Procentowy udział poszczególnych zgodności (0-6) z mapy zgodności #####
### AGLOMERACJA

setwd("E:/Projekty/projekty_22_23/LULC_stat/siatka/przyciete/ostateczne/koncowe")
result <- rast("mapa_zgodnosci.tif")
result_df <- as.data.frame(result)

value_counts <- table(result_df)
value_counts <- as.data.frame(round(( value_counts / 20681815) * 100, 2))
colors <- c("#d73027", "#f46d43", "#fdae61", "#ffffbf", "#d9ef8b", "#66bd63", "#1a9850")

styl_inter <- function(){
  style <<- theme(
    text = element_text(family = "Calibri"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#222222"),
    plot.background = element_rect(fill = "#eff6f1", color = NA),
    panel.background = element_rect(fill = "#eff6f1"), 
    axis.title = element_text(size = 18,
                              color = "#222222"), 
    axis.title.y = element_text(vjust= 2.2),
    plot.title = element_text(size = 18,
                              color = "#222222",
                              vjust = 2,
                              hjust = 0.5), 
    legend.margin = margin(t = -15, r = 7, l = 7, b  = 7),
    legend.title = element_text(size = 18,  color = "#222222"),
    legend.text = element_text(size = 12, color = "#222222"),
    axis.text = element_text(size = 16, 
                             color = "#222222"),
    axis.text.x = element_text(angle = 0, vjust = 1.2, hjust=1),
    legend.spacing.y = unit(0.01, 'cm'),
    legend.position = "right",
    legend.key.size = unit(0.01, "cm"),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin = margin(10, 10, -40, 10))
}
styl_inter()

wykres <- ggplot(data = value_counts, aes(x = result_df, y = Freq, fill = result_df, text = paste('<span style = " font-weight:bold">Liczba zgodności:</span>',
                                                                                                  '<span>',result_df ,'</span>',
                                                                                                  '</br></br><span style = " font-weight:bold">Udział procentowy:</span>',
                                                                                                  '<span>',paste(Freq, "%") ,'</span>'))) +
          geom_bar(stat = "identity", color = "#222222") + 
          scale_fill_manual(values = colors) +
          xlab("Liczba źródeł pokrycia terenu zgodnych z BDOT10k") + 
          ylab("Procent obszaru aglomeracji [%]") + style + theme(legend.position = "none")

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne/poprawione")
svg("udzial_procentowy_aglo.svg", width = 9, height = 6.5, family = "Calibri")     
wykres + theme(axis.title = element_text(size = 24,
                                         color = "#222222"), 
               axis.text = element_text(size = 24, 
                                        color = "#222222"))
invisible(dev.off())
ggplotly(wykres, tooltip = "text")%>% config(displayModeBar = FALSE)


### Poznań

setwd("E:/Projekty/projekty_22_23/LULC_stat/dane")
poznan <- vect("poznan.gpkg")

result_poznan <- crop(result, poznan)
result_poznan <- mask(result_poznan, poznan)

result_poznan_df <- as.data.frame(result_poznan)

value_counts_poznan <- table(result_poznan_df)
value_counts_poznan <- as.data.frame(round(( value_counts_poznan / 2624464) * 100, 2))
colors <- c("#d73027", "#f46d43", "#fdae61", "#ffffbf", "#d9ef8b", "#66bd63", "#1a9850")

wykres <- ggplot(data = value_counts_poznan, aes(x = result_poznan_df, y = Freq, fill = result_poznan_df, text = paste('<span style = " font-weight:bold">Liczba zgodności:</span>',
                                                                                                                       '<span>',result_poznan_df ,'</span>',
                                                                                                                       '</br></br><span style = " font-weight:bold">Udział procentowy:</span>',
                                                                                                                       '<span>',paste(Freq, "%") ,'</span>'))) +
  geom_bar(stat = "identity", color = "#222222") + 
  scale_fill_manual(values = colors) +
  xlab("Liczba źródeł pokrycia terenu zgodnych z BDOT10k") + 
  ylab("Procent obszaru miasta Poznań [%]") + style + theme(legend.position = "none") +
  scale_y_continuous(limits = c(0, 50))

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne/poprawione")
svg("udzial_procentowy_poznan.svg", width = 9, height = 6.5, family = "Calibri")     
wykres + theme(axis.title = element_text(size = 24,
                                         color = "#222222"), 
               axis.text = element_text(size = 24, 
                                        color = "#222222"))
invisible(dev.off())
ggplotly(wykres, tooltip = "text")%>%
  config(displayModeBar = FALSE)



















########## UDZIAŁ PROCENTOWY KLAS #############

library(ggplot2)
library(dplyr)

data <- data.frame()


for (i in 1:7) {
  data_frame <- as.data.frame(table(as.data.frame(lc_cropped[[i]])))##lc_list
  data_frame$LC <- lc_names[i]
  data_frame$proc <- round((data_frame$Freq / sum(data_frame$Freq)) * 100, 2)
  
  data <- rbind(data, data_frame)
}

data <- data %>%
  mutate(Var1 = case_when(
    Var1 == 1 ~ "Obszary podmokłe",
    Var1 == 2 ~ "Obszary wodne",
    Var1 == 3 ~ "Zabudowa",
    Var1 == 4 ~ "Lasy",
    Var1 == 5 ~ "Obszary trawiaste i krzewiaste",
    Var1 == 6 ~ "Pola uprawne",
    TRUE ~ as.character(Var1) 
  ),
  LC = case_when(
    LC == "bdot" ~ "BDOT10k",
    LC == "clc" ~ "CLC",
    LC == "esa" ~ "ESA",
    LC == "esri" ~ "Esri",
    LC == "lucas" ~ "LUCAS",
    LC == "sent" ~ "Sentinel-2",
    LC == "urban" ~ "Urban Atlas",
    TRUE ~ as.character(LC) 
  ))



## styl dla wykresu
styl <- function(){
  style <<- theme(
    text = element_text(family = "Calibri"),
    panel.grid.major.x = element_blank(),
    plot.background = element_rect(fill = "#ffffff", color = NA),
    panel.background = element_rect(fill = "#ffffff"), 
    axis.title = element_text(size = 36,
                              color = "#222222"), 
    axis.title.y = element_text(vjust= 2.2),
    plot.title = element_text(size = 16,
                              color = "#222222",
                              vjust = 2,
                              hjust = 0.5), 
    legend.title = element_text(size = 40,  color = "#222222"),
    legend.margin = margin(t = 30, r = 7, l = 7, b  = 7),
    legend.text = element_text(size = 40, color = "#222222"),
    axis.text = element_text(size = 36, 
                             color = "#222222"),
    axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
    legend.spacing.y = unit(0.2, 'cm'),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin = margin(10, 10, 40, 10),
    panel.grid = element_line(color = "#222222"))
}
styl()

styl_inter <- function(){
  style <<- theme(
    text = element_text(family = "Calibri"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#222222"),
    plot.background = element_rect(fill = "#eff6f1", color = NA),
    panel.background = element_rect(fill = "#eff6f1"), 
    axis.title = element_text(size = 18,
                              color = "#222222"), 
    axis.title.y = element_text(vjust= 2.2),
    plot.title = element_text(size = 18,
                              color = "#222222",
                              vjust = 2,
                              hjust = 0.5), 
    legend.margin = margin(t = -15, r = 7, l = 7, b  = 7),
    legend.title = element_text(size = 18,  color = "#222222"),
    legend.text = element_text(size = 12, color = "#222222"),
    axis.text = element_text(size = 16, 
                             color = "#222222"),
    axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1),
    legend.spacing.y = unit(0.01, 'cm'),
    legend.position = "right",
    legend.key.size = unit(0.01, "cm"),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin = margin(10, 10, -40, 10))
}
styl_inter()

wykres1 <- ggplot(data, aes(x = LC, y = proc, fill = Var1, text = paste('<span style = " font-weight:bold"> Źródło:</span>',
                                                                              '<span>',LC ,'</span>',
                                                                              '</br></br>',
                                                                       '<span style = " font-weight:bold">Klasa:</span>',
                                                                             '<span>',Var1 ,'</span>','</br>',
                                                                       '<span style = " font-weight:bold">Udział procentowy:</span>',
                                                                              '<span>',paste(proc, "%") ,'</span>'))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Procentowy udział klas [%]", fill = "") +
  scale_fill_manual(values = c("Obszary podmokłe" = "#122b0e", "Obszary wodne" = "#61b5ff", "Zabudowa" = "#b83938",
                               "Lasy" = "#2d5e47", "Obszary trawiaste i krzewiaste" = "#9dd3b6", "Pola uprawne" = "#f4f430")) +
  style + guides(fill = guide_legend(byrow = TRUE)) #+ theme(legend.position = "none")#+ theme(plot.margin = margin(0, 0, 0, 70))

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne")
svg("udzial_procentowy_aglo.svg", width = 9, height = 6.5, family = "Calibri")     
wykres
invisible(dev.off())
ggplotly(wykres1, tooltip = c("text")) %>%
  config(displayModeBar = FALSE) #%>% layout(legend = list(x = 0, y = -0.3, orientation = "h"))


library(patchwork)
combined <- wykres + wykres1 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne/poprawione")
svg("udzial_proc.svg", width = 24, height = 11.25, family = "Calibri")     
combined
invisible(dev.off())






 

####### wykres błędów klasyfikacji dla aglomeracji ########
error_factors <- list(lucas_factor, esa_factor, esri_factor, sent_factor, clc_factor, urban_factor)
error_aglo <- data.frame(matrix(nrow = 6, ncol = 6))
colnames(error_aglo) <- c("LUCAS", "ESA", "Esri", "Sentinel-2", "CLC", "Urban Atlas")


for (i in 1:6){
  matrix <- table(error_factors[[i]], bdot_factor)
  error_aglo[, i] <- round((1 - (diag(matrix) / apply(matrix, 1, sum))) * 100, 2)
}

#error_aglo$error <- round(apply(error_aglo, 1, mean), 2)
error_aglo$LC <- c("Obszary podmokłe", "Obszary wodne", "Zabudowa", "Lasy", "Obszary trawiaste i krzewiaste", "Pola uprawne")


## styl dla wykresu
styl <- function(){
  style <<- theme(
    text = element_text(family = "Calibri"),
    panel.grid.major.y = element_line(color = "#222222"),
    panel.grid.minor.y = element_line(color = "#222222"),
    plot.background = element_rect(fill = "#ffffff", color = NA),
    panel.background = element_rect(fill = "#ffffff"), 
    axis.title = element_text(size = 32,
                              color = "#222222"), 
    axis.title.y = element_text(vjust= 2.2, hjust = 0.6),
    plot.title = element_text(size = 16,
                              color = "#222222",
                              vjust = 2,
                              hjust = 0.5), 
    legend.margin = margin(t = 15, r = 7, l = 7, b  = 7),
    legend.title = element_text(size = 16,  color = "#222222"),
    legend.text = element_text(size = 36, color = "#222222"),
    axis.text = element_text(size = 36, 
                             color = "#222222"),
    axis.text.x = element_text(size = 36, angle = 45, vjust = 1, hjust=1),
    legend.spacing.y = unit(0.2, 'cm'),
    legend.key = element_blank(),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin = margin(10, 10, 0, 10))
}
styl()
styl_inter <- function(){
  style <<- theme(
    text = element_text(family = "Calibri"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#222222"),
    panel.grid.minor.y = element_line(color = "#222222"),
    plot.background = element_rect(fill = "#eff6f1", color = NA),
    panel.background = element_rect(fill = "#eff6f1"), 
    axis.title = element_text(size = 18,
                              color = "#222222"), 
    axis.title.y = element_text(vjust= 2.2),
    plot.title = element_text(size = 18,
                              color = "#222222",
                              vjust = 2,
                              hjust = 0.5), 
    legend.margin = margin(t = -15, r = 7, l = 7, b  = 7),
    legend.title = element_text(size = 18,  color = "#222222"),
    legend.text = element_text(size = 12, color = "#222222"),
    axis.text = element_text(size = 16, 
                             color = "#222222"),
    axis.text.x = element_text(angle = 45, vjust = 1.2, hjust=1),
    legend.spacing.y = unit(0.01, 'cm'),
    legend.position = "bottom",
    legend.key.size = unit(0.01, "cm"),
    legend.background = element_rect(fill = "transparent", color = "transparent"),
    plot.margin = margin(10, 10, -40, 10))
}
styl_inter()

ggplot(error_aglo, aes(x = reorder(LC, - error), y = error, fill = LC)) + geom_bar(stat = "identity") + style + 
  scale_fill_manual(values = c("Obszary podmokłe" = "#122b0e", "Obszary wodne" = "#61b5ff", "Zabudowa" = "#b83938",
                               "Lasy" = "#2d5e47", "Obszary trawiaste i krzewiaste" = "#9dd3b6", "Pola uprawne" = "#f4f430")) + 
  ylab("% powierzchni poszczególnych klas") + 
  xlab("") + style + theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 10, 10)) + 
  scale_y_continuous(limits = c(0, 80))
#Średni procent obszarów niepoprawnie zaklasyfikowanych względem BDOT w obu powiatach






### błędy klasyfikacji dla Poznania
error_factors <- list(lucas_crop_factor, esa_crop_factor, esri_crop_factor, sent_crop_factor, clc_crop_factor)
error_poznan <- data.frame(matrix(nrow = 6, ncol = 6))
colnames(error_poznan) <- c("LUCAS", "ESA", "Esri", "Sentinel-2", "CLC", "Urban Atlas")

for (i in 1:5){
  matrix <- table(error_factors[[i]], bdot_crop_factor)
  error_poznan[, i] <- round((1 - (diag(matrix) / apply(matrix, 1, sum))) * 100, 2)
}

UA_matrix <- function(){
  urban_atlas_matrix <- matrix(c(0, 0, 0, 0, 0, 0, 
                                 125,  67313,    129,   1020,    936,      2,
                                 279,   5636, 769551,  25311, 296850,  15208,
                                 2255,   4353,  31242, 399689,  66256,   9469,
                                 8857,   5045,  33300,  80616, 372198,  96743,
                                 3634,   1148,   2588,  13064,  66193, 237721),
                               nrow = 6, byrow = TRUE,
                               dimnames = list(c("1", "2", "3", "4", "5", "6"),
                                               c("1", "2", "3", "4", "5", "6")))
  
  error_poznan[, 6] <<- round((1 - (diag(urban_atlas_matrix) / apply(urban_atlas_matrix, 1, sum))) * 100, 2)
  error_poznan[1, 6] <<- 0
}
UA_matrix()

#error_poznan$error <- apply(error_poznan, 1, function(x) mean(x, na.rm = TRUE))
error_poznan$LC <- c("Obszary podmokłe", "Obszary wodne", "Zabudowa", "Lasy", "Obszary trawiaste i krzewiaste", "Pola uprawne")

ggplot(error_poznan, aes(x = reorder(LC, - error), y = error, fill = LC)) + geom_bar(stat = "identity") + style + 
  scale_fill_manual(values = c("Obszary podmokłe" = "#122b0e", "Obszary wodne" = "#61b5ff", "Zabudowa" = "#b83938",
                               "Lasy" = "#2d5e47", "Obszary trawiaste i krzewiaste" = "#9dd3b6", "Pola uprawne" = "#f4f430")) + 
  ylab("% powierzchni poszczególnych klas") + 
  xlab("") + style + theme(legend.position = "none") + theme(plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(limits = c(0, 80))






# wykres pogrupowany dla poznania
error_poznan_long <- error_poznan %>%
  pivot_longer(cols = -LC, names_to = "LandCover", values_to = "error")

error_poznan_long <- error_poznan_long %>%
  mutate(LandCover = case_when(
    LandCover == "Urban Atlas" ~ "Urban Atlas*",
    TRUE ~ as.character(LandCover) 
  ))
  
plot <- function() {
  wykres <<- ggplot(error_poznan_long, aes(fill=LC, y=error, x=LandCover,  text = paste('<span style = " font-weight:bold"> Źródło:</span>',
                                                                            '<span>',LC ,'</span>',
                                                                            '</br></br>',
                                                                            '<span style = " font-weight:bold">Klasa:</span>',
                                                                            '<span>',LandCover ,'</span>','</br>',
                                                                            '<span style = " font-weight:bold">Udział procentowy:</span>',
                                                                            '<span>',paste(error, "%") ,'</span>'))) + 
    geom_bar(position=position_dodge(width = 0.7), stat="identity", color = "#555555") + 
    scale_fill_manual(values = c("Obszary podmokłe" = "#122b0e", "Obszary wodne" = "#61b5ff", "Zabudowa" = "#b83938",
                                 "Lasy" = "#2d5e47", "Obszary trawiaste i krzewiaste" = "#9dd3b6", "Pola uprawne" = "#f4f430")) +
    style + labs(x = "", y = "Procent powierzchni klas [%]", fill = "", caption = "*brak obszarów podmokłych\n na terenie Poznania") + 
    guides(fill = guide_legend(byrow = TRUE)) + theme(plot.caption = element_text(vjust = 30, hjust = 0.05, size = 20, color = "#222222"),
                                                      plot.caption.position = "plot")
    
}
plot()

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne")
svg("error_poznan.svg", width = 9, height = 6.5, family = "Calibri")     
wykres
invisible(dev.off())


ggplotly(wykres, tooltip = c("text")) %>%
  config(displayModeBar = FALSE) %>% 
  layout(annotations = list(x = 1.4, y = 0.05, text = "*brak obszarów podmokłych\n na terenie Poznania",
                     xref='paper', yref='paper', showarrow = F,
                     font = list(size = 14, color = "#222222"))) #%>% layout(legend = list(x = 0, y = -0.3, orientation = "h"))




# wykres pogrupowany dla aglomeracji
error_aglo_long <- error_aglo %>%
  pivot_longer(cols = -LC, names_to = "LandCover", values_to = "error")

plot <- function() {
  wykres1 <<-ggplot(error_aglo_long, aes(fill=LC, y=error, x=LandCover, text = paste('<span style = " font-weight:bold"> Źródło:</span>',
                                                                                    '<span>',LC ,'</span>',
                                                                                   '</br></br>',
                                                                                    '<span style = " font-weight:bold">Klasa:</span>',
                                                                                    '<span>',LandCover ,'</span>','</br>',
                                                                                    '<span style = " font-weight:bold">Udział procentowy:</span>',
                                                                                    '<span>',paste(error, "%") ,'</span>'))) + 
    geom_bar(position=position_dodge(width = 0.7), stat="identity", color = "#555555") + 
    scale_fill_manual(values = c("Obszary podmokłe" = "#122b0e", "Obszary wodne" = "#61b5ff", "Zabudowa" = "#b83938",
                                 "Lasy" = "#2d5e47", "Obszary trawiaste i krzewiaste" = "#9dd3b6", "Pola uprawne" = "#f4f430")) +
    style + labs(x = "", y = "Procent powierzchni klas [%]", fill = "") + 
    guides(fill = guide_legend(byrow = TRUE)) + scale_y_continuous(limits = c(0, 100)) + theme(legend.position = "none")#+ theme(plot.margin = margin(0, 0, 0, 70))
}
plot()

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne")
svg("error_aglo.svg", width = 9, height = 6.5, family = "Calibri")     
wykres
invisible(dev.off())


ggplotly(wykres1, tooltip = c("text")) %>%
  config(displayModeBar = FALSE)



#library()
combined <- wykres + wykres1 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

setwd("E:/Projekty/projekty_22_23/LULC_stat/wykresy_artykul/poprawne/poprawione")
svg("error.svg", width = 24, height = 11.25, family = "Calibri")     
combined
invisible(dev.off())

ggplotly(combined)

#subplot(style(a, showlegend = F), b) %>% 
 # layout(title = 'Procentowy udział klas po reklasyfikacji dla wszystkich źródeł danych \n
 #        Aglomeracja                                                                                 Poznań', legend = list(orientation = 'h', x = 0, y = -0.3))




