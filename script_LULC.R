setwd("E:/Projekty/projekty_22_23/LULC_stat/siatka/przyciete/ostateczne/koncowe")

##Producer and user accuracy:
  
# Producer accuracy
# uncomment lines to run them
# PA <- diag / colsums
# # User accuracy
# UA <- diag / rowsums
# (outAcc <- data.frame(producerAccuracy = PA, userAccuracy = UA))
# 
# 
# 
# (455/34500000) *100
# 
# ######
# 
# (a <- table(values(bdot_ut)))
# (b <- table(values(lucas)))
# 
# a<- as.data.frame(a)
# a <- as.factor(a)
# 
# b <- as.data.frame(clc)
# b<- as.factor(b)
# sum(a[1:6]) - sum(b[1:6])
# 
# lucas_2 <- rast("lucas.tif")
# bdot_2 <- rast("bdot.tif")
# plot(lucas)











######## MAPA ZGODNOSCI

setwd("E:/Projekty/projekty_22_23/LULC_stat/siatka/przyciete/ostateczne/koncowe")
library(raster)
library(terra)
library(tmap)
library(sf)
library(ggplot2)
library(tibble)


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


setwd("E:/Projekty/projekty_22_23/LULC_stat/mapy/mapa_zgodnosci/poznan_przyciete")
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
names(conf_matrix_poznan) <- c("accuracy", "kappa")
rownames(conf_matrix_poznan) <- lc_names[-1]
lc_factors <- list(lucas_crop_factor, esa_crop_factor, esri_crop_factor, sent_crop_factor, clc_crop_factor)

#uzupełnienie macierzy przejsc
for (i in 1:5) {
  matrix <- table(lc_factors[[i]], bdot_crop_factor)
  
  ##accuracy
  n <- sum(matrix)
  diag <- diag(matrix)
  OA <- sum(diag) / n
  conf_matrix_poznan[i, 1] <- OA
  
  
  ##kappa
  rowsums <- apply(matrix, 1, sum)
  p <- rowsums / n
  colsums <- apply(matrix, 2, sum)
  q <- colsums / n
  expAccuracy <- sum(p*q)
  kappa <- (OA - expAccuracy) / (1 - expAccuracy)
  conf_matrix_poznan[i, 2] <- kappa
}

### macierz przejść dla urban atlas
UA_matrix <- function(){
  urban_atlas_matrix <- matrix(c(0, 0, 0, 0, 0, 0, 
                              125, 67313, 129, 1020, 935, 3,
                              279, 5636, 769551, 25311, 277147, 34911,
                              2255, 4353, 31242, 399689, 62059, 13666,
                              8857, 5045, 33300, 80616, 297823, 171118,
                              3634, 1148, 2588, 13064, 61369, 242545),
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

conf_matrix_aglo <- data.frame(matrix(nrow = 6, ncol = 2))
names(conf_matrix_aglo) <- c("accuracy", "kappa")
rownames(conf_matrix_aglo) <- lc_names[-1]

#uzupełnienie acierzy przejsc
for (i in 1:6) {
  matrix <- table(lc_factors[[i]], bdot_factor)
  
  ##accuracy
  n <- sum(matrix)
  diag <- diag(matrix)
  OA <- sum(diag) / n
  conf_matrix_aglo[i, 1] <- OA
  
  
  ##kappa
  rowsums <- apply(matrix, 1, sum)
  p <- rowsums / n
  colsums <- apply(matrix, 2, sum)
  q <- colsums / n
  expAccuracy <- sum(p*q)
  kappa <- (OA - expAccuracy) / (1 - expAccuracy)
  conf_matrix_aglo[i, 2] <- kappa
}










#### STYL DLA WYKRESÓW ####

style <- theme(
  text = element_text(family = "Calibri"),
  panel.grid.major.x = element_blank(),
  plot.background = element_rect(fill = "#4b6865"),
  panel.background = element_rect(fill = "#4b6865"), 
  axis.title = element_text(size = 14,
                            color = "#ffffff"), 
  plot.title = element_text(size = 16,
                            color = "#ffffff",
                            vjust = 2,
                            hjust = 0.5), 
  legend.background = element_rect(color = "#ffffff", 
                                   fill = "#4b6865"),
  legend.title = element_text(size = 16,  color = "#ffffff"),
  legend.text = element_text(size = 14, color = "#ffffff"),
  axis.text = element_text(size = 14, 
                           color = "#ffffff"),
  axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1))










#### Procentowy udział poszczególnych zgodności (0-6) z mapy zgodności #####
### AGLOMERACJA

setwd("E:/Projekty/projekty_22_23/LULC_stat/siatka/przyciete/ostateczne/koncowe")
result <- rast("mapa_zgodnosci.tif")
result_df <- as.data.frame(result)

value_counts <- table(result_df)
value_counts <- as.data.frame(round(( value_counts / 20681815) * 100, 2))

colors <- c("#d73027", "#f46d43", "#fdae61", "#ffffbf", "#d9ef8b", "#66bd63", "#1a9850")
ggplot(data = value_counts, aes(x = result_df, y = Freq, fill = result_df)) +
  geom_bar(stat = "identity", color = "#222222") + 
  scale_fill_manual(values = colors) +
  xlab("Liczba źródeł pokrycia terenu zgodnych z BDOT") + 
  ylab("Procent udziału w całym obszarze") + style + theme(legend.position = "none") 



### Poznań

setwd("E:/Projekty/projekty_22_23/LULC_stat/mapy/mapa_zgodnosci/poznan_przyciete")
poznan <- vect("poznan.gpkg")

result_poznan <- crop(result, poznan)
result_poznan <- mask(result_poznan, poznan)

result_poznan_df <- as.data.frame(result_poznan)

value_counts <- table(result_poznan_df)
value_counts <- as.data.frame(round(( value_counts / 2624464) * 100, 2))

colors <- c("#d73027", "#f46d43", "#fdae61", "#ffffbf", "#d9ef8b", "#66bd63", "#1a9850")
ggplot(data = value_counts, aes(x = result_poznan_df, y = Freq, fill = result_poznan_df)) +
  geom_bar(stat = "identity", color = "#222222") + 
  scale_fill_manual(values = colors) +
  xlab("Liczba źródeł pokrycia terenu zgodnych z BDOT") + 
  ylab("Procent udziału w całym obszarze") + style + theme(legend.position = "none") 





















########## UDZIAŁ PROCENTOWY KLAS #############

library(ggplot2)
library(dplyr)

data <- data.frame()

for (i in 1:7) {
  data_frame <- as.data.frame(table(as.data.frame(lc_list[[i]])))
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
  ))


wykres <- ggplot(data, aes(x = LC, y = proc, fill = Var1)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "", y = "Procentowy udział", fill = "Klasa pokrycia") +
  scale_fill_manual(values = c("Obszary podmokłe" = "#122b0e", "Obszary wodne" = "#61b5ff", "Zabudowa" = "#b83938",
                               "Lasy" = "#2d5e47", "Obszary trawiaste i krzewiaste" = "#9dd3b6", "Pola uprawne" = "#f4f430")) +
  style


print(wykres)

