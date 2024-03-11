



## This R code corresponds to the section "Application to real data" of the paper 
## Small area estimation of labour force indicators under unit-level multinomial mixed models.
## This is the script to obtain maps of the estimations and its mse's for the Multinomial Mixed Model.


## AUTHORS: Bugallo M., Esteban M.D., Hozba T., Morales D., Perez A.


library(maptools) # Maptools in now archived. If you don't have any version, try install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
library(RColorBrewer)


GroupClassification <- function(data, datacompare, intervals){
  n <- length(data)
  group <- matrix(0,nrow=n) 
  ninterv <- length(intervals)
  for (i in 1:n){
    for (j in 1:ninterv)
      if (datacompare[i]<intervals[j]){
        group[i] <- intervals[j]
        break
      }
  }
  result <- split(data,group)
  return(result)
}

PrintSpainMap <- function(pathmap, datos, colors, titlemap, textlegend){
  m <- matrix(c(1,1,1,2),2,2)
  layout(m, widths=c(1.5, 1), heights=c(1.5, 1), respect=F)
  xName <- readShapePoly(pathmap, IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))     
  xName$datos <- NA
  for (i in 1:length(colors))
    xName$datos[datos[[i]]] <- colors[i]
  
  xSC <- xName[xName$ESP_PROV_I < 35 | xName$ESP_PROV_I >38 | xName$ESP_PROV_I==36 | xName$ESP_PROV_I ==37,]
  plot(xSC,  xlab="",  col=xSC$datos, axes=F)
  title(titlemap, line=-1.2, cex.main=1.7)
  legend( "topright", textlegend, col = 1, pt.bg = colors, pch=21, bty = "n", cex=1.5, pt.cex = 2.5) 
  xC <- xName[xName$ESP_PROV_I == 35 | xName$ESP_PROV_I == 38,]
  plot(xC,  xlab="",  col=xC$datos)
  box()
}

pathmap <- "spainmap/esp_prov.shp"


## UNEMPLOYMENT RATES

datos <- read.table(file="UNEMP.RATE.txt", header=TRUE, sep=';', dec='.')
dom <- datos$prov

## Young men: sex1 and age4-1
estML <- datos$UNEMPRate_11
intervals_prop <- c(0,10,15,20, 25, Inf)  # Intervals  
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7,8,9)] ) # Colors
legend_prop <- expression("< 10 %", "10-15 %", "15-20 %", "20-25 %", "25-30 %") # Legend 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young men: age4-1", legend_prop)

## Young women: sex6 and age4-1
estML <- datos$UNEMPRate_16
intervals_prop <- c(0,15,20, 25, 30,40, Inf) # Intervals  
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8,9,10)] ) # Colors
legend_prop <- expression("10-15 %", "15-20 %", "20-25 %", "25-30 %", "30-40 %") # Legend 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young women: age4-1", legend_prop)

## Middle age men: sex1 and age4-2
estML <- datos$UNEMPRate_21
intervals_prop <- c(0,10,15,20, Inf)  # Intervals  
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7,8)] ) # Colors
legend_prop <- expression("< 10 %", "10-15 %", "15-20 %", "20-25 %") # Legend 
result <- GroupClassification(dom,estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop,  "Middle age men: age4-2", legend_prop)

## Middle age women: sex6 and age4-2
estML <- datos$UNEMPRate_26
intervals_prop <- c(0,10, 15,20, 25, 30, Inf) # Intervals  
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7,8,9)] ) # Colors
legend_prop <- expression("< 10 %", "10-15 %", "15-20 %", "20-25 %", "25-30 %") # Legend 
result <- GroupClassification(dom,estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Middle age women age4-2", legend_prop)

## Adult men: sex1 and age4-3
estML <- datos$UNEMPRate_31
intervals_prop <- c(0,10,15,20, Inf)  # Intervals  
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7)] ) # Colors
legend_prop <- expression("< 10 %", "10-15 %", "15-20 %") # Legend 
result <- GroupClassification(dom,estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult men: age4-3", legend_prop)

## Adult women: sex6 and age4-3
estML <- datos$UNEMPRate_36
intervals_prop <- c(0,10, 15,20, 25, 30, Inf) # Intervals 
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7,8,9)] ) # Colors
legend_prop <- expression("< 10 %", "10-15 %", "15-20 %", "20-25 %", "25-30 %") # Legend 
result <- GroupClassification(dom,estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult women age4-3", legend_prop)


## INACTIVE PROPORTIONS 

datos <- read.table(file="INACTIVE.txt", header=TRUE, sep=';', dec='.')
dom <- datos$prov

## Young men: sex1 and age4-1
estML <- datos$inac_11
intervals_prop <- c(0.10, 0.20, 0.30,  Inf) # Intervals 
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7)] ) # Colors
legend_prop <- expression("10-20 %", "20-30 %", "30-40 %") # Legend 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young men: age4-1", legend_prop)

### Young women: sex6 and age4-1
estML <- datos$inac_16
intervals_prop <- c(0.20, 0.30,  Inf) 
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7)] ) 
legend_prop <- expression("20-30 %", "30-40 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young women: age4-1", legend_prop)

## Middle age men: sex1 and age4-2
estML <- datos$inac_21
intervals_prop <- c(0.10, 0.20,  Inf) 
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6)] ) 
legend_prop <- expression("10-20 %", "20-30 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Middle age men: age4-2", legend_prop)

## Middle age women: sex6 and age4-2
estML <- datos$inac_26
intervals_prop <- c(0.10, 0.20, 0.30,  Inf) 
colorsprop <- c("white", brewer.pal(10,"RdYlBu")[c(6,7)] ) 
legend_prop <- expression("10-20 %", "20-30 %", "30-40 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Middle age women: age4-2", legend_prop)

## Adult men: sex1 and age4-3
estML <- datos$inac_31
intervals_prop <- c(0.30, 0.40, 0.50,  Inf) 
colorsprop <- c( brewer.pal(10,"RdYlBu")[c(7, 8, 9)] ) 
legend_prop <- expression("30-40 %", "40-50 %", "50-60 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult men: age4-3", legend_prop)

## Adult women: sex6 and age4-3
estML <- datos$inac_36
intervals_prop <- c(0.40, 0.50,  0.60,  Inf) 
colorsprop <- c( brewer.pal(10,"RdYlBu")[c(8, 9, 10)] ) 
legend_prop <- expression("40-50 %", "50-60 %", "60-70 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult women: age4-3", legend_prop)


## EMPLOYMENT PROPORTIONS

datos <- read.table(file="EMPLOYED.txt", header=TRUE, sep=';', dec='.')
dom <- datos$prov

## Young men: sex1 and age4-1
estML <- datos$ocup_11
intervals_prop <- c(0.50, 0.60, 0.70,  Inf)
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(7,8,10)] )
legend_prop <- expression("50-60 %", "60-70 %", "70-80 %")
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young men: age4-1", legend_prop)

## Middle age men: sex1 and age4-2
estML <- datos$ocup_21
intervals_prop <- c(0.60, 0.70,  Inf)
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(8,10)] )
legend_prop <- expression("60-70 %", "70-80 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Mid adult men: age4-2", legend_prop)

## Adult men: sex1 and age4-3
estML <- datos$ocup_31
intervals_prop <- c(0.3, 0.4, 0.50,  Inf) 
colorsprop <- c('white',brewer.pal(10,"RdYlBu")[c(6,7)] )
legend_prop <- expression("< 40 %", "40-50 %", "50-60 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult men: age4-3", legend_prop)

## Young women: sex6 and age4-1
estML <- datos$ocup_16
intervals_prop <-Inf 
colorsprop <- 'white'
legend_prop <- expression("< 40 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young women: age4-1", legend_prop)

## Middle age women: sex6 and age4-2
estML <- datos$ocup_26
intervals_prop <- c(0.4, 0.5, 0.60, 0.70,  Inf)  
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8,10)] )
legend_prop <- expression( "40-50 %", "50-60 %", "60-70 %", "70-80 %")
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Mid adult women: age4-2", legend_prop)

## Adult women: sex6 and age4-3
estML <- datos$ocup_36
intervals_prop <- c(0.25, 0.4, 0.50,  Inf) 
colorsprop <- c('white',brewer.pal(10,"RdYlBu")[c(6,7)] )
legend_prop <- expression("< 40 %", "40-50 %", "50-60 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult women: age4-3", legend_prop)


## UNEMPLOYMENT PROPORTIONS 

datos <- read.table(file="UNEMPLOYED.txt", header=TRUE, sep=';', dec='.')
dom <- datos$prov

## Young men: sex1 and age4-1
estML <- datos$parado_11
intervals_prop <- c(0.05, 0.10, 0.15,  Inf) 
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8)] )
legend_prop <- expression("5-10 %", "10-15 %", '15-20') 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young men: age4-1", legend_prop)

## Middle age men: sex1 and age4-2
estML <- datos$parado_21
intervals_prop <- c(0.05, 0.10, 0.15,  Inf) 
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8)] )
legend_prop <- expression("5-10 %", "10-15 %", '15-20') 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Mid adult men: age4-2", legend_prop)

## Adult men: sex1 and age4-3
estML <- datos$parado_31
intervals_prop <- c(0, 0.05,  Inf) 
colorsprop <- c('white',brewer.pal(10,"RdYlBu")[6] )
legend_prop <- expression("< 5 %", "5-10 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult men: age4-3", legend_prop)

## Young women: sex6 and age4-1
estML <- datos$parado_16
intervals_prop <- c(0.05, 0.10, 0.15,  Inf)
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8)] )
legend_prop <- expression("5-10 %", "10-15 %", '15-20')
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Young women: age4-1", legend_prop)

## Middle age women: sex6 and age4-2
estML <- datos$parado_26
intervals_prop <- c(0.05, 0.10, 0.15,  Inf)
colorsprop <- c(brewer.pal(10,"RdYlBu")[c(6,7,8)] )
legend_prop <- expression("5-10 %", "10-15 %", '15-20') 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Mid adult women: age4-2", legend_prop)

## Adult women: sex6 and age4-3
estML <- datos$parado_36
intervals_prop <- c(0, 0.05,  Inf) 
colorsprop <- c('white',brewer.pal(10,"RdYlBu")[6] )
legend_prop <- expression("< 5 %", "5-10 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "Adult women: age4-3", legend_prop)


## RRMSE UNEMPLOYMENT RATES

datos <- read.table(file="mse_UNEMPRate_500.txt", header=TRUE, sep=';', dec='.')$rrmse_UNEMPRatePI

## Young age men: sex1 and age4-1
estML <- datos[1:50]
dom <- 1:50
intervals_prop <- c(0.08, 0.15,  Inf)
colorsprop <- c(brewer.pal(5,"YlOrRd")[c(1, 2, 3)] )
legend_prop <- expression('< 8%', "8-15 %", "15-22 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "RRMSE Young men: age4-1", legend_prop)

## Middle age men: sex1 and age4-2
estML <- datos[105:154]
dom <- 1:50
intervals_prop <- c(0.08, 0.15, 0.22, Inf)  
colorsprop <- c(brewer.pal(5,"YlOrRd")[c(1, 2, 3, 4)] )
legend_prop <- expression('< 8%', "8-15 %", "15-22 %", "22-30 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "RRMSE Mid adult men: age4-2", legend_prop)

## Adult men: sex1 and age4-3
estML <- datos[209:258]
dom <- 1:50
intervals_prop <- c(0.08, 0.15, 0.22, Inf) 
colorsprop <- c(brewer.pal(5,"YlOrRd")[c(1, 2, 3, 4)] )
legend_prop <- expression('< 8%', "8-15 %", "15-22 %", "22-30 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "RRMSE Adult men: age4-3", legend_prop)

## Young age women: sex6 and age4-1
estML <- datos[53:102]
dom <- 1:50
intervals_prop <- c(0.08, 0.15,  Inf) 
colorsprop <- c(brewer.pal(5,"YlOrRd")[c(1, 2, 3)] )
legend_prop <- expression('< 8%', "8-15 %", "15-22 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "RRMSE Young women: age4-1", legend_prop)

## Middle age women: sex6 and age4-2
estML <- datos[157:206]
dom <- 1:50
intervals_prop <- c(0.08, 0.15, Inf)
colorsprop <- c(brewer.pal(5,"YlOrRd")[c(1, 2, 3)] )
legend_prop <- expression('< 8%', "8-15 %", "15-22 %") 
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "RRMSE Mid adult women: age4-2", legend_prop)

## Adult women: sex6 and age4-3
estML <- datos[261:310]
dom <- 1:50
intervals_prop <- c(0.08, 0.15, Inf)
colorsprop <- c(brewer.pal(5,"YlOrRd")[c(1, 2, 3)] )
legend_prop <- expression('< 8%', "8-15 %", "15-22 %")
result <- GroupClassification(dom, estML, intervals_prop)
PrintSpainMap(pathmap, result, colorsprop, "RRMSE Adult women: age4-3", legend_prop)
