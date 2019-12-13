#https://malfaro2.github.io/SP1649/SP1649_Semana13_2019#1
#https://malfaro2.github.io/SP1649/SP1649_Semana13.1_2019#1
#https://rspatial.org/raster/analysis/3-spauto.html
#https://malfaro2.github.io/SP1649/lab14.html
#https://malfaro2.github.io/SP1649/SP1649_Semana14_2019#1
#https://www.sciencedirect.com/science/article/pii/S2211675317302725
library(rgdal)
distr_shp <- readOGR("distritos_2011.shp")
plot(distr_shp)
distr_shp@data$HOG_NBI

library(raster)
if (!require("rspatial")) devtools::install_github('rspatial/rspatial')
library(rspatial)
library(RColorBrewer)
datos <- data.frame(distr_shp@data)
str(datos)

library(tmap)
distr_shp@data$NBI <- (100 - (distr_shp@data$HOG_NBI))
tmap_mode("view")
tm_shape(distr_shp) + tm_polygons("HOG_NBI", palette= "RdYlBu", n = 10, title = "Porcentaje de hogares con NBI") +
  tm_layout(main.title = "Costa Rica: Distritos", main.title.position = "center", main.title.size = 1.3) + 
  tm_credits("Fuente: Instituto Nacional de Estadística y Censos (INEC), 2011", 
             position = c("center", "BOTTOM"), size = 5)  
tmap_mode("plot")
?tm_polygons
?palette
spplot(distr_shp,"HOG_NBI", 
       at=seq(min(distr_shp@data$HOG_NBI,na.rm=TRUE),
              max(distr_shp@data$HOG_NBI,na.rm=TRUE),
              length=10),col.regions=rev(brewer.pal(10,"RdYlBu")))
display.brewer.all()
library(maptools)
library(spdep)
library(leaflet)

###############modelo lineal simple
###############modelo lineal simple
###############modelo lineal simple
#
mod <- lm(datos$HOG_NBI ~ datos$REL_SEC12 + datos$IND_SSJEF  + datos$REL_DEPDEM)
summary(mod)
hist(mod$residuals, main = "Histograma de Residuos", xlab = "Residuos", ylab = "Frecuencia")
mod$residuals
library(lmtest)
dwtest( mod, alternative = "two.sided") #se asume correlación de los 
#residuos estudentizados
bptest(mod) #se asume homocedasticidad
plot(mod$residuals,datos$HOG_NBI)
plot(datos$HOG_NBI,mod$fitted.values)
plot(mod$residuals, mod$fitted.values)


#######Matriz de vecindades
#######Matriz de vecindades
#######Matriz de vecindades
coords<-coordinates(distr_shp)
IDs <- row.names(distr_shp)
####Criterio de la reina
mfrow=c(1,2)
list.queen<-poly2nb(distr_shp, queen=TRUE, snap = 1350)
plot(distr_shp, border="grey60", axes=TRUE, main="Criterio de la Reina")
plot(list.queen, coordinates(distr_shp), pch=19, cex=0.6, add=TRUE)
###Criterio de la torre
list.torre<-poly2nb(distr_shp, queen=FALSE)
plot(distr_shp, border="grey60", axes=TRUE, main="Criterio de la Torre")
plot(list.torre, coordinates(distr_shp), pch=19, cex=0.6, add=TRUE)
####Se elige el criterio de la reina debido a la naturaleza del NBI
#####Escoger la estructura espacial


W<-nb2listw(list.queen, style="W", zero.policy=TRUE)
B<-nb2listw(list.queen, style="B", zero.policy=TRUE)
S<-nb2listw(list.queen, style="S", zero.policy=TRUE)
dist <- nbdists(list.queen, coords)
?nbdists
idw <- lapply(dist, function(x) 1/(x/1000))
#idw <- lapply(dsts, function(x) 1/(x/1000))
W1<-nb2listw(list.queen, style="W", glist=idw, zero.policy=TRUE)
#plot(W1,coords)


pal <- brewer.pal(9, "Reds")
oopar <- par(mfrow=c(1,3), mar=c(1,1,3,1)+0.1)
z <- t(listw2mat(W))
brks <- c(0,0.1,0.143,0.167,0.2,0.5,1)
nbr3 <- length(brks)-3
image(1:472, 1:472, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
      main="Estilo W", axes=FALSE, cex=10)
box()
z <- t(listw2mat(B))
image(1:472, 1:472, z[,ncol(z):1], col=pal[c(1,9)], main="Estilo B",
      axes=FALSE, cex=15)
box()
z <- t(listw2mat(W1))
brks <- c(0.001,0.035,0.073,0.093,0.12,0.26,0.5)
image(1:472, 1:472, z[,ncol(z):1], col=pal[c(1,(9-nbr3):9)], breaks=brks,
      main="Estilo IDW W", axes=FALSE)
box()


#########Test de Moran
moran.lmW<-lm.morantest(mod, W, alternative="two.sided")
print(moran.lmW)

moran.lmW1<-lm.morantest(mod, W1, alternative="two.sided")
print(moran.lmW1)

moran.lmB<-lm.morantest(mod, B, alternative="two.sided")
print(moran.lmB)

moran.lmS<-lm.morantest(mod, S, alternative="two.sided")
print(moran.lmS)
####si hay autocorrelación espacial

#Test del Multiplicador de Lagrange
LM<-lm.LMtests(mod, W1, test="all")
print(LM)
#esta prueba sugiere que se debe usar el modelo de Error Espacial

#####Modelo de error Espacial
##modelo de máxima verosimilitud
#
Emvmod<-errorsarlm(HOG_NBI ~ REL_SEC12 + IND_SSJEF  + REL_DEPDEM,
                   data= distr_shp@data, W1)
summary(Emvmod)
?errorsarlm
#Modelo con estimación de mínimos cuadrados Generalizados
# 
Emcgmod <- GMerrorsar(HOG_NBI ~ REL_SEC12 + REL_DEPDEM+ IND_SSJEF,
                        data= distr_shp@data, W1)
summary(Emcgmod)

class(Emvmod$residuals)
hist(Emvmod$residuals)
plot(Emvmod$residuals)
max((Emvmod$residuals))
shapiro.test(Emvmod$residuals)
#dwtest(Emvmod, alternative = "two.sided") #se asume correlación de los 
#residuos estudentizados
#bptest(Emvmod) #se asume homocedasticidad
plot(datos$HOG_NBI,Emvmod$residuals)


######probar modelos con rezagos
#MODELOS SAR AUTO REGRESIVOS
#y= pWy+Xß+e
?lagsarlm
Lagmod<-lagsarlm(HOG_NBI ~ REL_SEC12 + IND_SSJEF  + REL_DEPDEM,
                  data= distr_shp@data, W1)
summary(Lagmod)
#Modelo GAM ADITIVO
library(mgcv)
distr_shp@data$x1<-coordinates(distr_shp)[,1]/1000
distr_shp@data$y1<-coordinates(distr_shp)[,2]/1000
Gammod<- gam(HOG_NBI ~ REL_SEC12 + IND_SSJEF  + REL_DEPDEM +s(x1,y1),
             data=distr_shp@data, family ="poisson" )
?gam
summary(Gammod)

autocar<-spautolm(HOG_NBI ~ REL_SEC12 + IND_SSJEF  + REL_DEPDEM ,
                data=distr_shp@data, family="CAR",
                listw=W1)
summary(autocar)

autosar<-spautolm(HOG_NBI ~ REL_SEC12 + IND_SSJEF  + REL_DEPDEM ,
                  data=distr_shp@data, family="SAR",
                  listw=W1)
summary(autosar)
autocar$fit$residuals


######Residuos
######Guardar los residuos del modelo
distr_shp@data$mod <- mod$residuals     #residuales Gaussiano
distr_shp@data$emvmod <- Emvmod$residuals #residual Error
distr_shp@data$lagmod <- Lagmod$residuals     #Sar con rezago
distr_shp@data$autosar <- autosar$fit$residuals  #autosar 
distr_shp@data$autocar <- autocar$fit$residuals  #autocar

par(mfrow=c(3,2))
spplot(distr_shp,c("mod", "emvmod","lagmod","autosar","autocar"),
       at=seq(min(distr_shp@data$mod,na.rm=TRUE),
              max(distr_shp@data$mod,na.rm=TRUE),
              length=12),col.regions=rev(brewer.pal(11,"RdBu")))
par(mfrow=c(1,1))
######Normalidad
library(EnvStats)
par(mfrow=c(3,2))
qqPlot(distr_shp@data$mod     ,      main="Modelo LM", add.line = TRUE)
qqPlot(distr_shp@data$emvmod  , main="Modelo LMerr", add.line = TRUE)
qqPlot(distr_shp@data$lagmod  , main="Modelo LagSar", add.line = TRUE)
qqPlot(distr_shp@data$autosar , main="Modelo AutoSar", add.line = TRUE)
qqPlot(distr_shp@data$autocar , main="Modelo AutoCar", add.line = TRUE)
par(mfrow=c(1,1))
?qqPlot

###Shapiro
shapiro.test(distr_shp@data$emvmod)
shapiro.test(distr_shp@data$lagmod)
shapiro.test(distr_shp@data$autosar)
shapiro.test(distr_shp@data$autocar)
###################AIC
AIC(mod)
AIC(Emvmod)
AIC(Lagmod)
AIC(autosar)
AIC(autocar)
