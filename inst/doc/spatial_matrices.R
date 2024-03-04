## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SDPDmod)

## ----eval=T, echo=T, warning=FALSE, message=FALSE, results='hide'-------------
library("sf")
ger <- st_read(system.file(dsn = "shape/GermanyNUTS3.shp",
                         package = "SDPDmod"))

data(gN3dist, package = "SDPDmod")

## -----------------------------------------------------------------------------
    W_1 <- mOrdNbr(ger) ## first order neighbors

## -----------------------------------------------------------------------------
   W_2n <- mOrdNbr(sf_pol = ger, m = 2) ## second order neighbors
   W_3n <- mOrdNbr(ger, 3) ## third order neighbors

## -----------------------------------------------------------------------------
   ls <- ger[which(substr(ger$NUTS_CODE,1,3)=="DE9"),] ## Lower Saxony districts
   W_len_sh <- SharedBMat(ls)

## -----------------------------------------------------------------------------
    W_knn <- mNearestN(distMat = gN3dist, m = 5) ## 5 nearest neighbors

## -----------------------------------------------------------------------------
    ## inverse distance no cut-off
    W_inv1 <- InvDistMat(distMat = gN3dist) 
    ## inverse distance with cut-off 100000 meters
    W_inv2 <- InvDistMat(distMat = gN3dist, distCutOff = 100000) 
    gN3dist2 <- gN3dist/1000 ## convert to kilometers
    ## inverse distance with cut-off 100 km
    W_inv3 <- InvDistMat(distMat = gN3dist2, distCutOff = 100)  
    ## inverse distance with cut-off 200km and exponent 2
    W_inv4 <- InvDistMat(gN3dist2, 200, powr = 2) 

## -----------------------------------------------------------------------------
    ## Exponential distance no cut-off
    W_exp1 <- ExpDistMat(distMat = gN3dist) 
    ## Exponential distance with cut-off 100000 meters
    W_exp2 <- ExpDistMat(distMat = gN3dist, distCutOff = 100000) 
    gN3dist2 <- gN3dist/1000 ## convert to kilometers
    ## Exponential distance with cut-off 100 km 
    W_exp3 <- ExpDistMat(gN3dist2, 100) 
    ## Exponential distance with cut-off 100 km
    W_exp4 <- DistWMat(gN3dist2, 100, type = "expo") 
    all(W_exp3==W_exp4)
    ## Exponential distance with cut-off 200 km and exponent 0.001
    W_exp5 <- ExpDistMat(gN3dist2, 200, expn = 0.001) 

## -----------------------------------------------------------------------------
    ## Double-Power distance no cut-off, exponent 2
    W_dd1 <- DDistMat(distMat = gN3dist) 
    ## Double-Power distance with cut-off 100000 meters, exponent 2
    W_dd2 <- DDistMat(distMat = gN3dist, distCutOff=100000) 
    gN3dist2 <- gN3dist/1000 ## convert to kilometers
    ## Double-Power distance with cut-off 100 km 
    W_dd3 <- DDistMat(gN3dist2, 100) 
    ## Double-Power distance with cut-off 100 km
    W_dd4 <- DistWMat(gN3dist2, 100, type = "doubled") 
    all(W_dd3==W_dd4) 
    ## Double-Power distance with cut-off 200km and exponent 3
    W_dd5 <- DDistMat(gN3dist2, 200, powr = 3) 

## -----------------------------------------------------------------------------
   W_2n_norm <- mOrdNbr(sf_pol = ger, m = 2, rn = T) ## second order neighbors
   W_2n_norm2 <- rownor(W_2n)
   all(W_2n_norm==W_2n_norm2)

## -----------------------------------------------------------------------------
  W_inv1_norm <- InvDistMat(distMat = gN3dist, mevn = T) ## Inverse distance
  W_inv1_norm2 <- eignor(W_inv1)
  all(W_inv1_norm==W_inv1_norm2)

