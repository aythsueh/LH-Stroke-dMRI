# AZ stroke subjects
# (partial) Correlations b/n stroke survivors's brain measures indices and the behavioral cognitive performance
# DTI data: White matter tracts

rm(list=ls())


# DTI data read-in
rhTract <- read.table('AZstroke_rh_Tracts.csv', sep=",", header=TRUE)


# clean data
# adjustment: MD*1000
AZrhTract$MD <- AZrhTract$MD*1000

# adjustment: units of tract volume: TractVol*1.4*1.4*3.0 (voxel size:mm^3)
AZrhTract$TractVol <- AZrhTract$TractVol*1.4*1.4*3.0

# adjustment: normalized tract volume=volume/eTIV
AZrhTract$TractVol <- AZrhTract$TractVol/AZrhTract$eTIV

# adjustment: streamlines=WayTotal/eTIV
AZrhTract$WayTotal <- AZrhTract$WayTotal/AZrhTract$eTIV

# check the data 
head(rhTract)


# Lesion size data read-in
Lesion_size <- read.table('Lesion_size.csv', sep=",", header=TRUE)
names(Lesion_size)[1] <- "Subject3"

# Behavioral data read-in
CogBeh <- read.table('AZstroke_behavioral_fsl.csv', sep=",", header=TRUE)
names(CogBeh)[1] <- "Subject2"


# Combine all data into a table
RH <- cbind(rhTract, CogBeh, Lesion_size, row.names=NULL)

# adjustment: normalized tract volume=volume/eTIV
RH$Lesion_size <- RH$Lesion_size/RH$eTIV
RH$Lesion_size <- RH$Lesion_size*1000

RH[, c(9, 18)] <- NULL


# separate each ROI into different tables using list
ROI <- vector("list", 9)
for (i in 1:9){
  ROI[[i]] <- RH[RH$StructNum==i,]
}




#### 2021/1/17 Multiple regression analysis
# DVs: WMI,PSI,SingleWordComp,PicDesTask
# IVs: FA, MD, WayTotal
# Covariates: Age, Education, Lesion size
# also adjust the p value

library(ppcor)


## the effects of time post-stroke
names(CogBeh)
Covariates <- CogBeh[,c(2:4,9,5:8)]

library("PerformanceAnalytics")
chart.Correlation(Covariates)

Les_eff <- ROI[[1]][,c(17,12:15)]
chart.Correlation(Les_eff)


## FA calculated by probabilistic fiber tracking
## mean FA values of WM tracts
# r value
FA_r <- data.frame(matrix(0, nrow=9, ncol=5))
names(FA_r) <- c("ROI", "WMI","PSI","SingleWordComp","PicDesTask")
FA_r[,1] <- unique(rhTract$TractName)
FA_r

# p-value
FA_p <- data.frame(matrix(0, nrow=9, ncol=5))
names(FA_p) <- c("ROI", "WMI","PSI","SingleWordComp","PicDesTask")
FA_p[,1] <- unique(rhTract$TractName)
FA_p


for (i in 1:9) {
  FA_WMI <- na.omit(ROI[[i]][,c(4,12,9,17)])
  pcor_FA_WMI <- pcor.test(x=FA_WMI$FA, y=FA_WMI$WAIS_WMI, 
                           z=FA_WMI[,c("Age","Lesion_size")])
  FA_r[i,2] <- pcor_FA_WMI$estimate
  FA_p[i,2] <- pcor_FA_WMI$p.value  
}
for (i in 1:9) {
  FA_PSI <- na.omit(ROI[[i]][,c(4,13,9,17)])
  pcor_FA_PSI <- pcor.test(x=FA_PSI$FA, y=FA_PSI$WAIS_PSI, 
                           z=FA_PSI[,c("Age","Lesion_size")])
  FA_r[i,3] <- pcor_FA_PSI$estimate
  FA_p[i,3] <- pcor_FA_PSI$p.value
}
for (i in 1:9) {
  FA_Com <- na.omit(ROI[[i]][,c(4,14,9,17)])
  pcor_FA_Com <- pcor.test(x=FA_Com$FA, y=FA_Com$SingleWordComp, 
                           z=FA_Com[,c("Age","Lesion_size")])
  FA_r[i,4] <- pcor_FA_Com$estimate
  FA_p[i,4] <- pcor_FA_Com$p.value
}
for (i in 1:9) {
  FA_Pic <- na.omit(ROI[[i]][,c(4,15,9,17)])
  pcor_FA_Pic <- pcor.test(x=FA_Pic$FA, y=FA_Pic$PicDesTask, 
                           z=FA_Pic[,c("Age","Lesion_size")])
  FA_r[i,5] <- pcor_FA_Pic$estimate
  FA_p[i,5] <- pcor_FA_Pic$p.value
}
FA_r
FA_p

## Adjust P-values for Multiple Comparisons
# false discovery rate (FDR) 
FA_p_adj_1 <- FA_p[,2:5]
FA_p_adj_2 <- c(FA_p_adj_1[,1], FA_p_adj_1[,2], FA_p_adj_1[,3], FA_p_adj_1[,4])
FA_p_adj_3 <- p.adjust(FA_p_adj_2, method = "fdr")
FA_p_adj_4 <- data.frame(ROI = unique(rhTract$TractName),
                         WMI = FA_p_adj_3[1:9], 
                         PSI = FA_p_adj_3[10:18], 
                         SingleWordComp = FA_p_adj_3[19:27], 
                         PicDesTask = FA_p_adj_3[28:36])

# output as csv files
FA_r[,2:5] <- round(FA_r[,2:5], 4)
write.csv(FA_r, file="Strokes_pcrr_Tracts_FA.csv")
FA_p_adj_4[,2:5] <- round(FA_p_adj_4[,2:5], 3)
write.csv(FA_p_adj_4, file="Strokes_pcrr_Tracts_FA_p.csv")

  
##  I did the same data processing procedure for MD, WayTotal, and Tract Volume as well. 


## generate correlation chart for 9 ROIs 
library("PerformanceAnalytics")
ROI1 <- ROI[[1]][,c(6:9,11,16,17,20,25,26)]
chart.Correlation(ROI1, histogram=TRUE, pch=19)

chart.Correlation(ROI[[9]][,c(4:7)])


## check the number in each behaviral measure
length(na.omit(ROI[[1]]$WAIS_WMI))
length(na.omit(ROI[[1]]$WAIS_PSI))
length(na.omit(ROI[[1]]$SingleWordComp))
length(na.omit(ROI[[1]]$PicDesTask))
