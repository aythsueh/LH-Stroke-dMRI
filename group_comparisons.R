# AZ stroke subjects and Controls
# Comparisons in each ROI the brain measures between the stroke and control groups at group-level
# DTI data: White matter tracts

rm(list=ls())

## AZ Stroke Subjects

# MRI data read-in
AZrhTract <- read.table('AZstroke_rh_Tracts.csv', sep=",", header=TRUE)

# clean data
# adjustment: MD*1000
AZrhTract$MD <- AZrhTract$MD*1000

# adjustment: units of tract volume: TractVol*1.4*1.4*3.0 (voxel size:mm^3)
AZrhTract$TractVol <- AZrhTract$TractVol*1.4*1.4*3.0

# adjustment: normalized tract volume=volume/eTIV
AZrhTract$TractVol <- AZrhTract$TractVol/AZrhTract$eTIV

# adjustment: streamlines=WayTotal/eTIV
AZrhTract$WayTotal <- AZrhTract$WayTotal/AZrhTract$eTIV


# Behavioral data read-in
AZCogBeh <- read.table('AZstroke_behavioral_fsl.csv', sep=",", header=TRUE)
AZCogBeh[, c(5:16)] <- NULL
names(AZCogBeh)[1] <- "Subject2"


# Combine all data into a table
AZRH <- cbind(AZrhTract, AZCogBeh, row.names=NULL)
AZRH[, 9] <- NULL
head(AZRH)


# separate each ROI into different tables using list
AZROI <- vector("list", 9)
for (i in 1:9){
  AZROI[[i]] <- AZRH[AZRH$StructNum==i,]
}

# check number of subjects
length(unique(AZROI[[1]]$Subject))


## I did the same data processing procidure for our control subjects




library(psych)

# define a function for calculating standard error
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

## FA
Group_Comp_FA <- data.frame(ROI=unique(AZrhTract$TractName), 
                            Strokes=rep(0,9),
                            S.stderr=rep(0,9),
                            Controls=rep(0,9),
                            C.stderr=rep(0,9),
                            p.values=rep(0,9),
                            adj.p=rep(0,9))
Group_Comp_FA


## FA calculated by probabilistic fiber tracking
for (i in 1:9) {
  FA_ttest <- t.test(AZROI[[i]]$FA, NCROI[[i]]$FA)
  Group_Comp_FA[i,c(2,4)] <- FA_ttest$estimate
  Group_Comp_FA[i,3] <- stderr(AZROI[[i]]$FA)
  Group_Comp_FA[i,5] <- stderr(NCROI[[i]]$FA)
  Group_Comp_FA[i,6] <- FA_ttest$p.value
}
Group_Comp_FA


## Adjust P-values for Multiple Comparisons
# false discovery rate (FDR) 
Group_Comp_FA[,6]
Group_Comp_FA[,7] <- p.adjust(Group_Comp_FA[,6], method = "fdr")

Group_Comp_FA[,c(2:7)] <- round(Group_Comp_FA[,c(2:7)], 3)
Group_Comp_FA


## Output
write.csv(Group_Comp_FA, file="Group_Comp_FA.csv")


##  I did the same data processing procedure for MD, WayTotal, and Tract Volume
## so :) 
