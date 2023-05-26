# AZ stroke subjects & healthy normal controls
# Demographics

rm(list=ls())

# Behavioral data read-in
AZ <- read.table('AZstroke_behavioral.csv', sep=",", header=TRUE)
NC <- read.table('NC_behavioral.csv', sep=",", header=TRUE)

# check number of the subjects for each group
length(unique(AZ$Subject))
length(unique(NC$Subject))



library("psych")

## 1. Age
describe(AZ$Age, na.rm=TRUE)
describe(NC$Age, na.rm=TRUE)

age.ttest <- t.test(AZ$Age, NC$Age)
age.ttest


## 2. Gender
describe(AZ$Gender, na.rm=TRUE)
describe(NC$Gender, na.rm=TRUE)

length(AZ$Gender[AZ$Gender==0]) # female
length(AZ$Gender[AZ$Gender==1]) # male

length(NC$Gender[NC$Gender==0]) # female
length(NC$Gender[NC$Gender==1]) # male

## FreeSurfer
# AZ 12 13
# NC 28 16

## FSL
# AZ 12 13
# NC 16 8

gender <- data.rbind(c(12, 13), c(28, 16))
gender.chisqtest <- chisq.test(gender, correct=FALSE)
gender.chisqtest


## 3. Education level
describe(AZ$Education, na.rm=TRUE)
describe(NC$Education, na.rm=TRUE)

educ.ttest <- t.test(AZ$Education, NC$Education)
educ.ttest


## 4. WAIS-4
# WMI
describe(AZ$WAIS_WMI, na.rm=TRUE)
describe(NC$WMI, na.rm=TRUE)

wmi.ttest <- t.test(AZ$WAIS_WMI, NC$WMI)
wmi.ttest

# PSI
describe(AZ$WAIS_PSI, na.rm=TRUE)
describe(NC$PSI, na.rm=TRUE)

psi.ttest <- t.test(AZ$WAIS_PSI, NC$PSI)
psi.ttest
