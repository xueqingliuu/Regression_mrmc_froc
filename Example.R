library(MASS)
library(RJafroc)
library(abind)

setwd("C:/Users/lloollmay/OneDrive - National University of Singapore/Desktop/PKU/group/Code")
source("Others_regression_cov.R")
source("Others_regression.R")
source("two_part_bootstrap.R")
#source("DBM.R")

####Transform the dataset
initial_data <- DfReadDataFile("4mm_all_AUC_new.xlsx", format = "JAFROC", delimiter = ",")
K1 <- dim(initial_data$ratings$NL)[3]-dim(initial_data$ratings$LL)[3]
K2 <- sum(initial_data$lesions$perCase)
dataset0 <- gen.data.1(initial_data,K1,K2)

nR = 1
nN = 161
nD = 183
nI = 2

dataset <- dataset0[dataset0[,5]==12,]
dataset[,5] <-1
### Regression
twopart_regression(dataset, nN, nD, nR, nI, nBoots=5000)


  # twopart_regression_1(dataset, nN, nD, nR, nI,K2)
# ### DBM
dataset <- cbind(dataset, 1*(dataset[,1]!=-Inf))
dataset <- cbind(dataset, 1*(dataset[,2]!=-Inf))
# PS_Sen <- pseudovalues_Sen(dataset, nN, nD, nR)
# PS_FPR <- pseudovalues_FPR(dataset, nN, nD, nR)
# PS_cafroc <- pseudovalues_cAFROC(dataset, nN, nD, nR)
# #PS_cafroc[PS_cafroc[,4]==9,2] <- NA
# DBM_Sen <- DBM_func(PS_Sen, nN, nD, nR)
# DBM_FPR <- DBM_func(PS_FPR, nN, nD, nR)
# DBM_cafroc <- DBM_func(PS_cafroc, nN, nD, nR)

###Covariates
library(readxl)
cov <- read_xlsx("casecov.xlsx")
id <- read_xlsx("caseid.xlsx")


c <- merge(id, cov, by.x = "case", by.y = "PatientID", all.x = TRUE)
cov_new <- c[order(c$new_caseID),]
d <- substr(cov_new$PatientAge,1,3)
cov_new$PatientAge <- as.numeric(d)

#dataset <- cbind(dataset, rep(cov_new$PatientSex,nR))
#dataset <- cbind(dataset, rep(cov_new$PatientAge,nR))
sex <- rep(cov_new$PatientSex,nR)
age <- rep(cov_new$PatientAge,nR)

dataset <- gen.data.1(initial_data,K1,K2)
dataset <- cbind(dataset, sex)
twopart_regression_sex(dataset, nN, nD, nR, nI, nBoots=5000)

dataset <- gen.data.1(initial_data,K1,K2)
dataset <- cbind(dataset, age)
twopart_regression_age(dataset, nN, nD, nR, nI, nBoots=5000)


#Figure
#TEST1
library(plotROC)
library(ggplot2)

len <- c()
dataset.t1 <- NULL
Reader <- NULL
for (i in 1:12) {
  len <- length(dataset[dataset[,5]==i&dataset[,7]==1,c(1,3)][,1])
  dataset.t1 <- rbind(dataset.t1, dataset[dataset[,5]==i&dataset[,7]==1,c(1,3)])
  Reader <- c(Reader, rep(paste("Reader",i) , len))
}
dataset.t1 <- data.frame(dataset.t1)
dataset00 <- cbind(dataset.t1, Reader)


cafrocplot1 <- ggplot(dataset00, aes(m = X1, d = X2, color=Reader))+ geom_roc(n.cuts=20,labels=FALSE,linetype=1,size=0.5,pointsize = 0.2)
cafrocplot1 + style_roc(xlab = "1 - Specificity", ylab="Sensitivity") + theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = 'none')


len <- c()
dataset.t2 <- NULL
Reader <- NULL
for (i in 1:12) {
  len <- length(dataset[dataset[,5]==i&dataset[,8]==1,c(2,3)][,1])
  dataset.t2 <- rbind(dataset.t2, dataset[dataset[,5]==i&dataset[,8]==1,c(2,3)])
  Reader <- c(Reader, rep(paste("Reader",i) , len))
}
dataset.t2 <- data.frame(dataset.t2)
dataset000 <- cbind(dataset.t2, Reader)


cafrocplot1 <- ggplot(dataset000, aes(m = X1, d = X2, color=Reader))+ geom_roc(n.cuts=20,labels=FALSE,linetype=1,size=0.5,pointsize = 0.2)
cafrocplot1 + style_roc(xlab = "1 - Specificity", ylab="Sensitivity") + theme(panel.grid.major = element_line(colour = NA),
                                                                              panel.background = element_rect(fill = "transparent",colour = NA),
                                                                              plot.background = element_rect(fill = "transparent", colour = NA),
                                                                              panel.grid.minor = element_blank(),
                                                                              legend.position = 'none')

###Plot LDR and NSMR
LDR.1 <- c(0.428, 0.504, 0.654, 0.448, 0.390, 0.461, 0.495, 0.637, 0.321, 0.568, 0.654, 0.386)
LDR.2 <- c(0.666, 0.788, 0.843, 0.741, 0.608, 0.716, 0.636, 0.770, 0.486, 0.732, 0.746, 0.621)
NSMR.1 <- c(0.019, 0.161, 0.137, 0.068, 0.025, 0.043, 0.019, 0.012, 0.025, 0.049, 0.068, 0.025)
NSMR.2 <- c(0.037, 0.137, 0.149, 0.037, 0.025, 0.006, 0.093, 0.069, 0, 0.043, 0.037, 0.013)

a <- data.frame(LDR=c(LDR.1, LDR.2),NSMR=c(NSMR.1,NSMR.2), Test=c(rep("Without CAD",12), rep("With CAD",12)))
ggplot(data=a, aes(x=NSMR, y=LDR, group=Test, color=Test))+
  geom_point()+
  geom_line()+
  xlab("NSMR")+
  ylab("LDR")+labs(color = "Test")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = element_blank(),
        legend.position = c(.85,.15),
        legend.box.background = element_rect(colour = "black"))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2))
