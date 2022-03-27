rna.data <- read.table('gene_expression_rnaseq',sep='\t',header = T)
rownames(rna.data) <- rna.data[,1]
clin.data <- read.table('TCGA.LUAD.sampleMap%2FLUAD_clinicalMatrix',sep='\t',header = T)
rownames(clin.data) <- clin.data$sampleID
rownames(clin.data) <- gsub(rownames(clin.data),pattern="-",replace=".")
commonpatient <- intersect(rownames(clin.data),colnames(rna.data))
rna <- rna.data[,commonpatient]
clin <- clin.data[commonpatient,]
library(limma)
#normal tissue is 1, abnormal tissue is 0
design <- cbind(intercept=1,normal = as.numeric(clin$sample_type=='Solid Tissue Normal'))
rna.fit <- lmFit(rna,design)
rna.fit <- eBayes(rna.fit)
lm_fit_result <- topTable(rna.fit,coef = 2,number = 20530)
write.table(lm_fit_result,file = '/Users/wangrixin/Desktop/LUAD/R_data/lm_fit_results.csv',sep=',')

#survival analysis--rna expression vs os
surv.data <- read.table('survival%2FLUAD_survival.txt',sep='\t',header = T,row.names = 1)
rownames(surv.data) <- gsub(rownames(surv.data),pattern="-",replace=".")
os.event <- as.numeric(surv.data[colnames(rna.data),]$OS)
os.time <- surv.data[colnames(rna.data),]$OS.time
library(survival)
os <- Surv(os.time,os.event)
univariate.result <- array(NA,c(nrow(rna.data),4))
univariate.result <- as.data.frame(univariate.result)
colnames(univariate.result) <- c('Hazard_ratio','LCI','UCI','P.value')
rownames(univariate.result) <- rownames(rna.data)
for (i in 1:nrow(univariate.result)) 
{coxphmodel <- coxph(os~as.numeric(rna.data[i,]))
univariate.result$Hazard_ratio[i] <- summary(coxphmodel)$coef[1,2]
univariate.result$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
univariate.result$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
univariate.result$P.value[i] <- summary(coxphmodel)$coef[1,5]
}
univariate.result <- as.data.frame(univariate.result)
univariate.result$adj.pval <- p.adjust(univariate.result$P.value,method = "fdr")
univariate.result <- univariate.result[order(univariate.result$adj.pval,decreasing =FALSE),]
View(univariate.result[1:10,])
write.csv(univariate.result,file = '/Users/wangrixin/Desktop/LUAD/R_data/rna_os.csv')

#multivariate analysis
age <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
gender <- as.factor(clin$gender)
os.event2 <- as.numeric(surv.data[commonpatient,]$OS)
os.time2 <- surv.data[commonpatient,]$OS.time
os2 <- Surv(os.time2,os.event2)
multivariate.result <- array(NA,c(nrow(rna),4))
multivariate.result <- as.data.frame(multivariate.result)
colnames(multivariate.result) <- c('Hazard_ratio','LCI','UCI','P.value')
rownames(multivariate.result) <- rownames(rna)
for (i in 1:nrow(multivariate.result)) 
{coxphmodel2 <- coxph(os2~as.numeric(rna[i,])+age+as.numeric(gender))
multivariate.result$Hazard_ratio[i] <- summary(coxphmodel2)$coef[1,2]
multivariate.result$LCI[i] <- summary(coxphmodel2)$conf.int[1,3]
multivariate.result$UCI[i] <- summary(coxphmodel2)$conf.int[1,4]
multivariate.result$P.value[i] <- summary(coxphmodel2)$coef[1,5]
}
multivariate.result <- as.data.frame(multivariate.result)
multivariate.result$adj.pval <- p.adjust(multivariate.result$P.value,method = "fdr")
# Arrange the data frame by adjusted P.value
multivariate.result <- multivariate.result[order(multivariate.result$adj.pval,decreasing =FALSE),]
# View the top 10 genes associated with overall survival 
View(multivariate.result[1:10,])
write.csv(multivariate.result,file = '/Users/wangrixin/Desktop/LUAD/R_data/rna_age_gender_os.csv')


#survival analysis--cnv vs os
cnv.data <- read.table('CNV',sep = '\t',header = T,row.names = 1)


#rppa data
rppa.data <- read.table('TCGA.LUAD.sampleMap%2FRPPA_RBN',sep='\t',header=T,row.names = 1)
os.event2 <- as.numeric(surv.data[colnames(rppa.data),]$OS)
os.time2 <- surv.data[colnames(rppa.data),]$OS.time
os2 <- Surv(os.time2,os.event2)
univariate.result2 <- array(NA,c(nrow(rppa.data),4))
univariate.result2 <- as.data.frame(univariate.result2)
colnames(univariate.result2) <- c('Hazard_ratio','LCI','UCI','P.value')
rownames(univariate.result2) <- rownames(rppa.data)
for (i in 1:nrow(univariate.result2)) 
{coxphmodel <- coxph(os2~as.numeric(rppa.data[i,]))
univariate.result2$Hazard_ratio[i] <- summary(coxphmodel)$coef[1,2]
univariate.result2$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
univariate.result2$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
univariate.result2$P.value[i] <- summary(coxphmodel)$coef[1,5]
}
univariate.result2 <- as.data.frame(univariate.result2)
univariate.result2$adj.pval <- p.adjust(univariate.result2$P.value,method = "fdr")
univariate.result2 <- univariate.result2[order(univariate.result2$adj.pval,decreasing =FALSE),]
View(univariate.result2[1:10,])
write.csv(univariate.result2,file = '/Users/wangrixin/Desktop/LUAD/R_data/rppa_os.csv')



#mutation data
mutation.data <- read.csv('mc3_gene_level%2FLUAD_mc3_gene_level.txt',sep = '\t',header=T)
exclude <- as.numeric(is.na(mutation.data$sample))
mutation.data <- mutation.data[which(exclude==0),]
rownames(mutation.data) <- mutation.data[,1]
mutation.data <- mutation.data[,-1]
commonpatient2 <- intersect(rownames(clin.data),colnames(mutation.data))
clin2 <- clin.data[commonpatient2,]
mutation <- mutation.data[,commonpatient2]


