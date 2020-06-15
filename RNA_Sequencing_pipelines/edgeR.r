#This R script implements an edgeR pipeline (v.3.20.1), using the generalized linear model option,
#to identify genes differentially expressed with knockdown of 16p12.1 homologs in Drosophila.

#Quantified counts were generated from the Tophat2-HTSeq pipeline, and technical replicates were
#merged to produce a final set of quantified counts for each biological replicate.
#These data are available on NCBI GEO, accession number GSE GSE151330.

#Note that Sin (homolog of POLR3E) was sequenced in a separate batch with its own control (Elav-GAL4 at RT instead of 
#Elav-GAL4;UAS-Dicer2 at 25C), and therefore is analyzed separately from the other 3 homologs in the below pipeline.

#Load edgeR library
library("edgeR")

#Load raw counts file
biol_reps_data<-read.table("16p12_biol_reps.txt", sep="\t", stringsAsFactors=FALSE)
biol_reps_data_polr3e<-read.table("biol_reps_polr3e.txt", sep="\t", stringsAsFactors=FALSE)

#Move row names and column names out of data matrix
rownames(biol_reps_data)<-biol_reps_data[,1]
colnames(biol_reps_data)<-biol_reps_data[1,]
biol_reps_data<-biol_reps_data[,-1]
biol_reps_data<-biol_reps_data[-1,]

rownames(biol_reps_data_polr3e)<-biol_reps_data_polr3e[,1]
colnames(biol_reps_data_polr3e)<-biol_reps_data_polr3e[1,]
biol_reps_data_polr3e<-biol_reps_data_polr3e[,-1]
biol_reps_data_polr3e<-biol_reps_data_polr3e[-1,]

#Convert to matrix and convert characters to numeric type
biol_reps_matrix<-as.matrix(biol_reps_data)
class(biol_reps_matrix)<-"numeric"
biol_reps_matrix_polr3e<-as.matrix(biol_reps_data_polr3e)
class(biol_reps_matrix_polr3e)<-"numeric"

#Load data matrix into DGEList format to begin edgeR pipepline
groups<-c(rep("MOSMO",3),rep("UQCRC2",3),rep("CDR2",3),rep("Control",3))
y<-DGEList(biol_reps_matrix, group=groups)
groups_polr3e<-c(rep("POLR3E",3),rep("Control",3))
y_polr3e<-DGEList(biol_reps_matrix_polr3e, group=groups_polr3e)

#Filter variants: Recommended defaults here are 5 counts in smallest library, translated to counts/million 
#(here ~33 million=0.075 for POLR3E and ~16 million = 0.3 for rest),
#and in at least 3 samples (or expressed in all members of one group).
keep<-rowSums(cpm(y)>0.3)>=3
y<-y[keep, , keep.lib.sizes=FALSE] #Recalculate library sizes afterwards
keep_polr3e<-rowSums(cpm(y_polr3e)>0.075)>=3
y_polr3e<-y_polr3e[keep_polr3e, , keep.lib.sizes=FALSE]

#Calculate normalization factors
y<-calcNormFactors(y)
y_polr3e<-calcNormFactors(y_polr3e)

#Create MDS plots to display cliustering of samples relative to each other. Save as PDF.
plotMDS(y)
title(main="Multi-dimensional scaling plot of knockdown model sample clustering")
plotMDS(y_polr3e)
title(main="Multi-dimensional scaling plot of knockdown model sample clustering")

#Create the design matrix for all samples
design<-model.matrix(~ 0 + groups)
design_polr3e<-model.matrix(~ 0 + groups_polr3e)

#Estimate dispersion of variances along model
y <- estimateDisp(y, design, robust=TRUE)
y_polr3e <- estimateDisp(y_polr3e, design_polr3e, robust=TRUE)

#Model the quasi-lielihood dispersions of variation within model
fit <- glmQLFit(y, design, robust=TRUE)
fit_polr3e <- glmQLFit(y_polr3e, design_polr3e, robust=TRUE)

#Export normalized expression counts for downstream analysis
effectiveLibSize<-y$samples$lib.size*y$samples$norm.factors
normCounts<-log2(t(t(y$counts+0.5) / (effectiveLibSize+0.5)))
write.csv(normCounts,file="16p12_edgeR_norm_log_counts.csv")

effectiveLibSize_polr3e<-y_polr3e$samples$lib.size*y_polr3e$samples$norm.factors
normCounts_polr3e<-log2(t(t(y_polr3e$counts+0.5) / (effectiveLibSize_polr3e+0.5)))
write.csv(normCounts_polr3e,file="polr3e_edgeR_norm_log_counts.csv")

#Perform differential comparison tests. Controls are wild-type flies with same VDRC genetic background (GD or KK)
MOSMO<-glmQLFTest(fit, contrast=c(1,0,-1,0))
CDR2<-glmQLFTest(fit, contrast=c(0,1,-1,0))
UQCRC2<-glmQLFTest(fit, contrast=c(0,0,-1,1)))
POLR3E<-glmQLFTest(fit_polr3e, contrast=c(1,-1))

#Export results for each differential comparison to csv files
write.csv(topTags(MOSMO, n=nrow(MOSMO$table))$table, file="MOSMO_edgeR.csv")
write.csv(topTags(CDR2, n=nrow(CDR2$table))$table, file="CDR2_edgeR.csv")
write.csv(topTags(UQCRC2, n=nrow(UQCRC2$table))$table, file="UQCRC2_edgeR.csv")
write.csv(topTags(POLR3E, n=nrow(POLR3E$table))$table, file="POLR3E_edgeR.csv")
