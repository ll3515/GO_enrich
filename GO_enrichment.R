# R WebGestalt
library(colorout)
options(stringsAsFactors=F)
rm(list=ls())
library(R.helper)
source("~/Dropbox/tools/R_functions/Liisi_functions.R")

library("WebGestaltR")
setwd("/Users/ll3515/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/Rep_1p19q_rm/WGCNA/WebGestalt")

# Gene sets for which the enrichment is performed.
load("/Users/ll3515/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/Rep_1p19q_rm/WGCNA/objects/WGCNA_MODULES_GEXPR_.merge-height=0.15.genes=13284.samples=119_pw7_bicor.R")


module = "M4"

geneFile<-paste0("/Users/ll3515/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/Rep_1p19q_rm/WGCNA/modules/",names(module_ensg_list)[grep(module,names(module_ensg_list))],"__power7.txt"); geneFile
refFile<-"/Users/ll3515/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/Rep_1p19q_rm/WGCNA/modules/bkgrnd__power7.txt"
outputDirectory<-getwd()


# run gene set enrichment analysis

enrichResult<-WebGestaltR(enrichMethod="ORA",organism="hsapiens",
	 enrichDatabase="geneontology_Biological_Process",interestGeneFile=geneFile,
     interestGeneType="genesymbol",referenceGeneFile=refFile,
     referenceGeneType="genesymbol",is.output=TRUE,
     outputDirectory=outputDirectory,projectName=module)


# read in the enrichment results
path = paste0("/Users/ll3515/Dropbox/Glioma_project2/TCGA_AA_DA/IDH1_ASTROCYTOMAS/Rep_1p19q_rm/WGCNA/WebGestalt/Project_",module,"/enrichment_results_",module,".txt")

test = read.delim(path)

Head(test)


# Separate gene IDs into lists

table(duplicated(test$description))

rownames(test)<-test$description

GO = test$description

GOlist = list()

for (i in GO){

	GOlist[[i]] = unlist(strsplit(test[i,"OverlapGene_UserID"],";"))
}

# This is the key bit of analysis
#---------------------------------------------------------------------------
# Calculate distance between GO terms
#---------------------------------------------------------------------------

#https://www.statisticshowto.datasciencecentral.com/jaccard-index/
#Jaccard Index = (the number in both sets) / (the number in either set) * 100
#J(X,Y) = |X∩Y| / |X∪Y|

# Count the number of members which are shared between both sets.
# Count the total number of members in both sets (shared and un-shared).
# Divide the number of shared members (1) by the total number of members (2).
# Multiply the number you found in (3) by 100.

#Ja.b = length(intersect(a,b))/length(union(a,b))


JMatrix = matrix(ncol = length(names(GOlist)), nrow =length(names(GOlist)))
	rownames(JMatrix) <- names(GOlist); colnames(JMatrix) <- names(GOlist)

for(i in names(GOlist)){
	for(j in names(GOlist)){

		a = GOlist[[i]]
		b = GOlist[[j]]

		JMatrix[i,j] = length(intersect(a,b))/length(union(a,b))

	}
}

Head(JMatrix)

# Each enrichment needs individual evaluation in terms of cutting the tree and choosing the representative terms
# We can think of a way to automatise this, but coming up with a clever solution would have taken more time than
# analysing the results individually.

# M1

d <- dist(JMatrix)
hc <- hclust(d, method = "ward.D2")
plot(hc, horiz = TRUE, cex = 0.8)

groups = cutree(hc, h = 5)
table(groups)

groups[groups == 1]

categories = c("mitochondrial translation","cellular respiration","proteasomal protein catabolic process","ribonucleoside metabolic process")


M1 = test[rownames(test)%in%categories,c("description","FDR")]
write.delim(M1, file = "M1_GO.txt")

# M2

d <- dist(JMatrix)
hc <- hclust(d, method = "ward.D2")
plot(hc, horiz = TRUE, cex = 0.5)

groups = cutree(hc, h = 10)
table(groups)

groups[groups == 14]

categories = c("inflammatory response","adaptive immune response","leukocyte chemotaxis","lymphocyte activation","macrophage activation"
	,"response to interferon-gamm","T-helper cell differentiation","leukocyte degranulation","actin cytoskeleton organization"
	,"interleukin-1 production","immunological synapse formation","cellular calcium ion homeostasis","positive regulation of JAK-STAT cascade"
	,"actin filament depolymerizatio")

M2 = test[rownames(test)%in%categories,c("description","FDR")]

write.delim(M2, file = "M2_GO.txt")

# etc..

#---------------------------------------------------------------------------
# Plotting
#---------------------------------------------------------------------------

M1 = read.delim("M1_GO.txt")
M1 = M1[c(1:4),]
M1$FDR[M1$FDR == 0]<-"3.2e-09"
M1$FDR = as.numeric(M1$FDR)
M1$logFDR = as.numeric(-log10(M1$FDR))

pdf("enrich_plot_M9_noaxis.pdf", h = 2, w = 4)
ggplot(data = M1, aes(x = description, y = logFDR)) +
geom_bar(stat="identity", width = 0.3,fill="#1e2e78") +
coord_flip(ylim = c(0,14)) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank())
dev.off()

# you can add text to the plot direclty, or annotate later manually
# I chose to do the latter, as adding text to the plot was too fidly and I was under time pressure
