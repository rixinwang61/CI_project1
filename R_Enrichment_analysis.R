#gene enrichment analysis via ActivaPathways
setwd('/Users/wangrixin/Desktop/LUAD/Active_pathways')
library(ActivePathways)
GMT <- read.GMT('c5.go.bp.v7.5.1.symbols.gmt')
#in orginal GMT file, there is a http link in the column of name
#to replace the link with a name, run following code
gmt <- GMT
for (i in 1:(length(GMT))) {
  gmt[[i]]$name <- GMT[[i]]$id
}
#remove pathways (gene sets) that have less than 10 or more than 500 annotated genes
#gmt <- Filter(function(term) length(term$genes) >= 10, gmt)
#gmt <- Filter(function(term) length(term$genes) <= 500, gmt)

scores <- read.csv("scores.csv",header = T,row.names = 1)
scores <- as.matrix(scores)
enriched_pathways <- ActivePathways(scores, gmt) 
res <- ActivePathways(scores, GMT, cytoscape.file.tag = "enrichmentMap__")
export_as_CSV(enriched_pathways, "enriched_pathways.csv")
temp <- temp <- read.GMT('enrichmentMap__pathways.gmt')
for (i in 1:(length(temp))) {
  temp[[i]]$name <- temp[[i]]$id
}
write.GMT(temp,file = 'enrichmentMap__pathways.gmt')

#gene enrichment analysis via ReactomePA
library(org.Hs.eg.db)
library(clusterProfiler)
library(GOplot)
library(ggplot2)
library(stringr)
library(ReactomePA)
setwd('/Users/wangrixin/Desktop/LUAD')
gene <- read.csv('TCGA_gene_list.txt',header = F)
gene <- as.character(gene[,1])
gene <- bitr(gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
#enrichKEGG
ego <- enrichKEGG(gene = gene$ENTREZID,organism = 'hsa',keyType = 'kegg',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05)
ego2 <- setReadable(ego,OrgDb = 'org.Hs.eg.db',keyType = 'ENTREZID')
png('TCGA_kegg.png',width=15,height=12,units='in',res=600)
cnetplot(ego2,showCategory =10,categorySize='Count',circular=TRUE,colorEdge=TRUE)
dev.off()
#enrichGO
go <- enrichGO(gene = gene$ENTREZID,OrgDb = 'org.Hs.eg.db',ont = 'BP',pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,readable = T)
png('TCGA_go.png',width=15,height=12,units='in',res=600)
cnetplot(go,showCategory =10,categorySize='Count',circular=TRUE,colorEdge=TRUE)
dev.off()

#enrichPathway
pathway <- enrichPathway(gene$ENTREZID,organism = 'human',pAdjustMethod = 'BH',pvalueCutoff = 0.05,readable = TRUE)
png('Disgenet_pathway.png',width=20,height=14,units='in',res=600)
cnetplot(pathway,showCategory =10,categorySize='Count',circular=TRUE,colorEdge=TRUE)
dev.off()

