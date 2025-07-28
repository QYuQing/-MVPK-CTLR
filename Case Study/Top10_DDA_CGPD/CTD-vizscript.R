# Load packages
library(data.table)
library(circlize)

# Import data
tetramers <- read.csv("F:\\KronRLS\\CaseStudy\\DDA\\Trovafloxacin_liver neoplasms, experimental\\CTD_combined.csv") # add you file name in between the quotation marks

# Prep Data
## make data frames for each variable
chemicals <- data.frame(tetramers$Chemical)
genes <- data.frame(tetramers$Gene)
phenotypes <- data.frame(tetramers$Phenotype.ID)
diseases <- data.frame(tetramers$Disease)

## Combine into dimer relationships
chem_gene <- cbind(chemicals, genes)
colnames(chem_gene) <- c("start", "end")

gene_pheno <- cbind(genes, phenotypes)
colnames(gene_pheno) <- c("start", "end")

pheno_disease <- cbind(phenotypes, diseases)
colnames(pheno_disease) <- c("start", "end")

## Combine into one data frame
cgpd <- rbind.data.frame(chem_gene, gene_pheno, pheno_disease)

## Add count values
df_counts <- setDT(cgpd)[,list(Count=.N),names(cgpd)]

## Lists for setting colors
chem_unique <- unique(chemicals)
genes_unique <- unique(genes)
pheno_unique <- unique(phenotypes)
disease_unique <- unique(diseases)

chem_colors <- colorRampPalette(c("blue4", "cyan")) (n = nrow(chem_unique))
gene_colors <- colorRampPalette(c("darkolivegreen", "seagreen1")) (n = nrow(genes_unique))
phenotype_colors <- colorRampPalette(c("purple4", "plum1")) (n = nrow(pheno_unique))
disease_colors <- colorRampPalette(c("firebrick1", "firebrick4")) (n = nrow(disease_unique))

node_colors <- c(chem_colors, gene_colors, phenotype_colors, disease_colors)
# ---------------------
# 输出 PDF 文件
# ---------------------
pdf("F:\\KronRLS\\CaseStudy\\DDA\\Trovafloxacin_liver neoplasms, experimental\\CTD_combined.pdf", width = 12, height = 12)
# Generate Chord Diagram with circlize package
circos.clear() # use this to reset all parameters

circos.par(gap.degree = 1) # set gap spacing between two neighboring nodes

## Make chord diagram
chordDiagram(df_counts,
             grid.col = node_colors,
             annotationTrack = "grid",
             preAllocateTracks = 1,
             annotationTrackHeight = 0.02)

## Add labels to diagram
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.5)#cex修改字体大小
}, bg.border = NA) 



# 关闭 PDF 输出
dev.off()