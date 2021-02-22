
library(TCGAbiolinks)
library(dplyr)
library(RColorBrewer)
library(SummarizedExperiment)
library(ggplot2)
library(ggthemes)
library(ComplexHeatmap)

#-------------------------------
# Obtaining DNA methylation data
#-------------------------------
#Breast cancer
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 27",
  legacy = TRUE,
  sample.type = c("Solid Tissue Normal", "Primary Tumor")
)

GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

write.csv(assay(met), "TCGA-BRCA_methylation.csv")
write.csv(colData(met), "TCGA-BRCA_methylation_clinical data.csv")

#----------------------------
# PCA visualization
#----------------------------

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

df_pca <- prcomp(assay(met))
df_out <- as.data.frame(df_pca$rotation)[,1:3]
df_out = data.frame(cbind(df_out, met@colData))

percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste0(colnames(df_out), " (", paste(as.character(percentage), " %", ")"))

variables = c("sample_type", "ajcc_pathologic_stage")
plots = lapply(variables, function(p){
  ggplot(df_out, aes(PC1, PC2))+
    labs(title=p, x=percentage[1], y = percentage[2], col=p)+
    geom_point(size=2, aes(col=factor(df_out[,p])))+
    ggthemes::scale_color_colorblind()+
    stat_ellipse(aes(col=df_out[,p]))+
    theme_minimal()
  
}
)
print(gridExtra::marrangeGrob(plots, nrow=length(variables)/2, ncol=2, 
                              newpage = F)
)


#-----------------------------------
# Differential Methylation Analysis
#-----------------------------------
dmc <- TCGAanalyze_DMC(
  data = met,
  groupCol = "sample_type", 
  group1 = "Primary Tumor", 
  group2 = "Solid Tissue Normal", 
  alternative = "two.sided",
  p.cut = 0.05,
  diffmean.cut = 0.2,
  adj.method = "BH",
  save = T,
  legend = "Tissue type",
  plot.filename = "TCGA-BRCA_methylation_TvsN.png",
  cores = 2 
)

#Correct methylation status
summary(factor(dmc$status))
summary(dmc$mean.Primary.Tumor.minus.mean.Solid.Tissue.Normal)
dmc$status[which(dmc$mean.Primary.Tumor.minus.mean.Solid.Tissue.Normal<0&dmc$status=="Hypermethylated in Primary Tumor")] <- "Hypomethylated in Primary Tumor"

#--------------------------
# DNA Methylation heatmap
#-------------------------
# get the probes that are Hypermethylated or Hypomethylated
probes <- rownames(dmc)[grep("hypo|hyper",dmc$status,ignore.case = TRUE)]
sig.met <- met[probes,]

ta <- HeatmapAnnotation(
  df = sig.met@colData[, c("sample_type", "ajcc_pathologic_stage", "race", "gender")],
  col = list(
    sample_type = c("Solid Tissue Normal" = "grey", "Primary Tumor" = "darkorange2"),
    gender = c("male" = "blue", "female" = "pink")
  )
)

# row annotation: add the status for Tumor in relation to Normal
ra = rowAnnotation(
  df = dmc[probes, "status", drop = F],
  col = list(
    "status" =
      c("Hypomethylated in Primary Tumor" = "orange",
        "Hypermethylated in Primary Tumor" = "darkgreen")
  ),
  width = unit(1, "cm")
)

heatmap  <- Heatmap(
  matrix = assay(sig.met),
  name = "DNA methylation",
  col = matlab::jet.colors(200),
  show_row_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = T,
  show_column_names = FALSE,
  bottom_annotation = ta,
  left_annotation = ra,
  column_title = "DNA Methylation"
) 
# Save to image
png("heatmap.png",width = 10, height = 10, res = 600, unit = "in")
draw(heatmap, annotation_legend_side =  "top")
dev.off()


#--------------------------
# create csv file
#-------------------------
sig.dmc = dmc[probes, ]
if(all.equal(rownames(sig.dmc), rownames(sig.met))){
csv.data = data.frame(rowRanges(sig.met), sig.dmc)
write.csv(csv.data, "TCGA-BRCA_DMP_TvsN.csv")
}



