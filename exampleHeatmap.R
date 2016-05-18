library(ComplexHeatmap)
library(RColorBrewer)

my_palette <- colorRampPalette(c("darkblue","aliceblue","brown4"))(n = 299)
prefix <- file.path("C:","Users","djenk")

icbp_status <- read.table(file.path(prefix,"Dropbox","bild_signatures",
                                    "david_results","updated_fig2",
                                    "icbpstatus.txt"), row.names = 1, header=T,
                          sep="\t")
single_pathway_best_icbp <-read.table(file.path(prefix,"Dropbox","bild_signatures",
                                                "Datasets","optimized_single_pathway_icbp.txt"),
                                      stringsAsFactors=FALSE, header=1, row.names=1,sep='\t')

er_pr_her2_status2 <- readRDS(file.path(prefix,"Dropbox","bild_signatures",
                                        "david_results","updated_fig2",
                                        "tcga_er_pr_her2_formatted.rds"))

intrinsic <- read.table(file.path(prefix,"Dropbox","bild_signatures",
                                  "david_results","drug_response_heatmap",
                                  "intrinsic.txt"), row.names = 1, header=T,
                        sep="\t")

phenotype.annotation <- data.frame(Phenotype=c("HER2/IGF1R/AKT Phenotype",
                                               "BAD/EGFR/KRAS/RAF1 Phenotype",
                                               "BAD/EGFR/KRAS/RAF1 Phenotype",
                                               "HER2/IGF1R/AKT Phenotype",
                                               "HER2/IGF1R/AKT Phenotype",
                                               "BAD/EGFR/KRAS/RAF1 Phenotype",
                                               "BAD/EGFR/KRAS/RAF1 Phenotype",
                                               "BAD/EGFR/KRAS/RAF1 Phenotype"))

colnames(intrinsic) <- "Intrinsic.Subtype"
rownames(intrinsic)[50] <- "T47D.Kbluc"
intrinsic <- intrinsic[rownames(single_pathway_best_icbp),,drop=F]
intrinsic <- cbind(icbp_status, intrinsic)
set.seed(123)
kmeans.4 <- kmeans(scale(single_pathway_best_icbp)[,c(1,2,3,4)], 4,
                   iter.max=100000, nstart = 100)
kmeans.4.clusters <-data.frame(kmeans.4$cluster)

kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "3",] <- "EGFR/BAD low"
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "1",] <- "EGFR/BAD high"
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "4",] <- "AKT/HER2 high"
kmeans.4.clusters[kmeans.4.clusters$kmeans.4.cluster == "2",] <- "AKT/HER2 low"
colnames(kmeans.4.clusters) <- "KMeans.Clusters"
kmeans.4.clusters <- cbind(intrinsic, kmeans.4.clusters)

ha_row_icbp1 = HeatmapAnnotation(df = kmeans.4.clusters[,1:4,drop=F],
                                 name="Intrinsic Subtype",
                                 col = list(PR = c("Positive" =  "#4DAF4A",
                                                   "Negative" = "#984EA3",
                                                   "Unavailable" = "grey"),
                                            HER2 = c("Positive" =  "#FFFF33",
                                                     "Negative" = "#F781BF",
                                                     "Unavailable" = "grey"),
                                            ER = c("Positive" =  "#E41A1C",
                                                   "Negative" = "#377EB8",
                                                   "Unavailable" = "grey"),
                                            Intrinsic.Subtype = c("Basal" = brewer.pal(6, "Dark2")[1],
                                                                  "Claudin-low" = brewer.pal(6, "Dark2")[2],
                                                                  "HER2-Basal" = brewer.pal(6, "Dark2")[3],
                                                                  "HER2-Luminal" = brewer.pal(6, "Dark2")[4],
                                                                  "Luminal" = brewer.pal(6, "Dark2")[5],
                                                                  "Normal-like" = brewer.pal(6, "Dark2")[6],
                                                                  "Unavailable" = "grey")),
                                 which = "row", width = unit(1.333, "cm"), show_legend = F)

ha_row_icbp2 = HeatmapAnnotation(df = kmeans.4.clusters[,5,drop=F],
                                 name="K-Means Clusters",
                                 col = list(KMeans.Clusters = c("AKT/HER2 high" = brewer.pal(5, "Set1")[1],
                                                                "AKT/HER2 low" = brewer.pal(5, "Set1")[2],
                                                                "EGFR/BAD high" = brewer.pal(5, "Set1")[3],
                                                                "EGFR/BAD low" = brewer.pal(5, "Set1")[5])),
                                 which = "row", width = unit(0.333, "cm"))

topha_4 <- HeatmapAnnotation(df = phenotype.annotation[c(1:6,8),,drop=F],
                             col = list(Phenotype = c("HER2/IGF1R/AKT Phenotype" =  "coral3",
                                                      "BAD/EGFR/KRAS/RAF1 Phenotype" = "aquamarine4")),
                             height = unit(0.333, "cm"), show_legend = F)

h1 <- Heatmap(scale(single_pathway_best_icbp[,c(1:6,8)]), cluster_rows = T,
              cluster_columns = T, show_row_names = F, show_column_names = T,
              row_title_gp = gpar(fontsize =10), combined_name_fun = NULL,
              name="Scaled\nPathway\nActivity",col=my_palette,
              top_annotation = topha_4, column_title = "ICBP",
              show_heatmap_legend = F,
              column_dend_reorder = c(1,100,100,10,1,100,100))

draw(h1+ha_row_icbp1,row_dend_side = "left", annotation_legend_side = "bottom")





geneListBAD=as.matrix(c(bad_gene_list_DOWN, bad_gene_list_UP))

geneListAKT=as.matrix(c(akt_gene_list_DOWN, akt_gene_list_UP))

geneListHER2=as.matrix(c(her2_gene_list_DOWN, her2_gene_list_UP))

geneListIGF1R=as.matrix(c(igf1r_gene_list_DOWN, igf1r_gene_list_UP))

geneListKRASGV=as.matrix(c(krasgv_gene_list_DOWN, krasgv_gene_list_UP))

geneListRAF=as.matrix(c(raf_gene_list_DOWN, raf_gene_list_UP))

geneListEGFR=as.matrix(c(egfr_gene_list_DOWN, egfr_gene_list_UP))

geneListKRADQH=as.matrix(c(krasqh_gene_list_DOWN, krasqh_gene_list_UP))

