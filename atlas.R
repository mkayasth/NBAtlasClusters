library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(copykat)

# needs to be set for large dataset analysis
options(future.globals.maxSize = 6 * 1024^3) # setting max to 6 GB.

# loading the rds file from NBAtlas -- lightweight version.
atlas <- readRDS("NBAtlas/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds")

# first, lets run supervised data cluster from the scVI_umap that came along with the rds file.
dim(Embeddings(atlas, reduction = "scvi_umap"))
atlas <- FindNeighbors(atlas, reduction = "scvi_umap", dims = 1:2)  # scVI UMAP has 2 dimensions for 50k cells.
atlas <- FindClusters(atlas, resolution = 0.5)

table(atlas$orig.ident)
table(atlas$Cell_type)
# orig.ident and cell_type data match exactly in the .rds file (original annotation and annotation from scVI umap are the same).

# dim plot for the scvi umap obtained from the paper itself.
plot_supervised <- DimPlot(atlas, reduction = "scvi_umap", group.by = "orig.ident")
plot_supervised

# Now, running PCA, umap and seeing if the clusters still match.
atlas <- readRDS("NBAtlas/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds")
atlas2 <- NormalizeData(atlas)
atlas2 <- FindVariableFeatures(atlas2)

# unsupervised clustering of the cells.
atlas2 <- ScaleData(atlas2)
atlas2 <- RunPCA(atlas2)

dim(Embeddings(atlas2, reduction = "pca"))
atlas2 <- FindNeighbors(atlas2, dims = 1:50)
atlas2 <- FindClusters(atlas2, resolution = 0.5)
atlas2 <- RunUMAP(atlas2, dims = 1:50, return.model = T)


# Visualize UMAP based on supervised clustering.
plot_unsupervised <- DimPlot(atlas2, reduction = "umap", group.by = "orig.ident")

# clusters as per umap (unsupervised).
DimPlot(atlas2, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()

# unsupervised clusters but data labeled later from orig.identities.
plot_unsupervised


plot_supervised + plot_unsupervised

# Maybe batch correction is the issue here. We have data from 7 different studies.
atlas <- readRDS("NBAtlas/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds")
atlas2 <- NormalizeData(atlas)
atlas2 <- FindVariableFeatures(atlas2)
atlas2 <- ScaleData(atlas2)
atlas2 <- RunPCA(atlas2)

atlas2 <- RunHarmony(atlas2, group.by.vars = "Study")

# number of dimensions in this reduction.
atlas2@reductions$harmony


atlas2 <- FindNeighbors(atlas2, reduction = "harmony", dims = 1:50)
atlas2 <- FindClusters(atlas2, resolution = 0.5)

atlas2 <- RunUMAP(atlas2, reduction = "harmony", dims = 1:50)


plot_unsupervised_batchEffect <- DimPlot(atlas2, reduction = "umap", group.by = "orig.ident")
plot_unsupervised_batchEffect

# self label based on cluster.
DimPlot(atlas2, reduction = "umap", group.by = "seurat_clusters")

plot_supervised + plot_unsupervised_batchEffect

# just looking at specific CNAs inside neuroendocrine.
subset_neuroendocrine <- subset(atlas, subset = Cell_type == "Neuroendocrine")

subset_neuroendocrine <- NormalizeData(subset_neuroendocrine)
subset_neuroendocrine <- FindVariableFeatures(subset_neuroendocrine)


# sketching subset_neuroendocrine for copyKat. Cant load dense matrix into memory yet.
subset_neuroendocrine <- SketchData(
  object = subset_neuroendocrine,
  ncells = 500,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# laod the whole dense matrix into memory.
counts_matrix_neuroendocrine <- as.matrix(GetAssayData(subset_neuroendocrine, 
                                          layer = "counts"))

# run copycKAT.
copykat_results <- copykat(
  rawmat = counts_matrix_neuroendocrine,
  id.type = "S",        # S for gene symbols, E if Ensembl.
  sam.name = "MySample"
)


# output txt file has log2ratios for CNA. Greater than 0 = duplication; less than 0 = deletion.

# load excel file. My name for the copyKAT output is: MySample_copykat_CNA_raw_results_gene_by_cell.txt
copykat_result <- read.csv("MySample_copykat_CNA_raw_results_gene_by_cell.txt", sep = "\t", header = TRUE)


# Select gene name (hgnc_symbol) and all cell columns (everything after 'band')
subset_data <- copykat_result[, c("hgnc_symbol", names(copykat_result)[8:ncol(copykat_result)])]


# neuroblastoma in neuroendocrine cell has 1p loss and 17q gain.

# some chromosome 1p genes for loss calculation.

#1p locus genes.
genes_1p <- c("SAMD11", "NOC2L", "HES4", "SDF4", "B3GALT6", "UBE2J2", "ACAP3", "INTS11", "AURKAIP1", "CCNL2", "MRPL20", "SSU72", "MIB2", "GNB1", "PRKCZ", "FAAP20", "RER1", "LRRC47", "CHD5", "RPL22", "ACOT7", "CAMTA1", "PARK7", "RERE", "ENO1", "SLC25A33", "CLSTN1", "CTNNBIP1", "LZIC", "UBE4B", "KIF1B", "PGD", "DFFA", "CASZ1", "TARDBP", "SRM", "EXOSC10", "FBXO44", "MIIP", "VPS13D", "PRDM2", "SPEN", "FBXO42", "SZRD1", "NECAP2", "NBPF1", "ATP13A2", "SDHB", "ARHGEF10L", "UBR4", "MRTO4", "AKR7A2", "CAPZB", "CAMK2N1", "PINK1", "DDOST", "HP1BP3", "EIF4G3", "USP48", "CDC42", "KDM1A", "HNRNPR", "RPL11", "PITHD1", "LYPLA2", "PNRC2", "SRSF10", "SRRM1", "CLIC4", "RSRP1", "TMEM50A", "MTFR1L", "ZNF593", "SH3BGRL3", "HMGN2", "ARID1A", "TMEM222", "WASF2", "IFI6", "STX12", "DNAJC8", "PHACTR4", "TRNAU1AP", "SNHG12", "TAF12", "YTHDF2", "EPB41", "SRSF4", "SDC3", "PUM1", "NKAIN1", "SNRNP40", "ZCCHC17", "PEF1", "ADGRB2", "PTP4A2", "KPNA6", "TXLNA", "CCDC28B", "EIF3I", "HDAC1", "MARCKSL1", "BSDC1", "ZBTB8OS", "SYNC", "S100PBP", "AK2", "PHC2", "CSMD2", "SMIM12", "SFPQ", "ZMYM4", "KIAA0319L", "C1orf216", "AGO4", "AGO1", "AGO3", "ADPRHL2", "TRAPPC3", "MAP7D1", "THRAP3", "LSM10", "MRPS15", "MEAF6", "GNL2", "MANEAL", "YRDC", "C1orf122", "SF3A3", "UTP11", "AKIRIN1", "NDUFS5", "PABPC4", "PPIE", "CAP1", "PPT1", "RLF", "RIMS3", "NFYC", "CITED4", "CTPS1", "SCMH1", "HIVEP3", "FOXJ3", "PPCS", "PPIH", "YBX1", "SVBP", "EBNA1BP2", "MED8", "HYI", "PTPRF", "ST3GAL3", "ATP6V0B", "B4GALT2", "DMAP1", "ERI3", "RNF220", "RPS8", "EIF2B3", "UROD", "ZSWIM5", "MUTYH", "PRDX1", "AKR1A1", "GPBP1L1", "MAST2", "PIK3R3", "POMGNT1", "LRRC41", "UQCRH", "ATPAF1", "CMPK1", "SPATA6", "BEND5", "AGBL4", "ELAVL4", "FAF1", "RNF11", "EPS15", "OSBPL9", "NRDC", "RAB3B", "BTF3L4", "ZFYVE9", "PRPF38A", "ZYG11B", "SCP2", "MAGOH", "HSPB11", "LRRC42", "TMEM59", "MRPL37", "SSBP3", "USP24", "MYSM1", "JUN", "TM2D1", "PATJ", "USP1", "DOCK7", "JAK1", "DNAJC6", "LEPROT", "SGIP1", "MIER1", "SERBP1", "LRRC7", "LRRC40", "SRSF11", "ANKRD13C", "ZRANB2", "NEGR1", "TYW3", "SLC44A5", "ACADM", "RABGGTB", "ST6GALNAC5", "MIGA1", "FUBP1", "TTLL7", "PRKACB", "RPF1", "GNG5", "C1orf52", "DDAH1", "ODF2L", "SH3GLB1", "LMO4", "GTF2B", "ZNF326", "ZNF644", "BTBD8", "RPAP2", "RPL5", "MTF2", "DR1", "FNBP1L", "DNTTIP2", "ARHGAP29", "CNN3", "PTBP2", "DPYD", "MIR137HG", "SNX7", "PLPPR5", "PLPPR4", "TRMT13", "RTCA", "EXTL2", "DPH5", "RNPC3", "NTNG1", "PRPF38B", "WDR47", "TMEM167B", "GSTM3")
copykat_1p <- subset_data[subset_data$hgnc_symbol %in% genes_1p, ]
rownames(copykat_1p) <- copykat_1p$hgnc_symbol
copykat_1p <- copykat_1p[, -which(names(copykat_1p) == "hgnc_symbol")]

# some chromosome 17q genes for gain calculation.

# 17q locus genes.
genes_17q <- c("WSB1", "NLK", "TMEM97", "IFT20", "TNFAIP1", "POLDIP2", "TMEM199", "SARM1", "PIGS", "KIAA0100", "SDF2", "SUPT6H", "RPL23A", "TRAF4", "ERAL1", "FLOT2", "DHRS13", "PHF12", "NUFIP2", "GIT1", "ANKRD13B", "SSH2", "NSRP1", "BLMH", "CPD", "GOSR1", "ATAD5", "TEFM", "NF1", "COPRS", "UTP6", "SUZ12", "LRRC37B", "RHOT1", "C17orf75", "CDK5R1", "MYO1D", "TMEM98", "AP2B1", "TAF15", "ZNHIT3", "MYO19", "GGNBP2", "AATF", "ACACA", "TADA2A", "SYNRG", "DDX52", "MRPL45", "SOCS7", "ARHGAP23", "SRCIN1", "MLLT6", "CISD3", "PCGF2", "CWC25", "RPL23", "LASP1", "RPL19", "FBXL20", "MED1", "CDK12", "STARD3", "MIEN1", "MED24", "THRA", "MSL1", "CASC3", "WIPF2", "RARA", "SMARCE1", "KRT10", "TMEM99", "KRT19", "EIF1", "JUP", "P3H4", "NT5C3B", "ACLY", "CNP", "DNAJC7", "NKIRAS2", "RAB5C", "STAT5B", "STAT3", "ATP6V0A1", "COASY", "MLX", "EZH1", "RAMP2", "VPS25", "AARSD1", "RUNDC1", "RPL27", "VAT1", "RND2", "NBR1", "ARL4D", "DHX8", "DUSP3", "MPP3", "MPP2", "TMEM101", "LSM12", "G6PC3", "HDAC5", "TMUB2", "UBTF", "RUNDC3A", "SLC25A39", "GRN", "FAM171A2", "GPATCH8", "CCDC43", "GJC1", "EFTUD2", "C1QL1", "DCAKD", "NMT1", "MAPT", "KANSL1", "ARL17B", "ARL17A", "NSF", "GOSR2", "NPEPPS", "MRPL10", "SCRN2", "NFE2L1", "CBX1", "CALCOCO2", "UBE2Z", "SNF8", "IGF2BP1", "ZNF652", "PHB", "SPOP", "SLC35B1", "KAT7", "PDK2", "SAMD14", "MRPL27", "LRRC59", "RSAD1", "ANKRD40", "LUC7L3", "TOB1", "SPAG9", "NME1", "NME2", "MBTD1", "UTP18", "TOM1L1", "COX11", "STXBP4", "MMD", "ANKFN1", "DGKE", "COIL", "SCPEP1", "AKAP1", "MSI2", "MRPS23", "VEZF1", "SRSF1", "DYNLL2", "SUPT4H1", "MTMR4", "PPM1E", "TRIM37", "PRR11", "GDPD1", "DHX40", "PTRH2", "VMP1", "HEATR6", "USP32", "APPBP2", "BCAS3", "TBX2", "INTS2", "MED13", "METTL2A", "TANC2", "CYB561", "DCAF7", "TACO1", "MAP3K3", "LIMD2", "CCDC47", "DDX42", "FTSJ3", "POLG2", "DDX5", "CEP95", "SMURF2", "LRRC37A3", "AMZ2P1", "GNA13", "RGS9", "CACNG4", "HELZ", "PITPNC1", "NOL11", "BPTF", "KPNA2", "AMZ2", "ABCA5", "SLC39A11", "SSTR2", "COG1", "FAM104A", "C17orf80", "CDC42EP4", "RPL38", "NAT9", "HID1", "CDR2L", "NT5C", "SUMO2", "MRPS7", "MIF4GD", "GRB2", "TSEN54", "SAP30BP", "GALK1", "UNK", "WBP2", "ACOX1", "SRP68", "EXOC7", "RNF157", "UBALD2", "PRPSAP1", "UBE2O", "CYGB", "PRCD", "SNHG16", "MXRA7", "JMJD6", "METTL23", "SRSF2", "MFSD11", "SEC14L1", "TNRC6C", "TK1", "PGS1", "CYTH1", "TIMP2", "LGALS3BP", "CANT1", "CBX4", "TBC1D16", "GAA", "EIF4A3", "SLC26A11", "RNF213", "ENDOV", "CEP131", "ACTG1", "FAAP100", "NPLOC4", "OXLD1", "CCDC137", "ARL16", "HGS", "MRPL12", "P4HB", "ARHGDIA", "ALYREF", "NPB", "MAFG", "PYCR1", "ASPSCR1", "RAC3", "DCXR", "RFNG", "DUS1L", "FASN", "CCDC57", "NARF", "FOXK2", "WDR45B", "RAB40B", "FN3KRP", "TBCD")
copykat_17q <- subset_data[subset_data$hgnc_symbol %in% genes_17q, ]
rownames(copykat_17q) <- copykat_17q$hgnc_symbol
copykat_17q <- copykat_17q[, -which(names(copykat_17q) == "hgnc_symbol")]


# Computing average CNA values for 1p and 17q genes.
mean_1p <- apply(copykat_1p, 1, mean)  # higher value expectation.
mean_17q <- apply(copykat_17q, 1, mean)  # lower value expectation.

# t-test to compare the two groups
p_value <- t.test(mean_1p, mean_17q)$p.value
p_value

boxplot(mean_1p, mean_17q, names=c("1p Loss", "17q Gain"),
        col=c("blue", "red"), main="Copy Number Changes in Neuroblastoma")
