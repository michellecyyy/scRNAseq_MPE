Cell Clustering
MPE <- Read10X(data.dir = "MPE_filtered_feature_bc_matrix")
MPE <- CreateSeuratObject(counts = MPE.data, project = "MPE", min.cells = 62, min.features = 200)
MPE[["MPE.mt"]] <- PercentageFeatureSet(MPE, pattern = "^MT-")
MPE <- subset(MPE, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
MPE <- NormalizeData(MPE, normalization.method = "LogNormalize", scale.factor = 10000)
MPE <- FindVariableFeatures(MPE, selection.method = "vst", nfeatures = 2000,
                            mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.0125, x.high.cutoff = 3 and y.cutoff = 0.5)
all.genes <- rownames(MPE)
MPE <- ScaleData(MPE, features = all.genes)
MPE <- RunPCA(MPE, features = VariableFeatures(object = MPE))
MPE <- FindNeighbors(MPE, dims = 1:20)
MPE <- FindClusters(MPE, resolution = 0.8)
MPE <- RunTSNE(object = MPE, dims.use = 1:20, do.fast = TRUE)
MPE<- AddMetaData(object = MPE, metadata = Barcode[,2], col.name = 'Patient')
MPE<- AddMetaData(object = MPE, metadata = Barcode[,3], col.name = 'Sample')
MPE<- AddMetaData(object = MPE, metadata = Barcode[,4], col.name = 'Library')
MPE_NoLegend<-TSNEPlot(object = All.cell.tSNE,pt.size=1)+scale_color_npg()+NoLegend()
MPE_Sample_NoLegend<-TSNEPlot(object = All.cell.tSNE,pt.size=1,split.by="Sample")+scale_color_npg()+NoLegend()
MPE.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.25)
markergene<-DotPlot(MPE, features = rownames(Cell.marker),group.by = "order",
                    cols=c("yellow","navy"),dot.scale=10,scale.by="size",
                            scale.min=5)
markergene_1<-FeaturePlot(MPE,reduction="tsne",pt.size=0.5,features = 
                            c("CD3D","CD3E","CD3G","CD2"))
markergene_2<-FeaturePlot(MPE,reduction="tsne",pt.size=0.5,features = 
                            c("MS4A1","CD79A","CD79B","IGHM"))
markergene_3<-FeaturePlot(MPE,reduction="tsne",pt.size=0.5,features = 
                            c("NKG7","GNLY","KLRD1","NCR1"))
markergene_4<-FeaturePlot(MPE,reduction="tsne",pt.size=0.5,features = 
                            c("LYZ", "CD14","FCGR3A","CD163"))
markergene_Sample_1<-FeaturePlot(MPE,reduction="tsne",pt.size=0.1,features = 
                                   c("CD3D","CD3E","CD3G","CD2"),split.by="Sample")
markergene_Sample_2<-FeaturePlot(MPE,reduction="tsne",pt.size=0.1,features = 
                              c("MS4A1","CD79A","CD79B","IGHM"),split.by="Sample")
markergene_Sample_3<-FeaturePlot(MPE,reduction="tsne",pt.size=0.1,features = 
c("NKG7","GNLY","KLRD1","NCR1"),split.by="Sample")
markergene_Sample_4<-FeaturePlot(MPE,reduction="tsne",pt.size=0.1,features = 
c("LYZ","CD14","FCGR3A","CD163"),split.by="Sample")



Differentially Expressed Gene and Pathway Analysis
MPE <- CreateSeuratObject(counts = MPE, min.cells = 62, min.features = 200, 
project = " MPE ")
data <- as(as.matrix(MPE @assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = MPE @meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
MPE _monocle <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
MPE _monocle <- estimateSizeFactors(MPE _monocle)
MPE _monocle <- estimateDispersions(MPE _monocle)
MPE _monocle <- detectGenes(MPE _monocle, min_expr = 0.1)
expressed_genes<-row.names(subset(fData(MPE _monocle),num_cells_expressed >= 10))
start_time <- Sys.time()
function(MPE _monocle, qvalue=qvalue){
    diff_test_res <- differentialGeneTest(
        MPE _monocle,
        fullModelFormulaStr="~cellType",
        cores = 3)
    sig_genes_0.05 <- subset(diff_test_res, qval < 0.05)
    sig_genes_0.01 <- subset(diff_test_res, qval < 0.01)
    print(paste(nrow(sig_genes_0.05), " significantly DE genes (FDR<0.05).", sep=""))
    print(paste(nrow(sig_genes_0.01), " significantly DE genes (FDR<0.01).", sep=""))
    diff_test_res <- subset(diff_test_res, qval< qvalue)
    return(diff_test_res)
} 
end_time <- Sys.time()
end_time - start_time

MPE_go <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)

Trajectory Analysis
lognorm <- t(read.table('MPE_counts.txt', sep=" ", header=TRUE))
anno_table <- read.table('MPE_annotation.txt')
pDat <- data.frame(cell=colnames(lognorm), celltype='undefined', stringsAsFactors=FALSE)
rownames(pDat) <- pDat$cell
pDat[rownames(anno_table), 2] <- as.character(anno_table$celltype)
eset <- Biobase::ExpressionSet(lognorm, phenoData=Biobase::AnnotatedDataFrame(pDat))
dir.create('destiny', showWarnings=FALSE)
dmap <- DiffusionMap(eset)
plot.DiffusionMap(dmap, dims=c(1,2))
plot.DiffusionMap(dmap, dims=c(2,3))
plot.DiffusionMap(dmap, dims=c(1,3))

