rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
## 首先载入文章的数据
load(file='../input.Rdata')
counts=a
# using raw counts is the easiest way to process data through Seurat.
counts[1:4,1:4];dim(counts)
library(stringr) 
meta=metadata
head(meta) 
# 下面的基因是文章作者给出的
gs=read.table('top18-genes-in-4-subgroup.txt')[,1]
gs
library(pheatmap)


fivenum(apply(counts,1,function(x) sum(x>0) ))
boxplot(apply(counts,1,function(x) sum(x>0) ))
fivenum(apply(counts,2,function(x) sum(x>0) ))
hist(apply(counts,2,function(x) sum(x>0) ))
# 上面检测了 counts 和 meta 两个变量，后面需要使用

library(Seurat)
# https://satijalab.org/seurat/mca.html
# 构建 Seurat 需要的对象
# 其中 min.cells 和 min.genes 两个参数是经验值
sce <- CreateSeuratObject(counts = counts, 
                          meta.data =meta,
                          min.cells = 5, 
                          min.genes = 2000, 
                          project = "sce")
# 参考：https://github.com/satijalab/seurat/issues/668

sce
?seurat
table(apply(counts,2,function(x) sum(x>0) )>2000)
table(apply(counts,1,function(x) sum(x>0) )>4)
## 可以看到上面的过滤参数是如何发挥作用的
head(sce@meta.data) 
dim(sce@assays$RNA@data)

## 默认使用细胞名字字符串的切割第一列来对细胞进行分组
# 所以在这里只有一个SS2分组信息, 我们就增加一个 group.by 参数
VlnPlot(object = sce, 
        features = c("nFeature_RNA", "nCount_RNA" ), 
        group.by = 'plate',
        ncol = 2)
VlnPlot(object = sce, 
        features = c("nFeature_RNA", "nCount_RNA"), 
        group.by = 'g',
        ncol = 2)
### 同样的发现，普通的层次聚类得到的4组，很明显是检测到的基因数量的差异造成的。

## 可以给sce对象增加一个属性，供QC使用
ercc.genes <- grep(pattern = "^ERCC-", x = rownames(x = sce), value = TRUE)
percent.ercc <- Matrix::colSums(sce[ercc.genes, ]) / Matrix::colSums(sce)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
sce <- AddMetaData(object = sce, metadata = percent.ercc,
                   col.name = "percent.ercc")
VlnPlot(object = sce, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.ercc" ), 
        group.by = 'g',
        ncol = 3)
## 发现一个很有趣的现象，细胞能检测到的基因数量与其含有的ERCC序列反相关。


## 下面是一些可视化函数
VlnPlot(sce,group.by = 'plate',c("Gapdh","Bmp3","Brca1","Brca2","nFeature_RNA"))

FeatureScatter(object = sce,  feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(object = sce,  feature1 = "Brca1", feature2 = "Brca2")

CellScatter(sce,sce@assays$RNA@counts@Dimnames[2][[1]][3], 
         sce@assays$RNA@counts@Dimnames[2][[1]][4],
        )

sce <- NormalizeData(object = sce, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)
sce@assays$RNA@data[1:4,1:4]
pheatmap(as.matrix(sce@assays$RNA@data[gs,]))
# 这里需要理解  dispersion 值
sce <- FindVariableFeatures(object = sce, 
                         mean.function = ExpMean, 
                         dispersion.function = LogVMR )
## 默认值是：x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 1
# 需要认真读说明书：https://satijalab.org/seurat/pbmc3k_tutorial.html

# This function is unchanged from (Macosko et al.), 
# but new methods for variable gene expression identification are coming soon. 
# We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. 


## 根据经验阈值挑选的变化基因个数。
length( sce@assays$RNA@var.features)

# Scaling the data and removing unwanted sources of variation

## 去除一些技术误差，比如"nFeature_RNA", "nCount_RNA"或者ERCC
head(sce@meta.data) 
sce <- ScaleData(object = sce, 
                 vars.to.regress = c("nCount_RNA",'nFeature_RNA',"percent.ercc" ))
# 后面就不需要考虑ERCC序列了。
sce@assays$RNA@scale.data[1:4,1:4]
pheatmap(as.matrix(sce@scale.data[gs,])) 

## 普通的PCA降维
sce <- RunPCA(object = sce, pc.genes = sce@assays$RNA@var.features, 
              do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5)
tmp <- sce@reductions$pca@feature.loadings
head(tmp)
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
PCAPlot(sce, dims = 1:2,group.by = 'plate')
PCAPlot(sce, dims = 1:2,group.by = 'g')
#DimPlot(sce, reduction = "pca",group.by = "g")等价于上式
sce <- ProjectPCA(sce, do.print = FALSE)
DimHeatmap(sce, dims = 1, cells = 100, balanced = TRUE)
DimHeatmap(sce, dims = 1:10, cells = 100, balanced = TRUE)
PCHeatmap(object = sce, dims = 1, cells = 100, balanced = TRUE)
#DimHeatmap()与PCHeatmap()两个用法一致
### 根据参数来调整最后的分组个数
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce,resolution = 0.3)
table(sce@meta.data$RNA_snn_res.0.3)
#Seurat包更新了，所以以下是旧包代码
#sce1 <- FindClusters(object = sce, reduction.type = "pca", 
                    #dims = 1:20, force.recalc = T,
                    #resolution = 0.4, print.output = 0, 
                    #save.SNN = TRUE)
#PrintFindClustersParams(sce1)
#与之前的分组比较
table(sce1@meta.data$RNA_snn_res.0.3,sce@meta.data$g) 
## resolution 是最关键的参数
sce1 <- sce

sce <- RunTSNE(object = sce, dims = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = sce,label =T)
## firstly find marker one by one by change ident.1 
markers_df <- FindMarkers(object = sce, ident.1 = 2, 
                          min.pct = 0.25)
print(x = head(markers_df))
markers_genes =  rownames(head(x = markers_df, n = 5))
# 可视化最后找到的marker基因
VlnPlot(object = sce, features =markers_genes, 
        log = TRUE)
# 首先可视化找到的marker基因
FeaturePlot(object = sce, 
            features =markers_genes, 
            cols = c("grey", "blue"),
            reduction = "tsne")
# 然后可视化文献作者给出的基因
gs <- read.csv(file = "top18-genes-in-4-subgroup.txt",sep = "\t",header = F)
gs <- as.character(gs)
FeaturePlot(object = sce, 
            features =gs[1:18], 
            cols = c("grey", "blue"), 
            reduction = "tsne")
FeaturePlot(object = sce, 
            features =gs[19:36], 
            cols = c("grey", "blue"), 
            reduction = "tsne")
FeaturePlot(object = sce, 
            features =gs[37:54], 
            cols = c("grey", "blue"), 
            reduction = "tsne")
FeaturePlot(object = sce, 
            features =gs[55:72], 
            cols = c("grey", "blue"), 
            reduction = "tsne")


## 对每个类别细胞都找到自己的marker基因
# Then find markers for every cluster compared to all remaining cells, report # only the positive ones
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh = 0.25)
library(dplyr)
sce.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
intersect(top10$gene,gs$V1)
# setting slim.col.label to TRUE will print just the cluster IDS instead of# every cell name
DoHeatmap(object = sce, genes = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
# save(sce.markers,sce,file='sce_seurat.Rdata')
# load(file='sce_seurat.Rdata') 
FeaturePlot(object = sce, 
            features =top10$gene, 
            cols = c("grey", "blue"), 
            reduction = "tsne")

top20 <- sce.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
intersect(top20$gene,gs$V1)
DoHeatmap(object = sce, genes = top20$gene, 
          slim.col.label = TRUE, remove.key = TRUE)

 





