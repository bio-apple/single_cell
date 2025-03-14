# 单细胞测序学习笔记

## 1.单核测序（Single-nucleus RNA-seq）+单细胞RNA测序(Single-Cell RNA-seq)    
已有研究表明，细胞核能够准确反映细胞的所有转录模式在实验设计中选择单细胞还是单核测序主要由组织样本类型驱动，单细胞测序高度依赖新鲜组织，需保持细胞活性。

 **样本状态：是新鲜组织还是冷冻/存档样本？**
        冷冻 → 单核测序；
        新鲜 → 单细胞测序可行。
    
 **组织类型：解离是否会导致特定细胞类型丢失？**
        易损组织（如大脑） → 单核测序；
        易解离（如血液） → 单细胞测序。
    
 **研究目标：关注转录活性还是翻译产物？**
        转录调控 → 单核测序；
        蛋白质合成 → 单细胞测序。

|	| 单核测序 (snRNA-seq)	                 | 单细胞测序 (scRNA-seq)                  |
|----|-----------------------------------|------------------------------------|
|目标	| 分离和富集细胞核                          | 	分离和富集完整单细胞                        |
|样本起始状态| 	新鲜或冷冻组织                          | 	通常需要新鲜组织                          |
|主要步骤| 	组织破碎 → 核分离 → 核富集 → RNA提取 → 文库构建	 | 组织解离 → 单细胞悬液 → 细胞富集 → RNA提取 → 文库构建 |
|下游测序| 	与单细胞测序类似（如10X Genomics）	         | 与单核测序类似（如10X Genomics）             |
|内含子比例	| 高（>50%）| 	低（<10%）                           |


## 2.工作流程

包括分离单细胞（或细胞核）、将RNA转化为cDNA、制备测序文库（Illumina）并进行测序。 主要在通量（每次实验捕获多少细胞）、定量类型（全长或基于标签）以及成本方面有所不同。

SMART-seq2 是一种流行的低通量方法，提供全长转录本定量。它非常适合详细研究较小群体的细胞（例如，差异性异构体使用、低表达转录本的表征）。

Chromium 是一种流行的高通量方法，使用唯一分子标识符（UMI）进行转录本定量（从3’端或5’端）。它非常适合研究高度异质的组织并大规模采样大量细胞。

![library](./10X_Partition_2-1024x406.png)

## 3.bioinformatics

*Roadmap for typical single-cell RNA sequencing data analysis*

![Roadmap for typical single-cell RNA sequencing data analysis](./figure/Roadmap_for_typical_single-cell_RNA_sequencing_data_analysis.jpg)

*Overview of the analysis modules for single-cell RNA sequencing data analysis*

![Overview of the analysis modules for single-cell RNA sequencing data analysis](./figure/Overview_of_the_analysis_modules_for_single-cell_RNA_sequencing_data_analysis.jpg)

[Jovic D, Liang X, Zeng H, et al. Single‐cell RNA sequencing technologies and applications: A brief overview[J]. Clinical and translational medicine, 2022, 12(3): e694.](https://onlinelibrary.wiley.com/doi/full/10.1002/ctm2.694)

### 3-0:[Raw data processing](https://www.sc-best-practices.org/introduction/raw_data_processing.html)

![counts](./figure/overview_raw_data_processing.jpg)

### 3-1：mapping

Raw data processing pipelines such as Cell Ranger、STARsolo, Kallisto, Alevin, and Alevin-fry,Summary of the results for each evaluated section of interest and mapper. Good results are coloured in green, intermediate in yellow, and poor results in red.

![mapping tools](mapping/mapping.jpg)

[Brüning R S, Tombor L, Schulz M H, et al. Comparative analysis of common alignment tools for single-cell RNA sequencing[J]. Gigascience, 2022, 11: giac001.](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giac001/6515741?login=true)

### 3-2:considerations_in_Quality control

![counts](./considerations_in_quality_control/counts.png)

常见的细胞质量控制 (QC) 步骤主要基于以下三个指标：

        每个条形码的 UMI 计数深度 (count depth)
        每个条形码检测到的基因数 (gene count per barcode)
        每个条形码中线粒体基因 UMI 计数的比例 (fraction of mitochondrial counts per barcode)

这些异常条形码可能对应于：

        凋亡细胞 (dying cells)
        细胞膜破损的细胞 (cells with broken membranes)
        多胞 (doublets, 即两个或多个细胞的 mRNA 被误认为是单个细胞的表达数据)

![considerations_in_quality_control](considerations_in_quality_control/Overview_of_unimodal_analysis_steps_for_scRNA-seq.png)

[Kim G D, Lim C, Park J. A practical handbook on single-cell RNA sequencing data quality control and downstream analysis[J]. Molecules and Cells, 2024, 47(9): 100103.](https://www.sciencedirect.com/science/article/pii/S1016847824001286)

如果某个条形码的 UMI 计数低、检测到的基因数少、但线粒体基因占比高，则可能是细胞膜破损导致胞质mRNA泄漏，仅剩线粒体mRNA被保留下来。 

相反，UMI计数异常高、检测到的基因数异常多，可能代表多胞，通常需要设定高UMI计数阈值来过滤掉潜在的多胞。

线粒体基因比例较高的细胞，可能参与呼吸代谢 (respiratory processes)，并不一定是死细胞。

UMI计数较低、基因数较少的细胞，可能属于静息 (quiescent) 细胞群，并非一定是低质量细胞。

UMI 计数较高的细胞，可能只是细胞体积较大，并不一定是多胞。

因此，在设定QC阈值时，应该联合多个QC指标进行筛选，并尽量使用宽松的阈值，以免误删真正的细胞群体。未来，可以使用多变量QC依赖性模型来提高QC筛选的灵敏度。

在单细胞 RNA 测序 (scRNA-seq) 中，空胞指的是： 含有条形码 (barcode) 但没有真正的细胞，仅含有少量环境 RNA (ambient RNA)。 这些 RNA 可能来源于细胞裂解后的游离 RNA，在实验过程中随机进入微滴 (droplet) 或孔板 (well) 中。

### 3-3:Mitochondrial gene content cutoff
*A recent systematic survey of scRNA-seq data suggested that a mitochondrial proportion threshold of 10% is appropriate to distinguish between healthy and low-quality cells in most human tissues, while in mouse tissues, the recommended threshold is 5%.*

[Osorio D, Cai JJ. Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control[J]. Bioinformatics, 2021, 37(7): 963-967.](https://academic.oup.com/bioinformatics/article/37/7/963/5896986?login=false)

### 3-4:remove background noise(remove_ambient_RNA_contamination)

We use our genotype-based estimates to evaluate the performance of three methods (**CellBender, DecontX, SoupX**) that are designed to quantify and remove background noise. We find that CellBender provides the most precise estimates of background noise levels and also yields the highest improvement for marker gene detection.

[Janssen P, Kliesmete Z, Vieth B, et al. The effect of background noise and its removal on the analysis of single-cell expression data[J]. Genome biology, 2023, 24(1): 140.](https://link.springer.com/article/10.1186/s13059-023-02978-x)

### 3-5:Normalization

![normalization](./Normalization/Conceptual_differences_between_variance-stabilizing_transformations.png)

[Ahlmann-Eltze C, Huber W. Comparison of transformations for single-cell RNA-seq data[J]. Nature Methods, 2023, 20(5): 665-672.](https://www.nature.com/articles/s41592-023-01814-1)

### 3-6:highly_variable_gene

        pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

*HVG selection methods can be classified into four categories:*

![HVG_selection_methods](./highly_variable_gene/HVG_selection_methods.png)

*Classification of baseline HVG selection methods.*

![Classification_of_baseline_HVG_selection_methods](./highly_variable_gene/Classification_of_baseline_HVG_selection_methods.png)

[Zhao R, Lu J, Zhou W, et al. A systematic evaluation of highly variable gene selection methods for single-cell RNA-sequencing[J]. bioRxiv, 2024: 2024.08. 25.608519.](https://www.biorxiv.org/content/10.1101/2024.08.25.608519v1.abstract)

### 3-7:batch_effect

![batch effect](./experiments_batch_effect/10-Figure1-1.png)

[Hicks S C, Townes F W, Teng M, et al. Missing data and technical variability in single-cell RNA-sequencing experiments[J]. Biostatistics, 2018, 19(4): 562-578.](https://academic.oup.com/biostatistics/article/19/4/562/4599254?login=false#123896284)

### 3-8:Dimensionality Reduction:PCA+Cluster the cells:UMAP+tSNE

在scRNA-seq数据分析中，我们通过寻找与已知细胞状态或细胞周期阶段相关的细胞身份来描述数据集中的细胞结构。这一过程通常被称为细胞身份注释。
为此，我们将细胞组织成簇，以推断相似细胞的身份。聚类本身是一个常见的无监督机器学习问题。我们可以通过在降维后的表达空间中最小化簇内距离来得出簇。在这种情况下，表达空间决定了细胞在降维表示下的基因表达相似性。例如，这种低维表示可以通过主成分分析（PCA）确定.
PCA：用于前期降维、去噪或输入下游分析（如聚类），不建议直接用于最终可视化。在我们运行PCA之后，需要决定将哪些主成分（PCs）纳入下游分析。我们希望纳入足够多的PCs以保留生物学信号，但又要尽量减少PCs的数量，以避免数据中的噪声干扰。根据每个主成分解释的方差百分比对主成分进行排名
In this example, we can observe an ‘elbow’ around PC 9-10, suggesting that the majority of true signal is captured in the first 10 PCs.

        pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
        ElbowPlot(pbmc)

![PCA](./PCA_tSNE_UMAP/PCA.png)

KNN图由反映数据集中细胞的节点组成。我们首先在PCA降维后的表达空间中为所有细胞计算欧几里得距离矩阵，然后将每个细胞与其K个最相似的细胞连接起来。通常，K的值根据数据集的大小设置为5到100之间。KNN图通过在图中将表达空间中的密集区域表示为密集连接区域，反映了表达数据的底层拓扑结构
![KNN](./PCA_tSNE_UMAP/KNN.jpeg)

        pbmc <- FindNeighbors(pbmc, dims = 1:10)
        pbmc <- FindClusters(pbmc, resolution = 0.5)

UMAP（均匀流形逼近和投影）和t-SNE（t-随机邻域嵌入）是单细胞数据集常用的降维和可视化技术。UMAP最近已成为这类分析的黄金标准，因为它具有更高的计算效率并且能更好地保持全局结构；尽管与t-SNE一样，它在局部距离上的准确性可能更高。
Suggest a resolution of 0.4-1.2 for data sets of ~3,000 cells.

        pbmc <- RunUMAP(pbmc, dims = 1:10)

![UMAP](./PCA_tSNE_UMAP/UMAP.png)

tSNE is slow.tSNE doesn’t scale well to large numbers of cells (10k+)

– UMAP is quite a bit quicker than tSNE

– UMAP can preserve more global structure than tSNE*

– UMAP can run on raw data without PCA preprocessing*

– UMAP can allow new data to be added to an existing projection


 
## 4.资源链接

**A useful tool to estimate how many cells to sequence has been developed by the Satija Lab**:https://satijalab.org/howmanycells/

**It is difficult to precisely estimate how much an experiment will cost, although we point to this tool from the Satija Lab as a starting point**:https://satijalab.org/costpercell/
    
**Single-cell best practices**:https://www.sc-best-practices.org/preamble.html

**Analysis of single cell RNA-seq data**:https://www.singlecellcourse.org/index.html
