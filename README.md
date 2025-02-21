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

包括分离单细胞（或细胞核）、将RNA转化为cDNA、制备测序文库（Illumina）并进行测序。

主要在通量（每次实验捕获多少细胞）、定量类型（全长或基于标签）以及成本方面有所不同。

SMART-seq2 是一种流行的低通量方法，提供全长转录本定量。它非常适合详细研究较小群体的细胞（例如，差异性异构体使用、低表达转录本的表征）。

Chromium 是一种流行的高通量方法，使用唯一分子标识符（UMI）进行转录本定量（从3’端或5’端）。它非常适合研究高度异质的组织并大规模采样大量细胞。

![library](./10X_Partition_2-1024x406.png)

## 3.批次效应

使用不同技术对同一材料进行测序，也可能观察到这些效应（如下图所示），如果没有适当标准化，可能会导致错误的结论。

![batch effect](./batch_effect_protocols.png)

避免批次效应实验设计很重要

![balance](./balanced_design.png)

[Hicks S C, Townes F W, Teng M, et al. Missing data and technical variability in single-cell RNA-sequencing experiments[J]. Biostatistics, 2018, 19(4): 562-578.](https://academic.oup.com/biostatistics/article/19/4/562/4599254?login=false#123896284)

## 4.资源链接

**A useful tool to estimate how many cells to sequence has been developed by the Satija Lab**:https://satijalab.org/howmanycells/

**It is difficult to precisely estimate how much an experiment will cost, although we point to this tool from the Satija Lab as a starting point**:https://satijalab.org/costpercell/
    
**Single-cell best practices**:https://www.sc-best-practices.org/preamble.html

**Analysis of single cell RNA-seq data**:https://www.singlecellcourse.org/index.html

