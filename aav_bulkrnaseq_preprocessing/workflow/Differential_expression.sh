# 1. Install and load the edgeR package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)

# 2. Load your count data (replace with your actual data)
#   Assume 'counts' is a matrix with genes as rows and samples as columns
#   and 'group' is a vector indicating the experimental condition for each sample
counts <- matrix(rnbinom(100*6, size=10, prob=0.5), ncol=6)
group <- factor(c(1,1,2,2))

# 3. Create a DGEList object
y <- DGEList(counts=counts, group=group)

# 4. Filter low-expressed genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]

# 5. Normalize library sizes (TMM normalization)
y <- calcNormFactors(y)

# 6. Estimate dispersions
y <- estimateDisp(y)

# 7. Perform differential expression analysis (exact test)
et <- exactTest(y)

# 8. Get results and adjust p-values
tt <- topTags(et, n=Inf)
tt$table$padj <- p.adjust(tt$table$PValue, method="BH")

# 9.  Filter for significant results
significant_genes <- tt[tt$table$padj < 0.05, ]