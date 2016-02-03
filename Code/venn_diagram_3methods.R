#### I want to make a venn diagram showing overlaps of significant genes found between DE Control vs. Knockown, DE Zero vs. Nonzero, and nonparametric SAMSeq

library(VennDiagram)

# Get gene lists from each of the three analyses
genes.C.KD <- rownames(read.table("../Ratliff_DE_ControlKD/Results/DEresults.txt"))
genes.zero.nonzero <- read.table("../Ratliff_DE_ZeroNonzero/Results/DEresults_full_zerovnonzero.txt")
genes.zero.nonzero <- rownames(genes.zero.nonzero[genes.zero.nonzero$FDR<0.05,])
genes.SAMSeq <- read.table("Results/siggenes.txt", header=TRUE, colClasses=c(rep("character",3)))$Gene

# Names of categories
category <- c("ControlKD", "ZeroNonzero", "NonParametric")

# Areas of each
area1 <- length(genes.C.KD)
area2 <- length(genes.zero.nonzero)
area3 <- length(genes.SAMSeq)

# Areas of 2-group overlap
n12 <- length(intersect(genes.C.KD, genes.zero.nonzero))
n23 <- length(intersect(genes.zero.nonzero, genes.SAMSeq))
n13 <- length(intersect(genes.C.KD, genes.SAMSeq))

# Area of 3-group overlap
n123 <- length(intersect(intersect(genes.C.KD, genes.zero.nonzero), genes.SAMSeq))

draw.triple.venn(area1=area1, area2=area2, area3=area3, n1=n12, n23=n23, n13=n13, n123=n123, category=category, lty="blank", fill = c("skyblue", "pink1", "mediumorchid"))

