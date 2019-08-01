#!/usr/bin/env Rscript

### Notes and suggestions are broken up by section
# Data Loading:
# Phenotype association by LM:
# Hub gene network:
# Cell-Type enrichment
#
#
#
#
#
#
#

library(WGCNA)
library(nlme)
library(tidyverse)
library(igraph)
library(pSI)
library(gridExtra)
library(cowplot)

 

# Odds Ratio functions

OR <- function(q,k,m,t) {
    ## 2 x 2 table:
    ##         inTest   !inTest
    ## inRef     q        k
    ## !inRef    m        t
    
    fisher.out <- fisher.test(
        matrix(c(q, k - q, m - q, t - m - k + q), 2, 2), 
        conf.int = TRUE
    )
    OR <- fisher.out$estimate
    pval <- fisher.out$p.value
    upCI <- fisher.out$conf.int[1]
    downCI <- fisher.out$conf.int[2]
    
    output <- c(OR, pval, upCI, downCI)
    names(output) <- c("OR", "Fisher p", "-95%CI", "+95%CI")
    return(output)
}

ORA <- function(testpath,refpath,testbackground,refbackground) {
    q <- length(intersect(testpath,refpath)) ## overlapped pathway size
    k <- length(intersect(refpath,testbackground))  ## input gene set
    m <- length(intersect(testpath,refbackground)) ## input module
    t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)

    empvals <- OR(q,k,m,t)

    tmpnames <- names(empvals)
    empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
    names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
    return(empvals)
}
 
### Data loading

# Network
load("./finallme2/wgcna/ds2.mod50/organoid.lme.network.recut.RData")
geneTree = networks$tree
datExpr=networks$datExpr
merged = networks$merged
modules = merged$colors
MEs = networks$MEs
kMEtable = networks$kMEtable
annot=read.delim("F:/references/gencode.v19.gene.name.txt",header = T,sep = "\t")
annot=annot[match(rownames(datExpr),annot$geneid),]
datMeta=read.delim("./finallme2/1m_metadata_pc_SVA_after_hvariableByLab_removed.txt",header = T,sep = "\t")
datMeta=datMeta[match(colnames(datExpr),rownames(datMeta)),]
datMeta$Group=factor(datMeta$Group,levels = c("CTL","DEL","DUP"))
datMeta$seq=c(1:34)

# Load cell-type expression signatures for enrichment
pSI.zhang = read.csv(file="F:/datasets/psychencode/cellTypes/Zhang_Human_Neuro2015_pSI_GSE21653.csv",row.names=1,stringsAsFactors = F)
pSI.zeisel = read.csv(file="F:/datasets/psychencode/cellTypes/Zeisel_level1_Mouse_Science2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)
pSI.goldmann = read.csv(file="F:/datasets/psychencode/cellTypes/Goldman_levelHybrid_Mouse_NatImmunol2015_pSI.ensg.csv",row.names=1,stringsAsFactors = F)
markers.Habib = read.csv("F:/datasets/psychencode/cellTypes/Habib_NatMeth2017_Human_DroncSeq_clusters.csv",row.names=1)

# protein expression
proExpr=read.delim("./census/joint/union.onProteinLevel.forEigenProtein.txt")

# loop through each module except gray

for(i in 2:ncol(MEs$eigengenes)) {
    me <- MEs$eigengenes[,i]
    moduleNumber <- as.numeric(gsub("ME", "", colnames(MEs$eigengenes)[i]))
    moduleColor <- labels2colors(moduleNumber)
    moduleGenes <- annot$geneid[modules == moduleNumber]

    # Module page init
    pdf(file = paste0("./finallme2/figures/moduleAnnot.ds2mod50/", moduleNumber, moduleColor, ".pdf"), width = 16, height = 11)

    # Phenotype association by LM
    s <- summary(lm(me ~ Group, data = datMeta))$coefficients
    dat2 <- data.frame()
    if (!is.na(s)){
        for(grp in c("DEL", "DUP")) {
        rowID <- paste0("Group", grp)
        dat2 <- rbind(
            dat2,
            data.frame(
                Module = moduleColor, 
                moduleNumber = moduleNumber, 
                Group = grp, 
                beta = s[rowID, "Estimate"], 
                SE = s[rowID, "Std. Error"], 
                t = s[rowID, "t value"], 
                p = s[rowID, "Pr(>|t|)"]
            )
        )
    }
    dat2$fdr = p.adjust(dat2$p, method = "fdr")
    dat2$p.symbol = ""
    dat2$p.symbol[dat2$fdr < 0.05] = ""
    g1 <- ggplot(
        data = dat2, 
        mapping = aes(
            x = Group, y = beta, label = p.symbol
        )
    ) + 
        geom_bar(stat = "identity", fill = moduleColor) + 
        geom_errorbar(
            mapping = aes(
                ymin = beta - SE, ymax = beta + SE
            ), 
            position = position_dodge(width = 0.8), width = 0.25, size = 0.25
        ) +
        geom_text(
            mapping = aes(
                y = beta + 1.3 * sign(beta) * SE 
            )
            color = "red", size = 8
        ) +
        labs(
            x = "", y = "Linear Regression Beta", 
            title = "Module-Trait Association"
        ) +
        theme(
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 15),
            axis.title.y = element_text(size = 20)
        )
    } else {
        df <- data.frame()
        g1 <- ggplot(df) + 
            geom_point() + 
            xlim(0, 10) + 
            ylim(0, 10)
    }
    
    # Module eigengene vs Age trajectory

    dat <- data.frame(Eigengene = me, Age = datMeta$seq, Group = datMeta$Group)
    g2 <- ggplot(
        data = dat,
        mapping = aes(x = Age, y = Eigengene, color = Group)
    ) + 
        geom_point(size = 5)+
        labs(
            title = "Eigengene Expression",
            x = "", y = "Eigengene"
        ) +
        scale_color_manual(
            values = c("black","red","blue")
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 15),
            axis.title.y = element_text(size=20),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 17)
        )
 

    # Gene ontology 
    go <- read.delim(
        paste(
            "./finallme2/wgcna/go/ds2.mod50/bpmf.noiea.overlap3/go.M",
            moduleNumber,
            ".txt",
            sep = ""
        ),
        header = T, sep = "\t"
    )
    go <- go[complete.cases(go), ]
    if(nrow(go) == 0) {
        df <- data.frame()
        g3 <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10)
    } else {
        g3 <- ggplot(
            data = go, 
            mapping = aes(
                x = reorder(term.name, -log10(p.value)), 
                y = -log10(p.value)
            )
        ) + 
            coord_flip() +
            geom_bar(
                stat="identity", fill="royalblue"
            ) +
            labs(title = "Functional Enrichment") +
            geom_hline(yintercept = -log10(0.05), lty = 2, color="red") +
            theme(
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12, face = "bold"),
                axis.title.x = element_text(size = 20)
            )
    }

    # Hub gene network

    hubGenes <- moduleGenes[order(kMEtable[moduleGenes, i], decreasing = T)[1:20]]
    hubGene.symbols <- annot$genename[match(hubGenes, annot$geneid)]
    adjMat <- adjacency(
        t(datExpr[hubGenes, ]), type = "signed", corFnc = "bicor", power = 10
    )
    adjMat[adjMat < quantile(adjMat,0.1)] <- 0
    graph <- graph.adjacency(
        as.matrix(adjMat), mode = "undirected", weighted = T, diag = F
    )
    plotcord <- data.frame(layout_with_fr(graph))
    colnames(plotcord) <- c("X1", "X2")
    edgelist <- get.edgelist(graph, names = F)
    edges <- data.frame(plotcord[edgelist[, 1], ], plotcord[edgelist[, 2], ])
    colnames(edges) <- c("X1", "Y1", "X2", "Y2")
    plotcord <- cbind(plotcord, data.frame(gene = hubGene.symbols))
 
    g4 <- ggplot() + 
        geom_segment(
            data = edges,
            mapping = aes(x = X1, y = Y1, xend = X2, yend = Y2), 
            size = 0.5, colour = "grey"
        ) +
        geom_point(
            data = plotcord,
            mapping = aes(X1, X2), 
            size = 4,
            color = moduleColor
        ) + 
        geom_text(
            data = plotcord,
            mapping = aes(
                x = X1, y = X2 + 0.2, label = gene
            ),
            fontface = "bold",
            size = 6
        ) +
        labs(title = "Hub Gene Network") +
        theme_classic() + 
        theme(
            axis.text = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank()
        )

    # Cell-Type enrichment
    f.zhang = fisher.iteration(pSI.zhang, moduleGenes, p.adjust = T)
    f.zeisel = fisher.iteration(pSI.zeisel, moduleGenes, p.adjust = T)
    f.goldman = fisher.iteration(pSI.goldmann, moduleGenes, p.adjust = T)
    f.habib = data.frame(
        P = matrix(
            NA, 
            nrow = length(unique(na.omit(markers.Habib$Cluster.ID))), 
            ncol=1
        ), 
        row.names = unique(na.omit(markers.Habib$Cluster.ID))
    )

    for(j in 1:nrow(f.habib)) {

        f <- ORA(
            testpath = annot$genename[modules == moduleNumber],
            refpath = markers.Habib$Gene.ID[
                na.omit(markers.Habib$Cluster.ID == rownames(f.habib)[j])
            ],
            testbackground = annot$genename, 
            refbackground = markers.Habib$Gene.ID
        )
        
        if (as.numeric(f[[1]]) > 1) {
            f.habib$P[j] = as.numeric(f[[2]])
        } else {
            f.habib$P[j] = 1
        }
    }

    modCellType = rbind(
        data.frame(Dataset="Zhang", CellType=rownames(f.zhang), log10fdr=-log10(f.zhang[,2])),
        data.frame(Dataset="Zeisel", CellType=rownames(f.zeisel), log10fdr=-log10(f.zeisel[,2])),
        data.frame(Dataset="Goldman", CellType=rownames(f.goldman), log10fdr=-log10(f.goldman[,2])),
        data.frame(Dataset="Habib",  CellType=rownames(f.habib), log10fdr=-log10(f.habib[,1]))
    )
    modCellType$CellType = factor(modCellType$CellType, levels=unique(modCellType$CellType))
    g5 <- ggplot(
        modCellType,
        aes(x=CellType,y=log10fdr, fill=Dataset)
    ) + 
        geom_bar(stat="identity") + 
        coord_flip() + 
        geom_hline(yintercept = -log10(0.05),lty = 2) + 
        ggtitle("Cell-Type Enrichment") +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 10,face = "bold"),
        )

    grid.arrange(
        grobs = list(g1, g2, g3, g4, g5), 
        layout_matrix = rbind(
            c(1,1,2,2,2,5,5),
            c(3,3,3,4,4,4,4)
        ),
        top = paste0(
            "Module ", moduleNumber, " (", moduleColor, ")"
        ), 
        padding = unit(2, "line")
    )

    dev.off()
}
