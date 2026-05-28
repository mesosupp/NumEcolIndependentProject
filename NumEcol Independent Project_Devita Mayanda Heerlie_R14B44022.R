# NumEcol Independent Project, Author: Devita Mayanda Heerlie (李珊珊), Student ID: R14B44022
# --------------------------------------------------
# PART 1: dbRDA with selected variables
# --------------------------------------------------

# Load the library ----
library(vegan)
library(cluster)
library(dendextend)
set.seed(123) # For reproducibility


# Load & orient data ----
otu <- read.delim('https://raw.githubusercontent.com/mesosupp/NumEcolIndependentProject/main/GVMAGs_OTU.txt', row.names = 1)
env <- read.csv('https://raw.githubusercontent.com/mesosupp/NumEcolIndependentProject/main/GVMAGs_env.csv', row.names = 1 )


# Data Preprocessing ----
env[env == "Na"] <- NA
env[] <- lapply(env, function(x) as.numeric(as.character(x)))
env <- env[, colMeans(is.na(env)) <= 0.20] # Drop columns with >20% NA (TOC/Moisture)
env <- na.omit(env)


# Align and Log-Transform OTU Data ----
otu.t <- t(otu)
shared.samples <- intersect(rownames(otu.t), rownames(env))
env <- env[shared.samples, ]
otu.log <- log1p(otu.t[shared.samples, ])
otu.log <- otu.log[, colSums(otu.log) > 0] # Remove zero-abundance OTUs


# Calculate the Distance Matrix ----
bc.dist <- vegdist(otu.log, method = "bray")


# Grouping and Visualization settings ----
OTU.group <- as.factor(gsub("[0-9]+$", "", sub("_.*", "", rownames(env))))
pal <- c("#FF7F00", "#FC8D62", "#4DAF4A", "#A6D854", "#999999", "#984EA3", 
         "#E41A1C", "#990000", "#111BC3", "#377EB8", "#F781BF", "#A91B60", 
         "#000000", "#FFE533")
point.cols <- pal[OTU.group]

## Function to set the margin and shows the legend
apply_legend <- function() {
  par(xpd = TRUE)
  legend(par("usr")[2] + 0.01 * diff(par("usr")[1:2]), par("usr")[4],
         legend = levels(OTU.group), pch = 16, col = pal[1:nlevels(OTU.group)],
         cex = 0.65, bty = "n", title = "OTUs")
  par(xpd = FALSE)}


# Run The Analysis and Plot ----
## Run the PCoA
pcoa <- cmdscale(bc.dist, eig = TRUE)
pcoa.eig <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

## Plot 1: PCoA
par(mar = c(5, 5, 4, 10))
plot(pcoa$points[,1], pcoa$points[,2], pch = 16, col = point.cols, cex = 1.3,
     xlab = paste0("PCoA1 (", pcoa.eig[1], "%)"), ylab = paste0("PCoA2 (", pcoa.eig[2], "%)"),
     main = "PCoA (Bray-Curtis)")
apply_legend()

## Envfit
env.fit <- envfit(pcoa$points, env, permutations = 999)
print(env.fit)
### MAT, Salinity, pH, NO2, TN, TP, SO4, and TOC were highly significant (p<0.05)

## Plot 2: PCoA + Environmental Vectors
par(mar = c(5, 5, 4, 10))
plot(pcoa$points[,1], pcoa$points[,2], pch = 16, col = point.cols, cex = 1.3,
     xlab = paste0("PCoA1 (", pcoa.eig[1], "%)"), ylab = paste0("PCoA2 (", pcoa.eig[2], "%)"),
     main = "PCoA + Environmental Vectors (p < 0.05)")
plot(env.fit, p.max = 0.05, col = "black", lwd = 2, cex = 0.85)
### envfit drop SO4
apply_legend()

## Variable Selection & dbRDA
dbrda.all <- capscale(bc.dist ~ ., data = env)
dbrda.0 <- capscale(bc.dist ~ 1, data = env)
dbrda.sel <- ordiR2step(dbrda.0, scope = formula(dbrda.all), direction = "forward", 
                        R2scope = RsquareAdj(dbrda.all)$adj.r.squared, permutations = 999)
### ordiR2step appropriately dropped NH4 because they did not explain enough unique variation (p>0.10)
sel.vars <- attr(terms(dbrda.sel), "term.labels")
dbrda.final <- capscale(as.formula(paste("bc.dist ~", paste(sel.vars, collapse = " + "))), data = env)


# Plot 3: dbRDA with Selected Environmental Vectors
par(mar = c(5, 5, 4, 10))
ordiplot(dbrda.final, display = c("sites", "bp"), type = "n", main = "dbRDA (Selected Variables)")
points(dbrda.final, display = "sites", pch = 16, col = point.cols, cex = 1.2)
text(dbrda.final, display = "bp", col = "black", cex = 0.9, arrow.mul = 1)
apply_legend()




# --------------------------------------------------
# PART 2: Hierarchical Classification & Validation
# --------------------------------------------------
# Cluster Color ----
clust.pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")

# Compute distances and Ward's clustering tree ----
## Calculate the Ward's distance
bc.dist.ward <- sqrt(bc.dist)
clust.hclust <- as.hclust(agnes(bc.dist.ward, method = "ward"))

## Define cluster memberships for both k=2 and k=4 
groups_k2 <- cutree(clust.hclust, k = 2)
groups_k4 <- cutree(clust.hclust, k = 4)

## Convert baseline tree to a dendrogram object
dend <- as.dendrogram(clust.hclust)


# Shared Helper Functions ----
## Create function to color the cluster
apply_cluster_legend <- function(k = 4) {
  par(xpd = TRUE)
  legend(par("usr")[2] + 0.01 * diff(par("usr")[1:2]), par("usr")[4],
         legend = paste("Cluster", 1:k), pch = 16, col = clust.pal[1:k],
         cex = 0.75, bty = "n", title = "Cluster")
  par(xpd = FALSE)
}


# Ward's Clustering Evaluation (k = 2 vs k = 4) ----
## k = 2 Tree 
dend_k2 <- as.dendrogram(clust.hclust)
labels_colors(dend_k2) <- point.cols[clust.hclust$order] # Insert the OTU color into the labels
dend_k2 <- set(dend_k2, "labels_cex", 0.35)
plot(dend_k2, main = "Ward's Clustering (k = 2)", xlab = "", sub = "") ## Plot k = 2 Tree
cluster.order.k2 <- unique(groups_k2[clust.hclust$order])
rect.dendrogram(dend_k2, k = 2, border = "red", lwd = 2) # Create the cluster

## k = 4 Tree
dend_k4 <- as.dendrogram(clust.hclust)
labels_colors(dend_k4) <- point.cols[clust.hclust$order]
dend_k4 <- set(dend_k4, "labels_cex", 0.35)
plot(dend_k4, main = "Ward's Clustering (k = 4)", xlab = "", sub = "") # Plot k = 4 Tree
cluster.order.k4 <- unique(groups_k4[clust.hclust$order])
rect.dendrogram(dend_k4, k = 4, border = clust.pal[cluster.order.k4], lwd = 2)

# PCoA Spatial Map + Ward Clusters Overlay ----
clust.cols <- clust.pal[groups_k4]
plot(pcoa$points[,1], pcoa$points[,2], pch = 16, col = clust.cols, cex = 1.4,
     xlab = paste0("PCoA1 (", pcoa.eig[1], "%)"), ylab = paste0("PCoA2 (", pcoa.eig[2], "%)"),
     main = "PCoA Space + Ward Clusters (k = 4)")
apply_cluster_legend(k = 4)


# Alternative Linkage Method Validations ----
## Single Linkage Tree
h.single <- hclust(bc.dist, "single")
dend.single <- as.dendrogram(h.single)
labels_colors(dend.single) <- point.cols[h.single$order]
dend.single <- set(dend.single, "labels_cex", 0.35)
plot(dend.single, main = "d. Single Linkage", xlab = "", sub = "")

## Complete Linkage Tree
h.complete <- hclust(bc.dist, "complete")
dend.complete <- as.dendrogram(h.complete)
labels_colors(dend.complete) <- point.cols[h.complete$order]
dend.complete <- set(dend.complete, "labels_cex", 0.35)
plot(dend.complete, main = "e. Complete Linkage", xlab = "", sub = "")

## Average Linkage Tree
h.average <- hclust(bc.dist, "average")
dend.average <- as.dendrogram(h.average)
labels_colors(dend.average) <- point.cols[h.average$order]
dend.average <- set(dend.average, "labels_cex", 0.35)
plot(dend.average, main = "f. Average Linkage", xlab = "", sub = "")
