# Cross-lagged panel network (CLPN) analyses.
# Supply the prepared SPSS analysis file by setting CLPN_DATA_FILE.
# If CLPN_DATA_FILE is unset, place analysis_data.sav in the working directory.
# The script writes figures and summary tables to the working directory.

#====Load required packages====
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!require(glmnet)) { install.packages("glmnet"); library(glmnet) } # LASSO regression and cross-validation
if (!require(doParallel)) { install.packages("doParallel"); library(doParallel) } # Parallel cross-validation
if (!require(haven)) { install.packages("haven"); library(haven) } # SPSS data import
if (!require(dplyr)) { install.packages("dplyr"); library(dplyr) } # Data manipulation
if (!require(tidyr)) { install.packages("tidyr"); library(tidyr) } # Reshaping results
if (!require(qgraph)) { install.packages("qgraph"); library(qgraph) } # Network plots and centrality
if (!require(bootnet)) { install.packages("bootnet"); library(bootnet) } # Network estimation and bootstrap analyses
if (!require(lavaan)) { install.packages("lavaan"); library(lavaan) } # Model-based node predictability
if (!require(ggplot2)) { install.packages("ggplot2"); library(ggplot2) } # Plots and themes
if (!require(patchwork)) { install.packages("patchwork"); library(patchwork) } # Multi-panel figures
if (!require(scales)) { install.packages("scales"); library(scales) } # Numeric axis labels via scales::number_format()
if (!require(writexl)) { install.packages("writexl"); library(writexl) } # Excel export
if (!require(networktools)) { install.packages("networktools"); library(networktools) } # Bridge centrality

#====Import the analysis data====
data_file <- Sys.getenv("CLPN_DATA_FILE", unset = "analysis_data.sav")
if (!file.exists(data_file)) {
  stop("Set CLPN_DATA_FILE to the path of the prepared SPSS analysis file.")
}

# Read only the analysis columns. The first two slots can be blank placeholders
# and are never imported. The order below follows the prepared analysis file.
# In each wave: 2 smartphone usage items, 20 depression items, 6 positive
# suicidal ideation items, and 8 negative suicidal ideation items.
# The final 3 selected columns contain the numeric analysis covariates.
data1 <- read_sav(data_file, col_select = 3:113)
LBCPHT1 <- data1[, 1:36]
LBCPHT2 <- data1[, 37:72]
LBCPHT3 <- data1[, 73:108]
LBCPHCov <- data1[, 109:111]

# Combine adjacent waves and append the covariates.
LBCPHT1T2 <- as.data.frame(bind_cols(LBCPHT1,LBCPHT2))
LBCPHT1T2.cov <- as.data.frame(bind_cols(LBCPHT1T2, LBCPHCov))

LBCPHT2T3 <- as.data.frame(bind_cols(LBCPHT2,LBCPHT3))
LBCPHT2T3.cov <- as.data.frame(bind_cols(LBCPHT2T3, LBCPHCov))

# ====CLPN estimation with covariates====
# Set number of covariates
numCovar = 3
# Set number of symptoms
k = 36
k1 = 2  # Number of Variables in Construct 1
k2 = 20 # Number of Variables in Construct 2
k3 = 6  # Number of Variables in Construct 3
k4 = 8  # Number of Variables in Construct 3

# Estimate covariate-adjusted coefficients for each wave pair
## Set up empty matrix of coefficients
adjMatCov.T1T2 <- matrix(0, (numCovar+k), (numCovar+k))
adjMatCov.T2T3 <- matrix(0, (numCovar+k), (numCovar+k))

# Set the number of workers for cross-validation
available_cores <- parallel::detectCores()
num_cores <- if (is.na(available_cores)) 1L else max(1L, floor(available_cores / 2 - 1))
registerDoParallel(cores = num_cores)

# Fit one LASSO regression per later-wave node
for (i in 1:k) {
  set.seed(100)
  lassoregCov.T1T2 <- cv.glmnet(as.matrix(LBCPHT1T2.cov[,c(1:k,(k*2+1):(k*2+numCovar))]),  # Use all earlier-wave nodes and the covariates as predictors
                                LBCPHT1T2.cov[,(k+i)],  # Use the ith node at T2 as the outcome
                                nfolds = 10,  # Use 10-fold cross-validation
                                family = "gaussian",
                                alpha = 1,
                                standardize = TRUE,  # Standardize predictors internally for model fitting
                                parallel = TRUE)  # Use the registered parallel backend
  lambdaCov.T1T2 <- lassoregCov.T1T2$lambda.min  # Select the lambda with the lowest cross-validated error
  adjMatCov.T1T2[(1:(k+numCovar)), i] <- coef(lassoregCov.T1T2,  # Store the regression coefficients in column i
                                              s = lambdaCov.T1T2,
                                              exact = FALSE)[2:(numCovar+k+1)]  # Omit the intercept
}

for (i in 1:k) {
  set.seed(100)
  lassoregCov.T2T3 <- cv.glmnet(as.matrix(LBCPHT2T3.cov[,c(1:k,(k*2+1):(k*2+numCovar))]),  # Use all earlier-wave nodes and the covariates as predictors
                                LBCPHT2T3.cov[,(k+i)],  # Use the ith node at T3 as the outcome
                                nfolds = 10,  # Use 10-fold cross-validation
                                family = "gaussian",
                                alpha = 1,
                                standardize = TRUE,  # Standardize predictors internally for model fitting
                                parallel = TRUE)  # Use the registered parallel backend
  lambdaCov.T2T3 <- lassoregCov.T2T3$lambda.min  # Select the lambda with the lowest cross-validated error
  adjMatCov.T2T3[(1:(k+numCovar)), i] <- coef(lassoregCov.T2T3,  # Store the regression coefficients in column i
                                              s = lambdaCov.T2T3,
                                              exact = FALSE)[2:(numCovar+k+1)]  # Omit the intercept
}
# Release the parallel backend after fitting
stopImplicitCluster()
foreach::registerDoSEQ()

# Extract the node-by-node matrix from the covariate-adjusted fit
adjMat.T1T2 <- adjMatCov.T1T2[1:k, 1:k]
adjMat.T2T3 <- adjMatCov.T2T3[1:k, 1:k]

# ====Autoregressive edge analysis====
# Exclude autoregressive edges
adjMat.noAR.T1T2 <- adjMat.T1T2
diag(adjMat.noAR.T1T2) <- 0  # Set the diagonal to zero

adjMat.noAR.T2T3 <- adjMat.T2T3
diag(adjMat.noAR.T2T3) <- 0  # Set the diagonal to zero

# Count nonzero cross-lagged edges
sum(adjMat.noAR.T1T2!=0)
sum(adjMat.noAR.T2T3!=0)

# Create plot titles
LabwArrow.T1T2 <- expression(bold("T1" %->% "T2 Network"))
LabwArrow.T2T3 <- expression(bold("T2" %->% "T3 Network"))

# Extract autoregressive coefficients
ARedges.T1T2 <- diag(adjMat.T1T2)
ARedges.T2T3 <- diag(adjMat.T2T3)

# Define node labels for the autoregressive analysis
labelnames = c("H01","H02",
               "C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20",
               "SP01","SP02","SP03","SP04","SP05","SP06",
               "SN01","SN02","SN03","SN04","SN05","SN06","SN07","SN08")
shortnames = labelnames
longnames<-c(
  # Smartphone usage time
  "Time of smartphone usage per day in working days","Time of smartphone usage per day in weekends",
  # Depressive symptoms
  "Some things that usually don't bother me make me feel annoyed","I don't feel like eating; my appetite is poor",
  "I feel that I can't shake off the blues even with help from my family or friends","I feel I am not good as same as other people",
  "I have trouble keeping my mind on what I'm doing","I feel depressed","I feel that everything I do is an effort",
  "I feel depressed about the future","I think my life has been a failure","I feel fear","My sleep is restless","I feel unhappy",
  "I eat less than usual","I feel lonely","People are unfriendly","I do not enjoy life","I have crying spells",
  "I feel sad","I feel that people dislike me","I walk slowly",
  # Positive suicidal ideation
  "I have been very smooth in my studies and feel very happy","I feel that most of the situations in life are under my control",
  "I am very satisfied with life and feel hopeful for the future","Believe that you can effectively handle problem",
  "Feeling that life is worth living", "I feel confident in the future",
  # Negative suicidal ideation
  "Feeling hopeless in the future and having suicidal thoughts",
  "Feeling unhappy in their relationship with someone and wanting to commit suicide",
  "Unable to accomplish important things in life, wanting to be reborn",
  "Being helpless about one's own problems and wanting to commit suicide",
  "Feeling unable to meet the requirements of others,resulting in suicidal thoughts",
  "The problem is so difficult to solve, and suicide seems to be the only option",
  "Feeling very lonely and sad, wanting to end the pain and have a life", "Feeling frustrated in life and wanting to have a life of suicide")
# Assign nodes to constructs
groups <- c(rep("Smartphone Usage Time", k1), rep("CES-D", k2),
            rep("Positive Suicide Ideation",k3), rep("Negative Suicide Ideation", k4))

# Assemble the autoregressive results
ARedges.T1T2.Result <- data.frame(Node=rep(factor(shortnames, levels=rev(shortnames)), 1),
                                  Constructs=rep(groups, 1),
                                  Network=c(rep("T1 -> T2 Network", k)),
                                  ARedge=c(ARedges.T1T2))

ARedges.T2T3.Result <- data.frame(Node=rep(factor(shortnames, levels=rev(shortnames)), 1),
                                  Constructs=rep(groups, 1),
                                  Network=c(rep("T2 -> T3 Network", k)),
                                  ARedge=c(ARedges.T2T3))
ARedges.Result <- rbind(ARedges.T1T2.Result, ARedges.T2T3.Result)

ARedges.Result$Network <- factor(ARedges.Result$Network,  # Order the intervals chronologically
                                 levels = c("T1 -> T2 Network",
                                            "T2 -> T3 Network"))
# Plot autoregressive coefficients
ARedges.p <- ggplot(ARedges.Result, aes(x=Node, y=ARedge, group=Network)) +
  geom_path(aes(lty=Network),size = 0.7) +  # Connect points along each interval
  geom_point(aes(fill = Constructs), shape=21, size=3) +  # Fill points by construct
  coord_flip() +  # Flip the axes
  ylab("") + xlab("Nodes") +  # Set axis labels
  labs(title="Autoregressive edges") +
  scale_linetype_manual(values = c("solid", "82")) +
  scale_fill_manual(values=c("#a8e6cf","#74ADD1", "#FC8D62", "#FFFF33")) +  # Assign a color to each construct
  scale_color_manual(values=c("#a8e6cf","#74ADD1", "#FC8D62", "#FFFF33")) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0,2)) + # to make bars start at 0
  guides(lty=guide_legend(order=1, keywidth=unit(2,"cm")),  # Set the legend order
         fill=guide_legend(order=2),
         col=guide_legend(order = 3, title = "Cog")) +
  theme_bw() +  # Use a black-and-white theme
  theme(text=element_text(color="black"),  # Set text appearance
        plot.title=element_text(size=14, color="black", hjust=.5, face="bold"),
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black", face="bold"),
        legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=12, color="black"))
print(ARedges.p)
ggsave("LB_T1T3_CPH_ARedges.tiff",plot = ARedges.p,
       width = 10,height = 10,dpi = 300,compression = "lzw")
# Export autoregressive coefficients
AR_MasterTable <- data.frame(
  Node = shortnames,  # Node labels in matrix order
  Construct = groups,  # Construct membership
  AR_T1T2 = ARedges.T1T2,  # T1 to T2 autoregressive coefficients
  AR_T2T3 = ARedges.T2T3  # T2 to T3 autoregressive coefficients
)

write_xlsx(AR_MasterTable, "LB_CPH_Autoregressive_Coefficients.xlsx")

# Optional CSV export
# write.csv(AR_MasterTable, "LB_CPH_Autoregressive_Coefficients.csv", row.names = FALSE)
# Calculate mean signed autoregressive and cross-lagged coefficients
#T1 to T2
print("T1 to T2 Mean Autoregressive")
ARedge_mean_T1T2 <- mean(ARedges.T1T2)  # Mean signed autoregressive coefficient

print("T1 to T2 Mean Regressive without AR")
Edge_mean_T1T2 <- mean(c(adjMat.T1T2[lower.tri(adjMat.T1T2)],  # Extract the lower triangle
                         adjMat.T1T2[upper.tri(adjMat.T1T2)]))  # Extract the upper triangle

#T2 to T3
print("T2 to T3 Mean Autoregressive")
ARedge_mean_T2T3 <- mean(ARedges.T2T3)  # Mean signed autoregressive coefficient

print("T2 to T3 Mean Regressive without AR")
Edge_mean_T2T3 <- mean(c(adjMat.T2T3[lower.tri(adjMat.T2T3)],  # Extract the lower triangle
                         adjMat.T2T3[upper.tri(adjMat.T2T3)]))  # Extract the upper triangle

summary_table <- data.frame(
  Time_Interval = c("T1 -> T2", "T2 -> T3"),
  # Mean autoregressive coefficients
  Mean_Autoregressive = c(ARedge_mean_T1T2, ARedge_mean_T2T3),
  # Mean cross-lagged coefficients
  Mean_CrossLagged = c(Edge_mean_T1T2, Edge_mean_T2T3)
)
write.csv(summary_table, "Network_Edge_Mean.csv", row.names = F)

# ====CLPN visualization====
# Optional plots with covariate nodes
# Define node groups and a separate covariate group
groups_cov <- list(
  Main_Variables1 = 1:k1,  # Smartphone usage time
  Main_Variables2 = (k1+1):(k1+k2),  # Depressive symptoms
  Main_Variables3 = (k1+k2+1):(k1+k2+k3),  # Positive suicidal ideation
  Main_Variables4 = (k1+k2+k3+1):(k1+k2+k3+k4),  # Negative suicidal ideation
  Covariates = (k+1):(k+numCovar)  # Covariate nodes
)

# Define plot labels
main_labels1 <- paste0("TSU", 1:k1)  # Smartphone usage time labels
main_labels2 <- paste0("CES-D", 1:k2)  # CES-D labels
main_labels3 <- paste0("PSI", 1:k3)  # Positive suicidal ideation labels
main_labels4 <- paste0("NSI", 1:k4)  # Negative suicidal ideation labels
covar_labels <- paste0("Cov", 1:numCovar)  # Covariate labels
labels_cov <- c(main_labels1, main_labels2, main_labels3, main_labels4, covar_labels)

# Plot the network
qgraph(adjMatCov.T1T2,
       groups = groups_cov,
       labels = labels_cov,
       legend = TRUE,
       layout = "spring",
       theme = "colorblind",
       vsize = c(rep(4, k), rep(3, numCovar)),  # Draw covariate nodes slightly smaller
       shape = c(rep("circle", k), rep("square", numCovar)),  # Use distinct shapes for covariates
       title = "Cross-lagged Network with Covariates")

# Hide autoregressive edges
adjMatCov.noAR.T1T2 <- adjMatCov.T1T2
diag(adjMatCov.noAR.T1T2) <- 0

qgraph(adjMatCov.noAR.T1T2,
       groups = groups_cov,
       labels = labels_cov,
       legend = TRUE,
       layout = "spring",
       theme = "colorblind",
       vsize = c(rep(4, k), rep(3, numCovar)),
       shape = c(rep("circle", k), rep("square", numCovar)),
       title = "Cross-lagged Network with Covariates noAR")

# Plots without covariate nodes
# Calculate maximum edge weights and graph layouts
maxEdge.T1T2 <- max(c(adjMat.T1T2))
maxEdge.noAR.T1T2 <- max(c(adjMat.noAR.T1T2))  # Calculate the maximum edge weight
maxlayoutAR.noAR.T1T2 <- averageLayout(adjMat.noAR.T1T2)  # Calculate the graph layout

maxEdge.T2T3 <- max(c(adjMat.T2T3))
maxEdge.noAR.T2T3 <- max(c(adjMat.noAR.T2T3))  # Calculate the maximum edge weight
maxlayoutAR.noAR.T2T3 <- averageLayout(adjMat.noAR.T2T3)  # Calculate the graph layout

# Show autoregressive edges
legend_text <- paste(labelnames, longnames, sep = ": ")
# Define node colors
colors_palette <- c("#FFFF33", "#a8e6cf", "#FC8D62", "#74ADD1")
node_colors <- c(
  rep(colors_palette[1], k1),
  rep(colors_palette[2], k2),
  rep(colors_palette[3], k3),
  rep(colors_palette[4], k4)
)

# Draw the plot
tiff("LB_T1T3_CPH_Network.tiff", width = 24, height = 12, units = "in",res = 300, compression = "lzw", family = "serif")
lay_mat <- matrix(c(1, 2,
                    3, 3), nrow = 2, byrow = TRUE)
layout(lay_mat, heights = c(0.75, 0.25))
par(family = "serif", font.main = 2)

adjMat.T1T2.q <- qgraph(adjMat.T1T2, groups = groups, labels = labelnames,
                        legend = F,vsize=6,cut=0,
                        border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                        nodeNames=longnames,legend.cex=.3,GLratio = 2,
                        layout = "spring",
                        title = LabwArrow.T1T2, title.cex = 2)

adjMat.T2T3.q <- qgraph(adjMat.T2T3, groups = groups, labels = labelnames,
                        legend = F,vsize=6,cut=0,
                        border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                        nodeNames=longnames,legend.cex=.3,GLratio = 2,
                        layout = adjMat.T1T2.q$layout,
                        title = LabwArrow.T2T3, title.cex = 2)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Node legend title
text(x = 1, y = 1.35,  # Set the vertical position of the legend title
     labels = "Nodes",
     font = 2,  # Use bold text
     cex = 2)  # Set the text size
legend("center",
       legend = legend_text,  # Node legend text
       
       # Node legend markers
       pch = 21,  # Use filled circles
       pt.bg = node_colors,  # Match the fill to the node color
       pt.cex = 1.5,  # Set the marker size
       col = "black",  # Set the marker border color
       # -----------------------
       
       bty = "n",  # Omit the legend box
       cex = 1.2,  # Set the legend text size
       ncol = 4,  # Arrange the legend in columns
       x.intersp = 0.8,  # Set the spacing between symbols and labels
       text.width = 0.2 )
legend("bottom",
       legend = c("Smartphone Usage Time", "CES-D", "Positive Suicide Ideation", "Negative Suicide Ideation"),
       pt.bg = colors_palette,  # Use one color per construct
       pch = 22,  # Use square markers
       col = "black", bty = "n", horiz = TRUE, cex = 1.2, pt.cex = 2)
dev.off()

# Hide autoregressive edges
tiff("LB_T1T3_CPH_Network_noAR.tiff", width = 24, height = 12, units = "in",res = 300, compression = "lzw", family = "serif")
lay_mat <- matrix(c(1, 2,
                    3, 3), nrow = 2, byrow = TRUE)
layout(lay_mat, heights = c(0.75, 0.25))
par(family = "serif", font.main = 2)

adjMat.noAR.T1T2.q <- qgraph(adjMat.noAR.T1T2, groups = groups, labels = labelnames,
                             legend = F,vsize=6,cut=0, minimum = quantile(abs(adjMat.noAR.T1T2[adjMat.noAR.T1T2 != 0]), 0.80),
                             border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                             nodeNames=longnames,legend.cex=.3,GLratio = 2,
                             layout = adjMat.T1T2.q$layout,
                             title = LabwArrow.T1T2, title.cex = 2)

adjMat.noAR.T2T3.q <- qgraph(adjMat.noAR.T2T3, groups = groups, labels = labelnames,
                             legend = F,vsize=6,cut=0, minimum = quantile(abs(adjMat.noAR.T2T3[adjMat.noAR.T2T3 != 0]), 0.80),
                             border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                             nodeNames=longnames,legend.cex=.3,GLratio = 2,
                             layout = adjMat.T1T2.q$layout,
                             title = LabwArrow.T2T3, title.cex = 2)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Node legend title
text(x = 1, y = 1.35,  # Set the vertical position of the legend title
     labels = "Nodes",
     font = 2,  # Use bold text
     cex = 2)  # Set the text size
legend("center",
       legend = legend_text,  # Node legend text
       pch = 21,  # Use filled circles
       pt.bg = node_colors,  # Match the fill to the node color
       pt.cex = 1.5,  # Set the marker size
       col = "black",  # Set the marker border color
       bty = "n",  # Omit the legend box
       cex = 1.2,  # Set the legend text size
       ncol = 4,  # Arrange the legend in columns
       x.intersp = 0.8,  # Set the spacing between symbols and labels
       text.width = 0.2 )
# Construct legend
legend("bottom",
       legend = c("Smartphone Usage Time", "CES-D", "Positive Suicide Ideation", "Negative Suicide Ideation"),
       pt.bg = colors_palette,  # Use one color per construct
       pch = 22,  # Use square markers
       col = "black", bty = "n", horiz = TRUE, cex = 1.2, pt.cex = 2)
dev.off()

# ====CLPN estimator for node and edge analyses====
# Define the custom estimator
CLPN.net.cov <- function(data) {
  library(glmnet)
  numCovar <- 3
  k <- 36
  ## create empty adjacency matrix
  adjMatCovCLPN <- matrix(0, k+numCovar, k+numCovar)  # Initialize the coefficient matrix
  ## run CLPN loop to do series of nodewise regularized regressions
  for (i in 1:k) {
    #set.seed(100)
    lassoreg <- cv.glmnet(x=as.matrix(data[,c(1:k,(k*2+1):(k*2+numCovar))]),
                          y=data[,(k+i)], nfolds=10,
                          family="gaussian", alpha=1, standardize=TRUE, parallel=F)
    lambda <- lassoreg$lambda.min
    ## paste coefficients into adjacency matrix
    adjMatCovCLPN[1:(k+numCovar),i] <- coef(lassoreg, s=lambda, exact=FALSE)[2:(k+numCovar+1)]
  }
  ## remove covariates from adjacency matrix
  adjMatCLPN <- adjMatCovCLPN[1:k, 1:k]
  return(adjMatCLPN)
}

# ====Fit the CLPNs for bootstrap analyses====

set.seed(100)
adjMat.T1T2.e <- estimateNetwork(LBCPHT1T2.cov, fun=CLPN.net.cov, labels=labelnames, directed=T)
sink("LBCPHT1T2.txt")
summary(adjMat.T1T2.e)
sink()

plot(adjMat.T1T2.e)

set.seed(100)
adjMat.T2T3.e <- estimateNetwork(LBCPHT2T3.cov, fun=CLPN.net.cov, labels=labelnames, directed=T)
sink("LBCPHT2T3.txt")
summary(adjMat.T2T3.e)
sink()

plot(adjMat.T2T3.e)


#====bootnet====
set.seed(100)
adjMat.T1T2.eboot <- bootnet(adjMat.T1T2.e, type="nonparametric",
                             nBoots=1000, communities=groups,
                             directed=T,
                             statistics=c("edge","strength","outStrength","inStrength","expectedInfluence",
                                          "outExpectedInfluence","inExpectedInfluence","betweenness",
                                          "closeness"),  # Include both incoming and outgoing expected influence
                             nCores=num_cores+2)

write.csv(adjMat.T1T2.eboot$sampleTable, "adjMat.T1T2.eboot.sampleTable.csv")
sum_adjT1T2.edgeboot<- summary(adjMat.T1T2.eboot, statistics=c("edge","strength",
                                                               "outStrength","inStrength",
                                                               "expectedInfluence",
                                                               "outExpectedInfluence",
                                                               "inExpectedInfluence",
                                                               "betweenness",
                                                               "closeness"))
write.csv(sum_adjT1T2.edgeboot,"sum_adjT1T2.edgeboot.csv", row.names = F)

gc()

set.seed(100)
adjMat.T2T3.eboot <- bootnet(adjMat.T2T3.e, type="nonparametric",
                             nBoots=1000, communities=groups,
                             directed=T,
                             statistics=c("edge","strength","outStrength","inStrength","expectedInfluence",
                                          "outExpectedInfluence","inExpectedInfluence","betweenness",
                                          "closeness"),
                             nCores=num_cores+2)

write.csv(adjMat.T2T3.eboot$sampleTable, "adjMat.T2T3.eboot.sampleTable.csv")
sum_adjT2T3.edgeboot<- summary(adjMat.T2T3.eboot,statistics=c("edge","strength",
                                                              "outStrength","inStrength",
                                                              "expectedInfluence",
                                                              "outExpectedInfluence",
                                                              "inExpectedInfluence",
                                                              "betweenness",
                                                              "closeness"))
write.csv(sum_adjT2T3.edgeboot,"sum_adjT2T3.edgeboot.csv", row.names = F)

gc()

# ====Case-dropping bootstrap for node centrality====
set.seed(100)
adjMatT1T2.nboot <- bootnet(adjMat.T1T2.e, type="case", nCores = num_cores+2,
                            nBoots=1000,
                            directed=T, statistics=c("OutStrength", "InStrength", "OutExpectedInfluence", "InExpectedInfluence"))
write.csv(adjMatT1T2.nboot$sampleTable, "adjMat.T1T2.nboot.sampleTable.csv")
gc()

set.seed(100)
adjMatT2T3.nboot <- bootnet(adjMat.T2T3.e, type="case", nCores = num_cores+2,
                            nBoots=1000,
                            directed=T, statistics=c("OutStrength", "InStrength", "OutExpectedInfluence", "InExpectedInfluence"))
write.csv(adjMatT2T3.nboot$sampleTable, "adjMat.T2T3.nboot.sampleTable.csv")
gc()

# ====Bootstrap intervals for node centrality====
tiff("adjT1T2.inEIeboot.tiff", width = 10, height = 20, units = "in",res = 300,compression = "lzw", family = "serif")
plot(adjMat.T1T2.eboot, statistics = c("inExpectedInfluence"),
     plot = "area", order="sample", labels=T, legend=F, meanColor=NA) +
  xlab("inExpectedInfluence") +  # Label the x-axis
  ggtitle(LabwArrow.T1T2) +  # Set the plot title
  theme(plot.title=element_text(hjust=0.5, size=18, face="bold"),  # Style the plot title
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=18))
dev.off()

tiff("adjT1T2.outEIboot.tiff", width = 10, height = 20, units = "in",res = 300,compression = "lzw", family = "serif")
plot(adjMat.T1T2.eboot, statistics = c("outExpectedInfluence"),
     plot = "area", order="sample", labels=T, legend=F, meanColor=NA) +
  xlab("outExpectedInfluence") +  # Label the x-axis
  ggtitle(LabwArrow.T1T2) +  # Set the plot title
  theme(plot.title=element_text(hjust=0.5, size=18, face="bold"),  # Style the plot title
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=18))
dev.off()

tiff("adjT2T3.inEIeboot.tiff", width = 10, height = 20, units = "in",res = 300,compression = "lzw", family = "serif")
plot(adjMat.T2T3.eboot, statistics = c("inExpectedInfluence"),
     plot = "area", order="sample", labels=T, legend=F, meanColor=NA) +
  xlab("inExpectedInfluence") +  # Label the x-axis
  ggtitle(LabwArrow.T2T3) +  # Set the plot title
  theme(plot.title=element_text(hjust=0.5, size=18, face="bold"),  # Style the plot title
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=18))
dev.off()

tiff("adjT2T3.outEIboot.tiff", width = 10, height = 20, units = "in",res = 300,compression = "lzw", family = "serif")
plot(adjMat.T2T3.eboot, statistics = c("outExpectedInfluence"),
     plot = "area", order="sample", labels=T, legend=F, meanColor=NA) +
  xlab("outExpectedInfluence") +  # Label the x-axis
  ggtitle(LabwArrow.T2T3) +  # Set the plot title
  theme(plot.title=element_text(hjust=0.5, size=18, face="bold"),  # Style the plot title
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=18))
dev.off()

# ====Bootstrap intervals for edge weights====
tiff("adjT1T2.edgeboot.tiff", width = 20, height = 100, units = "in",res = 300,compression = "lzw", family = "serif")
plot(adjMat.T1T2.eboot, order="sample", labels=T, legend=F, meanColor=NA) +
  xlab("Edge Weight - T1 to T2") +  # Label the x-axis
  ggtitle(LabwArrow.T1T2) +  # Set the plot title
  theme(plot.title=element_text(hjust=0.5, size=18, face="bold"),  # Style the plot title
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=18))
dev.off()

tiff("adjT2T3.edgeboot.tiff", width = 20, height = 100, units = "in",res = 300,compression = "lzw", family = "serif")
plot(adjMat.T2T3.eboot, order="sample", labels=T, legend=F, meanColor=NA) +
  xlab("Edge Weight - T2 to T3") +  # Label the x-axis
  ggtitle(LabwArrow.T2T3) +  # Set the plot title
  theme(plot.title=element_text(hjust=0.5, size=18, face="bold"),  # Style the plot title
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=18))
dev.off()

# ====Bootstrap difference tests for edge weights====
tiff("adjT1T2.edgedif.tiff", width = 18, height = 15, units = "in",res = 300, compression = "lzw", family = "serif")
plot(adjMat.T1T2.eboot,
     "edge", plot="difference",
     onlyNonZero=T,  # Show only nonzero edges
     order="sample") +
  ggtitle(LabwArrow.T1T2) +
  theme(plot.title=element_text(hjust=0.5, size=10, face="bold"),
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=2.5),
        axis.text.y=element_text(size=2.5))
dev.off()
gc()

tiff("adjT2T3.edgedif.tiff", width = 18, height = 15, units = "in",res = 300, compression = "lzw", family = "serif")
plot(adjMat.T2T3.eboot,
     "edge", plot="difference",
     onlyNonZero=T,  # Show only nonzero edges
     order="sample") +
  ggtitle(LabwArrow.T2T3) +
  theme(plot.title=element_text(hjust=0.5, size=10, face="bold"),
        strip.text=element_blank(), # remove "edge" header
        axis.text.x=element_text(size=2.5),
        axis.text.y=element_text(size=2.5))
dev.off()
gc()

# ====Node centrality measures====
adjMat.T1T2.nboot <- adjMatT1T2.nboot
adjMat.T2T3.nboot <- adjMatT2T3.nboot
Cen.T1T2 <- centrality_auto(adjMat.T1T2.q)
write.csv(Cen.T1T2$node.centrality,"LB_CPH_T1T2_Centality.csv")

Cen.T2T3 <- centrality_auto(adjMat.T2T3.q)
write.csv(Cen.T2T3$node.centrality,"LB_CPH_T2T3_Centality.csv")

# Draw the plot
cen.title.T1T2 <- expression("T1" %->% "T2 Network Centrality")
cen.title.T2T3 <- expression("T2" %->% "T3 Network Centrality")

#T1 to T2
tiff(filename = "cenALL.T1T2.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
centralityPlot(adjMat.T1T2.q,include = c("all"),scale = "z-scores") +
  ggtitle(cen.title.T1T2) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  )

dev.off()


cen.outEI.T1T2<- centralityPlot(adjMat.T1T2.q,include = c("OutExpectedInfluence"),scale = "z-scores",orderBy = "OutExpectedInfluence") +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  )
cen.inEI.T1T2<- centralityPlot(adjMat.T1T2.q,include = c("InExpectedInfluence"),scale = "z-scores",orderBy = "InExpectedInfluence") +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  )
cen.EI.T1T2 <- cen.outEI.T1T2 + cen.inEI.T1T2 +
  plot_annotation(title = cen.title.T1T2,
                  theme = theme(
                    plot.title = element_text(family = "serif", face = "bold", size = 16, hjust = 0.5)
                  )
  )
tiff(filename = "cenEI.T1T2.tiff", width = 12,height = 10,units = "in",res = 300,compression = "lzw")
print(cen.EI.T1T2)
dev.off()

LB_CPH_T1T2_Centality_standardized <- centralityTable(adjMat.T1T2.q, standardized = T)
write.csv(LB_CPH_T1T2_Centality_standardized, "LB_CPH_T1T2_Centality_standardized.csv")

#T2 to T3
tiff(filename = "cenALL.T2T3.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
centralityPlot(adjMat.T2T3.q,include = c("all"),scale = "z-scores") +
  ggtitle(cen.title.T2T3) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  )

dev.off()


cen.outEI.T2T3<- centralityPlot(adjMat.T2T3.q,include = c("OutExpectedInfluence"),scale = "z-scores",orderBy = "OutExpectedInfluence") +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  )
cen.inEI.T2T3<- centralityPlot(adjMat.T2T3.q,include = c("InExpectedInfluence"),scale = "z-scores",orderBy = "InExpectedInfluence") +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  )
cen.EI.T2T3 <- cen.outEI.T2T3 + cen.inEI.T2T3 +
  plot_annotation(title = cen.title.T2T3,
                  theme = theme(
                    plot.title = element_text(family = "serif", face = "bold", size = 16, hjust = 0.5)
                  )
  )
tiff(filename = "cenEI.T2T3.tiff", width = 12,height = 10,units = "in",res = 300,compression = "lzw")
print(cen.EI.T2T3)
dev.off()

LB_CPH_T2T3_Centality_standardized <- centralityTable(adjMat.T2T3.q, standardized = T)
write.csv(LB_CPH_T2T3_Centality_standardized, "LB_CPH_T2T3_Centality_standardized.csv")
# Calculate correlation stability coefficients
sink("nodeStability.T1T2.txt")
print(corStability(adjMatT1T2.nboot))
sink()

sink("nodeStability.T2T3.txt")
print(corStability(adjMatT2T3.nboot))
sink()


# Plot node centrality stability
tiff(filename = "corStability.T1T2.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T1T2.nboot, statistics=c("outExpectedInfluence", "inExpectedInfluence")) +
  ggtitle(LabwArrow.T1T2) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8),
        axis.title = element_text(size=8),
        text = element_text(family = "serif"))
dev.off()

tiff(filename = "corStability.T2T3.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T2T3.nboot, statistics=c("outExpectedInfluence", "inExpectedInfluence")) +
  ggtitle(LabwArrow.T2T3) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8),
        axis.title = element_text(size=8),
        text = element_text(family = "serif"))
dev.off()

# Bootstrap difference tests for node centrality
##T1T2Network##
tiff(filename = "node.outSdif.T1T2.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T1T2.eboot, "outStrength", plot="difference", order="sample") +
  ggtitle(LabwArrow.T1T2) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()
tiff(filename = "node.inSdif.T1T2.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T1T2.eboot, "inStrength", plot="difference", order="sample")+
  ggtitle(LabwArrow.T1T2) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()
tiff(filename = "node.outEIdif.T1T2.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T1T2.eboot, "outExpectedInfluence", plot="difference", order="sample")+
  ggtitle(LabwArrow.T1T2) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()
tiff(filename = "node.inEIdif.T1T2.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T1T2.eboot, "inExpectedInfluence", plot="difference", order="sample")+
  ggtitle(LabwArrow.T1T2) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()

##T2T3Network##
tiff(filename = "node.outSdif.T2T3.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T2T3.eboot, "outStrength", plot="difference", order="sample") +
  ggtitle(LabwArrow.T2T3) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()
tiff(filename = "node.inSdif.T2T3.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T2T3.eboot, "inStrength", plot="difference", order="sample")+
  ggtitle(LabwArrow.T2T3) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()
tiff(filename = "node.outEIdif.T2T3.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T2T3.eboot, "outExpectedInfluence", plot="difference", order="sample")+
  ggtitle(LabwArrow.T2T3) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()
tiff(filename = "node.inEIdif.T2T3.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(adjMat.T2T3.eboot, "inExpectedInfluence", plot="difference", order="sample")+
  ggtitle(LabwArrow.T2T3) +
  theme(plot.title=element_text(hjust=0.5,face="bold",size=16),
        axis.text = element_text(size = 8), axis.title = element_text(size=8))
dev.off()

# ====CLPN predictability and influence====
# ====Calculate node predictability====
#T1 to T2
pred012 <- rep(0, k)  # Model 0 includes autoregression for T1 to T2
pred112 <- rep(0, k)  # Model 1 excludes autoregression
pred212 <- rep(0, k)  # Model 2 uses predictors from other constructs
pred312 <- rep(0, k)  # Model 3 uses other predictors from the same construct

myData = LBCPHT1T2  # Select the data for the corresponding wave pair
for (i in 1:k){  #loop through each variable
  name.i <- names(myData)[k+i]  # Select the ith later-wave node as the outcome
  group.i <- groups[i]  # Identify the construct of the outcome node
  # Define the predictor sets
  include.mod0 <- names(myData)[1:k]  # Use all earlier-wave nodes, including the autoregressive predictor
  include.mod1 <- names(myData)[1:k][-i]  # Use all earlier-wave nodes except the autoregressive predictor
  include.mod2 <- names(myData)[1:k][groups != group.i]  # Use only predictors from other constructs
  target_indices <- (groups == group.i) & (1:k != i)
  include.mod3 <- names(myData)[1:k][target_indices] #include predictors belonging to same construct
  # Extract the corresponding LASSO coefficients from the fitted matrix
  betas.mod0 <- adjMat.T1T2[1:k,i]  # Fix the coefficients at their glmnet estimates
  betas.mod1 <- adjMat.T1T2[1:k,i][-i]
  betas.mod2 <- adjMat.T1T2[1:k,i][groups != group.i]
  betas.mod3 <- adjMat.T1T2[1:k,i][target_indices]
  # Build the regression formulas
  mod0 <- paste(name.i, " ~ ", paste(paste(betas.mod0,'*', include.mod0, sep = ""), collapse = "+"))
  mod1 <- paste(name.i, " ~ ", paste(paste(betas.mod1,'*', include.mod1, sep = ""), collapse = "+"))
  mod2 <- paste(name.i, " ~ ", paste(paste(betas.mod2,'*', include.mod2, sep = ""), collapse = "+"))
  mod3 <- paste(name.i, " ~ ", paste(paste(betas.mod3,'*', include.mod3, sep = ""), collapse = "+"))
  # Fit the SEMs and extract R-squared
  fit0 <- sem(mod0, sample.cov = cor(myData), sample.nobs = nrow(myData)) #fit model with subset of predictors to get r-square
  pred012[i] <- inspect(fit0, "rsquare")
  fit1 <- sem(mod1, sample.cov = cor(myData), sample.nobs = nrow(myData)) #fit model with subset of predictors to get r-square
  pred112[i] <- inspect(fit1, "rsquare")
  fit2 <- sem(mod2, sample.cov = cor(myData), sample.nobs = nrow(myData)) #fit model with subset of predictors to get r-square
  pred212[i] <- inspect(fit2, "rsquare")
  fit3 <- sem(mod3, sample.cov = cor(myData), sample.nobs = nrow(myData))
  pred312[i] <- inspect(fit3, "rsquare")
}

print("Predictability: Include self at previous time point")
pred012
print("Predictability: Exclude self at previous time point")
pred112
print("Predictability: Exclude predictors belonging to same construct")
pred212
print("Predictability: Include predictors belonging to same construct")
pred312


#T2 to T3
pred023 <- rep(0, k)  # Model 0 includes autoregression for T2 to T3
pred123 <- rep(0, k)  # Model 1 excludes autoregression
pred223 <- rep(0, k)  # Model 2 uses predictors from other constructs
pred323 <- rep(0, k)  # Model 3 uses other predictors from the same construct

myData = LBCPHT2T3  # Select the data for the corresponding wave pair
for (i in 1:k){  #loop through each variable
  name.i <- names(myData)[k+i]  # Select the ith later-wave node as the outcome
  group.i <- groups[i]  # Identify the construct of the outcome node
  # Define the predictor sets
  include.mod0 <- names(myData)[1:k]  # Use all earlier-wave nodes, including the autoregressive predictor
  include.mod1 <- names(myData)[1:k][-i]  # Use all earlier-wave nodes except the autoregressive predictor
  include.mod2 <- names(myData)[1:k][groups != group.i]  # Use only predictors from other constructs
  target_indices <- (groups == group.i) & (1:k != i)
  include.mod3 <- names(myData)[1:k][target_indices] #include predictors belonging to same construct
  # Extract the corresponding LASSO coefficients from the fitted matrix
  betas.mod0 <- adjMat.T2T3[1:k,i]  # Fix the coefficients at their glmnet estimates
  betas.mod1 <- adjMat.T2T3[1:k,i][-i]
  betas.mod2 <- adjMat.T2T3[1:k,i][groups != group.i]
  betas.mod3 <- adjMat.T2T3[1:k,i][target_indices]
  # Build the regression formulas
  mod0 <- paste(name.i, " ~ ", paste(paste(betas.mod0,'*', include.mod0, sep = ""), collapse = "+"))
  mod1 <- paste(name.i, " ~ ", paste(paste(betas.mod1,'*', include.mod1, sep = ""), collapse = "+"))
  mod2 <- paste(name.i, " ~ ", paste(paste(betas.mod2,'*', include.mod2, sep = ""), collapse = "+"))
  mod3 <- paste(name.i, " ~ ", paste(paste(betas.mod3,'*', include.mod3, sep = ""), collapse = "+"))
  # Fit the SEMs and extract R-squared
  fit0 <- sem(mod0, sample.cov = cor(myData), sample.nobs = nrow(myData)) #fit model with subset of predictors to get r-square
  pred023[i] <- inspect(fit0, "rsquare")
  fit1 <- sem(mod1, sample.cov = cor(myData), sample.nobs = nrow(myData)) #fit model with subset of predictors to get r-square
  pred123[i] <- inspect(fit1, "rsquare")
  fit2 <- sem(mod2, sample.cov = cor(myData), sample.nobs = nrow(myData)) #fit model with subset of predictors to get r-square
  pred223[i] <- inspect(fit2, "rsquare")
  fit3 <- sem(mod3, sample.cov = cor(myData), sample.nobs = nrow(myData))
  pred323[i] <- inspect(fit3, "rsquare")
}

print("Predictability: Include self at previous time point")
pred023
print("Predictability: Exclude self at previous time point")
pred123
print("Predictability: Exclude predictors belonging to same construct")
pred223
print("Predictability: Include predictors belonging to same construct")
pred323

# ====Calculate outgoing influence====
#T1 to T2
adjMat.noGroup.T1T2 <- adjMat.noAR.T1T2  # Start from the matrix with autoregression excluded
# Matrix rows represent earlier-wave predictors
# Matrix columns represent later-wave outcomes
adjMat.noGroup.T1T2[groups[row(adjMat.noGroup.T1T2)] == groups[col(adjMat.noGroup.T1T2)]] <- 0  # Set within-construct coefficients to zero

adjMat.within.T1T2 <- adjMat.noAR.T1T2  # Retain only within-construct coefficients after masking
adjMat.within.T1T2[groups[row(adjMat.within.T1T2)] != groups[col(adjMat.within.T1T2)]] <- 0

infl012 <- apply(adjMat.T1T2^2, 1, sum) #all outgoing nodes including AR
infl112 <- apply(adjMat.noAR.T1T2^2, 1, sum) #all outgoing nodes excluding AR
infl212 <- apply(adjMat.noGroup.T1T2^2, 1, sum)  # Sum squared cross-construct coefficients for each earlier-wave predictor
infl312 <- apply(adjMat.within.T1T2^2, 1, sum)

print("Influence: all outgoing nodes including AR")
infl012
print("Influence: all outgoing nodes excluding AR")
infl112
print("Influence: rows = T1 predictors; columns = T2 DVs")
infl212
print("Influence: rows = T1 predictors; columns = T2 DVs")
infl312

#T2 to T3
adjMat.noGroup.T2T3 <- adjMat.noAR.T2T3  # Start from the matrix with autoregression excluded
# Matrix rows represent earlier-wave predictors
# Matrix columns represent later-wave outcomes
adjMat.noGroup.T2T3[groups[row(adjMat.noGroup.T2T3)] == groups[col(adjMat.noGroup.T2T3)]] <- 0  # Set within-construct coefficients to zero

adjMat.within.T2T3 <- adjMat.noAR.T2T3  # Retain only within-construct coefficients after masking
adjMat.within.T2T3[groups[row(adjMat.within.T2T3)] != groups[col(adjMat.within.T2T3)]] <- 0

infl023 <- apply(adjMat.T2T3^2, 1, sum) #all outgoing nodes including AR
infl123 <- apply(adjMat.noAR.T2T3^2, 1, sum) #all outgoing nodes excluding AR
infl223 <- apply(adjMat.noGroup.T2T3^2, 1, sum)  # Sum squared cross-construct coefficients for each earlier-wave predictor
infl323 <- apply(adjMat.within.T2T3^2, 1, sum)

print("Influence: all outgoing nodes including AR")
infl023
print("Influence: all outgoing nodes excluding AR")
infl123
print("Influence: rows = T1 predictors; columns = T2 DVs")
infl223
print("Influence: rows = T1 predictors; columns = T2 DVs")
infl323

# Export predictability and influence statistics
predinfl.T1T2 <- data.frame(
  Node = labelnames,  # Node label
  Group = groups,  # Construct membership
  # --- Predictability  ---
  # Variance explained in each later-wave node
  Pred_Total = pred012,  # Includes autoregression
  Pred_CrossLag = pred112,  # Excludes autoregression
  Pred_CrossConstruct = pred212,  # Cross-construct predictability
  Pred_WithinConstruct = pred312,  # Within-construct predictability
  
  # --- Influence  ---
  # Sum of squared outgoing coefficients for each earlier-wave node
  Infl_Total_AR = infl012,  # Includes autoregression
  Infl_Spillover = infl112,  # Excludes autoregression
  Infl_CrossConstruct = infl212,  # Cross-construct outgoing influence
  Infl_WithinConstruct = infl312  # Within-construct outgoing influence
)
write.csv(predinfl.T1T2, "predinfl.T1T2.csv", row.names = FALSE)

predinfl.T2T3 <- data.frame(
  Node = labelnames,  # Node label
  Group = groups,  # Construct membership
  # --- Predictability  ---
  # Variance explained in each later-wave node
  Pred_Total = pred023,  # Includes autoregression
  Pred_CrossLag = pred123,  # Excludes autoregression
  Pred_CrossConstruct = pred223,  # Cross-construct predictability
  Pred_WithinConstruct = pred323,  # Within-construct predictability
  
  # --- Influence  ---
  # Sum of squared outgoing coefficients for each earlier-wave node
  Infl_Total_AR = infl023,  # Includes autoregression
  Infl_Spillover = infl123,  # Excludes autoregression
  Infl_CrossConstruct = infl223,  # Cross-construct outgoing influence
  Infl_WithinConstruct = infl323  # Within-construct outgoing influence
)
write.csv(predinfl.T2T3, "predinfl.T2T3.csv", row.names = FALSE)

# ====Plot node predictability====
#create a plot of predictability and influence:
# Define construct colors
construct_colors <- c("#FFFF33", "#a8e6cf", "#FC8D62", "#74ADD1")
raw_colors <- c(rep(construct_colors[1], k1),
                rep(construct_colors[2], k2),
                rep(construct_colors[3], k3),
                rep(construct_colors[4], k4))
group_names <-    c("Smartphone Usage Time",
                    "CES-D",
                    "Positive Suicide Ideation",
                    "Negative Suicide Ideation")

names(construct_colors) <- group_names
# Set the display order
#T1 to T2
pred12.title <- expression("T1" %->% "T2 Network Predictability")
pred012_order <- order(pred012)
pred012_sorted <- pred012[pred012_order]
labelnames012_sorted <- labelnames[pred012_order]
pred012_color <- raw_colors[pred012_order]
df_pred012 <- data.frame(
  Node = labelnames,  # Node label
  Group = groups,  # Construct membership
  Value = pred012  # Predictability estimate
)

pred112_order <- order(pred112)
pred112_sorted <- pred112[pred112_order]
labelnames112_sorted <- labelnames[pred112_order]
pred112_color <- raw_colors[pred112_order]
df_pred112 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred112
)

pred212_order <- order(pred212)
pred212_sorted <- pred212[pred212_order]
labelnames212_sorted <- labelnames[pred212_order]
pred212_color <- raw_colors[pred212_order]
df_pred212 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred212
)

pred312_order <- order(pred312)
pred312_sorted <- pred312[pred312_order]
labelnames312_sorted <- labelnames[pred312_order]
pred312_color <- raw_colors[pred312_order]
df_pred312 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred312
)
#T2 to T3
pred23.title <- expression("T2" %->% "T3 Network Predictability")
pred023_order <- order(pred023)
pred023_sorted <- pred023[pred023_order]
labelnames023_sorted <- labelnames[pred023_order]
pred023_color <- raw_colors[pred023_order]
df_pred023 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred023
)

pred123_order <- order(pred123)
pred123_sorted <- pred123[pred123_order]
labelnames123_sorted <- labelnames[pred123_order]
pred123_color <- raw_colors[pred123_order]
df_pred123 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred123
)

pred223_order <- order(pred223)
pred223_sorted <- pred223[pred223_order]
labelnames223_sorted <- labelnames[pred223_order]
pred223_color <- raw_colors[pred223_order]
df_pred223 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred223
)

pred323_order <- order(pred323)
pred323_sorted <- pred323[pred323_order]
labelnames323_sorted <- labelnames[pred323_order]
pred323_color <- raw_colors[pred323_order]
df_pred323 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = pred323
)
# Scale the plotting opacity within each construct
# Define the opacity scaling for predictability and squared influence
calc_group_alpha <- function(df) {
  df %>%
    group_by(Group) %>%
    mutate(
      # Calculate the within-construct minimum and maximum
      min_v = min(Value),
      max_v = max(Value),
      # Rescale to the interval from 0 to 1
      # Handle groups with identical values
      Alpha_Score = ifelse(max_v == min_v,
                           1,
                           (Value - min_v) / (max_v - min_v))
    ) %>%
    ungroup()  # Remove grouping
}
# Define separate opacity scaling for positive and negative EI values
calc_EI_alpha <- function(df) {
  
  stopifnot(
    all(c("Group", "Value") %in% names(df)),
    is.numeric(df$Value),
    all(is.finite(df$Value))
  )
  
  df %>%
    group_by(Group, .EI_sign = sign(Value)) %>%
    mutate(
      Alpha_Score = {
        max_abs <- max(abs(Value))
        
        if (max_abs == 0) {
          rep(0, length(Value))
        } else {
          abs(Value) / max_abs
        }
      }
    ) %>%
    ungroup() %>%
    select(-.EI_sign)
}

# Scale the plotting opacity within each construct
plot_df012_pred <- calc_group_alpha(df_pred012)
plot_df112_pred <- calc_group_alpha(df_pred112)
plot_df212_pred <- calc_group_alpha(df_pred212)
plot_df312_pred <- calc_group_alpha(df_pred312)

plot_df023_pred <- calc_group_alpha(df_pred023)
plot_df123_pred <- calc_group_alpha(df_pred123)
plot_df223_pred <- calc_group_alpha(df_pred223)
plot_df323_pred <- calc_group_alpha(df_pred323)
# Draw the plot
# T1 -> T2
# Set the panel layout
common_theme.T1T2 <- list(
  geom_col(width = 1, color = "black", linewidth = 0.1),
  coord_flip(),
  scale_fill_manual(values = construct_colors),
  scale_alpha_continuous(range = c(0.4, 1), guide = "none"),
  scale_y_continuous(
    limits = c(0, 0.25),
    expand = c(0, 0),
    breaks = seq(0, 0.25, 0.05),
    guide = guide_axis(cap = "both")  # Cap the axis line at both ends
  ),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    axis.line.x = element_line(size = 0.8, color = "black", lineend = "square"),
    axis.ticks.x = element_line(size = 0.8, color = "black"),
    axis.ticks.length.x = unit(-0.15, "cm")
  )
)

pred012.p <- ggplot(plot_df012_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T1T2 +
  labs(title = "Total Predictability", y = "", x = "") +
  theme(
    axis.text.x = element_blank(),  # Hide tick labels
    axis.ticks.x = element_blank(),  # Hide tick marks
    axis.line.x = element_blank()  # Hide the axis line
  )

pred112.p <- ggplot(plot_df112_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T1T2 +
  labs(title = "Predictability without AR", y = "", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )
pred212.p <- ggplot(plot_df212_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T1T2 +
  labs(title = "Cross-construct Predictability", y = "", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )
pred312.p <- ggplot(plot_df312_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T1T2 +
  labs(title = "Within-construct Predictability",
       y = "", x = "") +
  theme(
    axis.text.x = element_text(size = 10, color = "black", margin = margin(t = 5)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t=10))
  )
pred.T1T2.p <- (pred012.p / pred112.p / pred212.p / pred312.p) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = pred12.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0.5))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 14),
    legend.text = element_text(family = "serif", size = 12),
    # Set spacing between panels
    plot.margin = margin(5, 5, 5, 5)
  )

# Save the figure
ggsave("pred.T1T2.tiff", plot = pred.T1T2.p,
       width = 20, height = 20,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")


#T2 -> T3
# Set the panel layout
common_theme.T2T3 <- list(
  geom_col(width = 1, color = "black", linewidth = 0.1),
  coord_flip(),
  scale_fill_manual(values = construct_colors),
  scale_alpha_continuous(range = c(0.4, 1), guide = "none"),
  scale_y_continuous(
    limits = c(0, 0.35),
    expand = c(0, 0),
    breaks = seq(0, 0.35, 0.05),
    guide = guide_axis(cap = "both")  # Cap the axis line at both ends
  ),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    axis.line.x = element_line(size = 0.8, color = "black", lineend = "square"),
    axis.ticks.x = element_line(size = 0.8, color = "black"),
    axis.ticks.length.x = unit(-0.15, "cm")
  )
)

# P1 Total Predictability
pred023.p <- ggplot(plot_df023_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T2T3 +
  labs(title = "Total Predictability", y = "", x = "") +
  theme(
    axis.text.x = element_blank(),  # Hide tick labels
    axis.ticks.x = element_blank(),  # Hide tick marks
    axis.line.x = element_blank()  # Hide the axis line
  )

# P2: Predictability withoutAR
pred123.p <- ggplot(plot_df123_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T2T3 +
  labs(title = "Predictability without AR", y = "", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

# P3: Cross-construct Predictability
pred223.p <- ggplot(plot_df223_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T2T3 +
  labs(title = "Cross-construct Predictability", y = "", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank()
  )

# P4: Within-construct Predictability
pred323.p <- ggplot(plot_df323_pred, aes(x = reorder(Node, Value), y = Value, fill = Group, alpha = Alpha_Score)) +
  common_theme.T2T3 +
  labs(title = "Within-construct Predictability",
       y = "", x = "") +
  theme(
    axis.text.x = element_text(size = 10, color = "black", margin = margin(t = 5)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t=10))
  )

pred.T2T3.p <- (pred023.p / pred123.p / pred223.p / pred323.p) +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = pred23.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0.5))
  ) &
  theme(
    # Place the legend on the right
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 14),
    legend.text = element_text(family = "serif", size = 12),
    # Set spacing between panels
    plot.margin = margin(5, 5, 5, 5)
  )

# Save the figure
ggsave("pred.T2T3.tiff", plot = pred.T2T3.p,
       width = 20, height = 20,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")

# ====Plot predictability in a rectangular layout====
# T1 to T2
common_theme.pred.T1T2.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(range = c(0.4, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

pred012.m <- ggplot(plot_df012_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.25),  # Set the axis limits
    breaks = seq(0, 0.25, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total Predictability", y = "", x = "")


pred112.m <- ggplot(plot_df112_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.2),  # Set the axis limits
    breaks = seq(0, 0.2, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Predictability withoutAR", y = "", x = "")

pred212.m <- ggplot(plot_df212_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.05),  # Set the axis limits
    breaks = seq(0, 0.05, 0.01),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct Predictability", y = "", x = "")

pred312.m <- ggplot(plot_df312_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.15),  # Set the axis limits
    breaks = seq(0, 0.15, 0.03),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct Predictability", y = "", x = "")

pred.T1T2.Matrix <- (pred012.m + pred112.m + pred212.m + pred312.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = pred12.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("pred.T1T2.Matrix.tiff", plot = pred.T1T2.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")

# T2 to T3
common_theme.pred.T2T3.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(range = c(0.4, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

pred023.m <- ggplot(plot_df023_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.35),  # Set the axis limits
    breaks = seq(0, 0.35, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total Predictability", y = "", x = "")


pred123.m <- ggplot(plot_df123_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.20),  # Set the axis limits
    breaks = seq(0, 0.20, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Predictability withoutAR", y = "", x = "")

pred223.m <- ggplot(plot_df223_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.06),  # Set the axis limits
    breaks = seq(0, 0.06, 0.01),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct Predictability", y = "", x = "")

pred323.m <- ggplot(plot_df323_pred, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.pred.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.15),  # Set the axis limits
    breaks = seq(0, 0.15, 0.03),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct Predictability", y = "", x = "")

pred.T2T3.Matrix <- (pred023.m + pred123.m + pred223.m + pred323.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = pred23.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("pred.T2T3.Matrix.tiff", plot = pred.T2T3.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")

# ====Plot CLPNs with predictability rings====
tiff("LB_T1T3_CPH_Network_withPred.tiff", width = 24, height = 12, units = "in",res = 300, compression = "lzw", family = "serif")
lay_mat <- matrix(c(1, 2,
                    3, 3), nrow = 2, byrow = TRUE)
layout(lay_mat, heights = c(0.75, 0.25))
par(family = "serif", font.main = 2)

adjMat.T1T2.q <- qgraph(adjMat.T1T2, groups = groups, labels = labelnames,
                        legend = F,vsize=6,cut=0,
                        border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                        nodeNames=longnames,legend.cex=.3,GLratio = 2,
                        layout = "spring",
                        title = LabwArrow.T1T2, title.cex = 2, pie = pred012)

adjMat.T2T3.q <- qgraph(adjMat.T2T3, groups = groups, labels = labelnames,
                        legend = F,vsize=6,cut=0,
                        border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                        nodeNames=longnames,legend.cex=.3,GLratio = 2,
                        layout = adjMat.T1T2.q$layout,
                        title = LabwArrow.T2T3, title.cex = 2, pie = pred023)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Node legend title
text(x = 1, y = 1.35,  # Set the vertical position of the legend title
     labels = "Nodes",
     font = 2,  # Use bold text
     cex = 2)  # Set the text size
legend("center",
       legend = legend_text,  # Node legend text
       
       # Node legend markers
       pch = 21,  # Use filled circles
       pt.bg = node_colors,  # Match the fill to the node color
       pt.cex = 1.5,  # Set the marker size
       col = "black",  # Set the marker border color
       # -----------------------
       
       bty = "n",  # Omit the legend box
       cex = 1.2,  # Set the legend text size
       ncol = 4,  # Arrange the legend in columns
       x.intersp = 0.8,  # Set the spacing between symbols and labels
       text.width = 0.2 )
legend("bottom",
       legend = c("Smartphone Usage Time", "CES-D", "Positive Suicide Ideation", "Negative Suicide Ideation"),
       pt.bg = colors_palette,  # Use one color per construct
       pch = 22,  # Use square markers
       col = "black", bty = "n", horiz = TRUE, cex = 1.2, pt.cex = 2)
dev.off()

# Hide autoregressive edges
tiff("LB_T1T3_CPH_Network_noAR_withPred.tiff", width = 24, height = 12, units = "in",res = 300, compression = "lzw", family = "serif")
lay_mat <- matrix(c(1, 2,
                    3, 3), nrow = 2, byrow = TRUE)
layout(lay_mat, heights = c(0.75, 0.25))
par(family = "serif", font.main = 2)

adjMat.noAR.T1T2.q <- qgraph(adjMat.noAR.T1T2, groups = groups, labels = labelnames,
                             legend = F,vsize=6,cut=0, minimum = quantile(abs(adjMat.noAR.T1T2[adjMat.noAR.T1T2 != 0]), 0.80),
                             border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                             nodeNames=longnames,legend.cex=.3,GLratio = 2,
                             layout = adjMat.T1T2.q$layout,
                             title = LabwArrow.T1T2, title.cex = 2, pie = pred012)

adjMat.noAR.T2T3.q <- qgraph(adjMat.noAR.T2T3, groups = groups, labels = labelnames,
                             legend = F,vsize=6,cut=0, minimum = quantile(abs(adjMat.noAR.T2T3[adjMat.noAR.T2T3 != 0]), 0.80),
                             border.color="black",color = c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                             nodeNames=longnames,legend.cex=.3,GLratio = 2,
                             layout = adjMat.T1T2.q$layout,
                             title = LabwArrow.T2T3, title.cex = 2, pie = pred023)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Node legend title
text(x = 1, y = 1.35,  # Set the vertical position of the legend title
     labels = "Nodes",
     font = 2,  # Use bold text
     cex = 2)  # Set the text size
legend("center",
       legend = legend_text,  # Node legend text
       pch = 21,  # Use filled circles
       pt.bg = node_colors,  # Match the fill to the node color
       pt.cex = 1.5,  # Set the marker size
       col = "black",  # Set the marker border color
       bty = "n",  # Omit the legend box
       cex = 1.2,  # Set the legend text size
       ncol = 4,  # Arrange the legend in columns
       x.intersp = 0.8,  # Set the spacing between symbols and labels
       text.width = 0.2 )
# Construct legend
legend("bottom",
       legend = c("Smartphone Usage Time", "CES-D", "Positive Suicide Ideation", "Negative Suicide Ideation"),
       pt.bg = colors_palette,  # Use one color per construct
       pch = 22,  # Use square markers
       col = "black", bty = "n", horiz = TRUE, cex = 1.2, pt.cex = 2)
dev.off()




# ====Plot outgoing influence in a rectangular layout====
# T1 to T2
common_theme.infl.T1T2.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(range = c(0.4, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

infl012.m <- ggplot(plot_df012_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.2),  # Set the axis limits
    breaks = seq(0, 0.2, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total Influence", y = "", x = "")


infl112.m <- ggplot(plot_df112_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.1),  # Set the axis limits
    breaks = seq(0, 0.1, 0.02),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Influence withoutAR", y = "", x = "")

infl212.m <- ggplot(plot_df212_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.05),  # Set the axis limits
    breaks = seq(0, 0.05, 0.01),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct Influence", y = "", x = "")

infl312.m <- ggplot(plot_df312_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T1T2.m +
  scale_y_continuous(
    limits = c(0, 0.1),  # Set the axis limits
    breaks = seq(0, 0.1, 0.02),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct Influence", y = "", x = "")

infl.T1T2.Matrix <- (infl012.m + infl112.m + infl212.m + infl312.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = infl12.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("infl.T1T2.Matrix.tiff", plot = infl.T1T2.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")

# T2 to T3
common_theme.infl.T2T3.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(range = c(0.4, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10, color = "black", margin = margin(r = 5)),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

infl023.m <- ggplot(plot_df023_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.35),  # Set the axis limits
    breaks = seq(0, 0.35, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total Influence", y = "", x = "")


infl123.m <- ggplot(plot_df123_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.15),  # Set the axis limits
    breaks = seq(0, 0.15, 0.03),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Influence withoutAR", y = "", x = "")

infl223.m <- ggplot(plot_df223_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.05),  # Set the axis limits
    breaks = seq(0, 0.05, 0.01),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct Influence", y = "", x = "")

infl323.m <- ggplot(plot_df323_infl, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.infl.T2T3.m +
  scale_y_continuous(
    limits = c(0, 0.1),  # Set the axis limits
    breaks = seq(0, 0.1, 0.02),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct Influence", y = "", x = "")

infl.T2T3.Matrix <- (infl023.m + infl123.m + infl223.m + infl323.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = infl23.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("infl.T2T3.Matrix.tiff", plot = infl.T2T3.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")
# ====Calculate incoming and outgoing expected influence====
# T1 to T2
adjMat.T1T2.EI <- adjMat.T1T2
# Exclude autoregressive edges
adjMat.T1T2.EI.noAR <- adjMat.T1T2.EI
diag(adjMat.T1T2.EI.noAR) <- 0

# Cross-construct coefficient matrix
adjMat.T1T2.EI.cross <- adjMat.T1T2.EI.noAR
adjMat.T1T2.EI.cross[groups[row(adjMat.T1T2.EI.cross)] == groups[col(adjMat.T1T2.EI.cross)]] <- 0

# Within-construct coefficient matrix
adjMat.T1T2.EI.within <- adjMat.T1T2.EI.noAR
adjMat.T1T2.EI.within[groups[row(adjMat.T1T2.EI.within)] != groups[col(adjMat.T1T2.EI.within)]] <- 0

# Calculate outgoing expected influence
outEI0.T1T2  <- rowSums(adjMat.T1T2.EI)  # Includes autoregression
outEI1.T1T2   <- rowSums(adjMat.T1T2.EI.noAR)  # Excludes autoregression
outEI2.T1T2  <- rowSums(adjMat.T1T2.EI.cross)  # Outgoing expected influence across constructs
outEI3.T1T2 <- rowSums(adjMat.T1T2.EI.within)  # Outgoing expected influence within constructs

# Calculate incoming expected influence
inEI0.T1T2  <- colSums(adjMat.T1T2.EI)
inEI1.T1T2   <- colSums(adjMat.T1T2.EI.noAR)
inEI2.T1T2  <- colSums(adjMat.T1T2.EI.cross)
inEI3.T1T2 <- colSums(adjMat.T1T2.EI.within)

# Export results to Excel
T1T2_allEI_raw <- data.frame(
  Node = labelnames,
  Group = groups,
  inEI_Total  = inEI0.T1T2,  # 0 = total
  inEI_noAR   = inEI1.T1T2,  # 1 = excluding autoregression
  inEI_Cross  = inEI2.T1T2,  # 2 = cross-construct
  inEI_Within = inEI3.T1T2,  # 3 = within-construct
  outEI_Total  = outEI0.T1T2,  # 0 = total
  outEI_noAR   = outEI1.T1T2,  # 1 = excluding autoregression
  outEI_Cross  = outEI2.T1T2,  # 2 = cross-construct
  outEI_Within = outEI3.T1T2  # 3 = within-construct
)
write_xlsx(T1T2_allEI_raw, "T1T2_Network_All_EI_Indices.xlsx")


# Set the display order
inEI012.title <- expression(bold("T1" %->% "T2 Network In-Expected Influence"))

df_inEI0.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI0.T1T2
)

df_inEI1.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI1.T1T2
)

df_inEI2.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI2.T1T2
)

df_inEI3.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI3.T1T2
)

outEI012.title <- expression(bold("T1" %->% "T2 Network Out-Expected Influence"))

df_outEI0.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI0.T1T2
)

df_outEI1.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI1.T1T2
)

df_outEI2.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI2.T1T2
)

df_outEI3.T1T2 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI3.T1T2
)
# Scale the plotting opacity within each construct
plot_df012_outEI <- calc_EI_alpha(df_outEI0.T1T2)
plot_df112_outEI <- calc_EI_alpha(df_outEI1.T1T2)
plot_df212_outEI <- calc_EI_alpha(df_outEI2.T1T2)
plot_df312_outEI <- calc_EI_alpha(df_outEI3.T1T2)

plot_df012_inEI <- calc_EI_alpha(df_inEI0.T1T2)
plot_df112_inEI <- calc_EI_alpha(df_inEI1.T1T2)
plot_df212_inEI <- calc_EI_alpha(df_inEI2.T1T2)
plot_df312_inEI <- calc_EI_alpha(df_inEI3.T1T2)
# T2 to T3
adjMat.T2T3.EI <- adjMat.T2T3

adjMat.T2T3.EI.noAR <- adjMat.T2T3.EI
diag(adjMat.T2T3.EI.noAR) <- 0

adjMat.T2T3.EI.cross <- adjMat.T2T3.EI.noAR
adjMat.T2T3.EI.cross[groups[row(adjMat.T2T3.EI.cross)] == groups[col(adjMat.T2T3.EI.cross)]] <- 0

adjMat.T2T3.EI.within <- adjMat.T2T3.EI.noAR
adjMat.T2T3.EI.within[groups[row(adjMat.T2T3.EI.within)] != groups[col(adjMat.T2T3.EI.within)]] <- 0

outEI0.T2T3  <- rowSums(adjMat.T2T3.EI)  # Includes autoregression
outEI1.T2T3   <- rowSums(adjMat.T2T3.EI.noAR)  # Excludes autoregression
outEI2.T2T3  <- rowSums(adjMat.T2T3.EI.cross)  # Outgoing expected influence across constructs
outEI3.T2T3 <- rowSums(adjMat.T2T3.EI.within)  # Outgoing expected influence within constructs

inEI0.T2T3  <- colSums(adjMat.T2T3.EI)
inEI1.T2T3   <- colSums(adjMat.T2T3.EI.noAR)
inEI2.T2T3  <- colSums(adjMat.T2T3.EI.cross)
inEI3.T2T3 <- colSums(adjMat.T2T3.EI.within)


# Export results to Excel
T2T3_allEI_raw <- data.frame(
  Node = labelnames,
  Group = groups,
  inEI_Total  = inEI0.T2T3,  # 0 = total
  inEI_noAR   = inEI1.T2T3,  # 1 = excluding autoregression
  inEI_Cross  = inEI2.T2T3,  # 2 = cross-construct
  inEI_Within = inEI3.T2T3,  # 3 = within-construct
  outEI_Total  = outEI0.T2T3,  # 0 = total
  outEI_noAR   = outEI1.T2T3,  # 1 = excluding autoregression
  outEI_Cross  = outEI2.T2T3,  # 2 = cross-construct
  outEI_Within = outEI3.T2T3  # 3 = within-construct
)
write_xlsx(T2T3_allEI_raw, "T2T3_Network_All_EI_Indices.xlsx")

# Set the display order
inEI023.title <- expression(bold("T2" %->% "T3 Network In-Expected Influence"))

df_inEI0.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI0.T2T3
)

df_inEI1.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI1.T2T3
)

df_inEI2.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI2.T2T3
)

df_inEI3.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = inEI3.T2T3
)

outEI023.title <- expression(bold("T2" %->% "T3 Network Out-Expected Influence"))

df_outEI0.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI0.T2T3
)

df_outEI1.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI1.T2T3
)

df_outEI2.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI2.T2T3
)

df_outEI3.T2T3 <- data.frame(
  Node = labelnames,
  Group = groups,
  Value = outEI3.T2T3
)
# Scale the plotting opacity within each construct
plot_df023_outEI <- calc_EI_alpha(df_outEI0.T2T3)
plot_df123_outEI <- calc_EI_alpha(df_outEI1.T2T3)
plot_df223_outEI <- calc_EI_alpha(df_outEI2.T2T3)
plot_df323_outEI <- calc_EI_alpha(df_outEI3.T2T3)

plot_df023_inEI <- calc_EI_alpha(df_inEI0.T2T3)
plot_df123_inEI <- calc_EI_alpha(df_inEI1.T2T3)
plot_df223_inEI <- calc_EI_alpha(df_inEI2.T2T3)
plot_df323_inEI <- calc_EI_alpha(df_inEI3.T2T3)


# ====Plot incoming expected influence in a rectangular layout====
# T1 to T2
common_theme.inEI.T1T2.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  geom_text(aes(y = 0, label = Node,
                # Place labels on opposite sides of zero according to the sign
                hjust = ifelse(Value > 0, 1.1, -0.1)),
            size = 3.5,  # Set the text size
            family = "serif",
            color = "black",
            alpha = 1),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(limits = c(0, 1), range = c(0.15, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

inEI012.m <- ggplot(plot_df012_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T1T2.m +
  scale_y_continuous(
    limits = c(-0.6, 0.6),  # Set the axis limits
    breaks = seq(-0.6, 0.6, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total In-Expected Influence", y = "", x = "")


inEI112.m <- ggplot(plot_df112_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T1T2.m +
  scale_y_continuous(
    limits = c(-0.5, 0.5),  # Set the axis limits
    breaks = seq(-0.5, 0.5, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "In-Expected Influence withoutAR", y = "", x = "")

inEI212.m <- ggplot(plot_df212_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T1T2.m +
  scale_y_continuous(
    limits = c(-0.25, 0.25),  # Set the axis limits
    breaks = seq(-0.25, 0.25, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct In-Expected Influence", y = "", x = "")

inEI312.m <- ggplot(plot_df312_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T1T2.m +
  scale_y_continuous(
    limits = c(-0.4, 0.4),  # Set the axis limits
    breaks = seq(-0.4, 0.4, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct In-Expected Influence", y = "", x = "")

inEI.T1T2.Matrix <- (inEI012.m + inEI112.m + inEI212.m + inEI312.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = inEI012.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("inEI.T1T2.Matrix.tiff", plot = inEI.T1T2.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")

# T2 to T3
common_theme.inEI.T2T3.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  geom_text(aes(y = 0, label = Node,
                # Place labels on opposite sides of zero according to the sign
                hjust = ifelse(Value > 0, 1.1, -0.1)),
            size = 3.5,  # Set the text size
            family = "serif",
            color = "black",
            alpha = 1),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(limits = c(0, 1), range = c(0.15, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

inEI023.m <- ggplot(plot_df023_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T2T3.m +
  scale_y_continuous(
    limits = c(-0.7, 0.7),  # Set the axis limits
    breaks = seq(-0.7, 0.7, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total In-Expected Influence", y = "", x = "")


inEI123.m <- ggplot(plot_df123_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T2T3.m +
  scale_y_continuous(
    limits = c(-0.6, 0.6),  # Set the axis limits
    breaks = seq(-0.6, 0.6, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "In-Expected Influence withoutAR", y = "", x = "")

inEI223.m <- ggplot(plot_df223_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T2T3.m +
  scale_y_continuous(
    limits = c(-0.35, 0.35),  # Set the axis limits
    breaks = seq(-0.35, 0.35, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct In-Expected Influence", y = "", x = "")

inEI323.m <- ggplot(plot_df323_inEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.inEI.T2T3.m +
  scale_y_continuous(
    limits = c(-0.4, 0.4),  # Set the axis limits
    breaks = seq(-0.4, 0.4, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct In-Expected Influence", y = "", x = "")

inEI.T2T3.Matrix <- (inEI023.m + inEI123.m + inEI223.m + inEI323.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = inEI023.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("inEI.T2T3.Matrix.tiff", plot = inEI.T2T3.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")



# ====Plot outgoing expected influence in a rectangular layout====
# T1 to T2
common_theme.outEI.T1T2.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  geom_text(aes(y = 0, label = Node,
                # Place labels on opposite sides of zero according to the sign
                hjust = ifelse(Value > 0, 1.1, -0.1)),
            size = 3.5,  # Set the text size
            family = "serif",
            color = "black",
            alpha = 1),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(limits = c(0, 1), range = c(0.15, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

outEI012.m <- ggplot(plot_df012_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T1T2.m +
  scale_y_continuous(
    limits = c(-1.2, 1.2),  # Set the axis limits
    breaks = seq(-1.2, 1.2, 0.2),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total Out-Expected Influence", y = "", x = "")


outEI112.m <- ggplot(plot_df112_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T1T2.m +
  scale_y_continuous(
    limits = c(-1.0, 1.0),  # Set the axis limits
    breaks = seq(-1.0, 1.0, 0.2),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Out-Expected Influence withoutAR", y = "", x = "")

outEI212.m <- ggplot(plot_df212_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T1T2.m +
  scale_y_continuous(
    limits = c(-0.4, 0.4),  # Set the axis limits
    breaks = seq(-0.4, 0.4, 0.05),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct Out-Expected Influence", y = "", x = "")

outEI312.m <- ggplot(plot_df312_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T1T2.m +
  scale_y_continuous(
    limits = c(-0.7, 0.7),  # Set the axis limits
    breaks = seq(-0.7, 0.7, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct Out-Expected Influence", y = "", x = "")

outEI.T1T2.Matrix <- (outEI012.m + outEI112.m + outEI212.m + outEI312.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = outEI012.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("outEI.T1T2.Matrix.tiff", plot = outEI.T1T2.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")

# T2 to T3
common_theme.outEI.T2T3.m <- list(
  geom_col(aes(alpha = Alpha_Score),width = 1, color = "black",linewidth = 0.1),
  coord_flip(),
  geom_text(aes(y = 0, label = Node,
                # Place labels on opposite sides of zero according to the sign
                hjust = ifelse(Value > 0, 1.1, -0.1)),
            size = 3.5,  # Set the text size
            family = "serif",
            color = "black",
            alpha = 1),
  scale_fill_manual(values = construct_colors, name = "Construct"),
  scale_alpha_continuous(limits = c(0, 1), range = c(0.15, 1), guide = "none"),
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))),
  theme_classic(),
  theme(
    text = element_text(family = "serif", color = "grey20"),
    plot.title.position = "panel",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    
    # Style the axes
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    
    axis.line.x = element_line(size = 0.5, color = "grey", lineend = "square"),
    axis.ticks.x = element_line(size = 1, color = "grey"),
    axis.ticks.length.x = unit(-0.15, "cm"),
    axis.text.x = element_text(size = 11, color = "grey20", margin = margin(t = 5))
  )
)

outEI023.m <- ggplot(plot_df023_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T2T3.m +
  scale_y_continuous(
    limits = c(-1.5, 1.5),  # Set the axis limits
    breaks = seq(-1.5, 1.5, 0.3),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Total Out-Expected Influence", y = "", x = "")


outEI123.m <- ggplot(plot_df123_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T2T3.m +
  scale_y_continuous(
    limits = c(-1.2, 1.2),  # Set the axis limits
    breaks = seq(-1.2, 1.2, 0.2),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Out-Expected Influence withoutAR", y = "", x = "")

outEI223.m <- ggplot(plot_df223_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T2T3.m +
  scale_y_continuous(
    limits = c(-0.5, 0.5),  # Set the axis limits
    breaks = seq(-0.5, 0.5, 0.1),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Cross-construct Out-Expected Influence", y = "", x = "")

outEI323.m <- ggplot(plot_df323_outEI, aes(x = reorder(Node, Value), y = Value, fill = Group)) +
  common_theme.outEI.T2T3.m +
  scale_y_continuous(
    limits = c(-0.8, 0.8),  # Set the axis limits
    breaks = seq(-0.8, 0.8, 0.2),  # Set the tick positions
    expand = c(0, 0),
    guide = guide_axis(cap = "both"),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(title = "Within-construct Out-Expected Influence", y = "", x = "")

outEI.T2T3.Matrix <- (outEI023.m + outEI123.m + outEI223.m + outEI323.m) +
  plot_layout(ncol = 2,
              guides = "collect") +
  plot_annotation(
    title = outEI023.title,  # Set the overall figure title
    theme = theme(plot.title = element_text(family = "serif", size = 20, face = "bold", hjust = 0))
  ) &
  theme(
    # Set the legend appearance
    legend.position = "right",
    legend.justification = "center",  # Center the legend vertically
    legend.title = element_text(family = "serif", face = "bold", size = 16, color = "black"),
    legend.text = element_text(family = "serif", size = 12, color = "black"),
    # Set spacing between panels
    plot.margin = margin(10, 20, 10, 20)
  )

# Save the figure
ggsave("outEI.T2T3.Matrix.tiff", plot = outEI.T2T3.Matrix,
       width = 18, height = 14,  # Set figure dimensions
       dpi = 300, compression = "lzw", bg = "white")


# ====Compare the networks across intervals====
#Net12
# Correlate cross-lagged edge weights across intervals
adj_T1T2 <- adjMat.T1T2.e$graph
adj_T2T3 <- adjMat.T2T3.e$graph
compare.Edge.Net12 <- cor.test(c(adj_T1T2[lower.tri(adj_T1T2,diag=FALSE)],
                                 adj_T1T2[upper.tri(adj_T1T2, diag=FALSE)]), # excludes AR
                               c(adj_T2T3[lower.tri(adj_T2T3,diag=FALSE)],
                                 adj_T2T3[upper.tri(adj_T2T3, diag=FALSE)]), # excludes AR
                               method="pearson")
compare.Edge.Net12
# Correlate incoming expected influence across intervals
test_inEI_Total  <- cor.test(inEI0.T1T2, inEI0.T2T3, method="pearson")
test_inEI_noAR   <- cor.test(inEI1.T1T2, inEI1.T2T3, method="pearson")
test_inEI_Cross  <- cor.test(inEI2.T1T2, inEI2.T2T3, method="pearson")
test_inEI_Within <- cor.test(inEI3.T1T2, inEI3.T2T3, method="pearson")

compare.outEI.Net12 <- cor.test(centrality(adjMat.noAR.T1T2.q)$OutExpectedInfluence,
                                centrality(adjMat.noAR.T2T3.q)$OutExpectedInfluence, method="pearson")
compare.outEI.Net12
# Correlate outgoing expected influence across intervals
test_outEI_Total  <- cor.test(outEI0.T1T2, outEI0.T2T3, method="pearson")
test_outEI_noAR   <- cor.test(outEI1.T1T2, outEI1.T2T3, method="pearson")
test_outEI_Cross  <- cor.test(outEI2.T1T2, outEI2.T2T3, method="pearson")
test_outEI_Within <- cor.test(outEI3.T1T2, outEI3.T2T3, method="pearson")

compare.inEI.Net12 <- cor.test(centrality(adjMat.noAR.T1T2.q)$InExpectedInfluence,
                               centrality(adjMat.noAR.T2T3.q)$InExpectedInfluence,method="pearson")
compare.inEI.Net12
# Calculate signed Jaccard overlap
jac.Net12 <- (sum(adj_T1T2>0 & adj_T2T3>0) +
                sum(adj_T1T2<0 & adj_T2T3<0)) /
  sum(adj_T1T2 !=0 | adj_T2T3!=0)
jac.Net12

shared_edges <- (sum(adj_T1T2>0 & adj_T2T3>0) +
                   sum(adj_T1T2<0 & adj_T2T3<0))

edges_T1T2_total <- sum(adj_T1T2 != 0)
edges_T2T3_total <- sum(adj_T2T3 != 0)

retention_rate <- shared_edges / edges_T1T2_total
replication_rate <- shared_edges / edges_T2T3_total

Master_Stability_Table <- data.frame(
  Comparison = rep("T1 -> T2 vs T2 -> T3", 12),
  
  Metric = c(
    "1. Network Edges (no AR)",
    "2. Jaccard Index",
    "3. Edge Retention Rate",
    "4. Edge Replication Rate",
    "5. In-EI (Total)",
    "6. In-EI (no AR)",
    "7. In-EI (Cross-construct)",
    "8. In-EI (Within-construct)",
    "9. Out-EI (Total)",
    "10. Out-EI (no AR)",
    "11. Out-EI (Cross-construct)",
    "12. Out-EI (Within-construct)"
  ),
  
  Estimate_r_or_Value = c(
    compare.Edge.Net12$estimate,
    jac.Net12,
    retention_rate,  # Report the edge retention rate as a proportion
    replication_rate,  # Report the edge replication rate as a proportion
    test_inEI_Total$estimate,
    test_inEI_noAR$estimate,
    test_inEI_Cross$estimate,
    test_inEI_Within$estimate,
    test_outEI_Total$estimate,
    test_outEI_noAR$estimate,
    test_outEI_Cross$estimate,
    test_outEI_Within$estimate
  ),
  
  p_value = c(
    format.pval(compare.Edge.Net12$p.value, eps = 0.001, digits = 3),
    NA,
    NA,
    NA,
    format.pval(test_inEI_Total$p.value, eps = 0.001, digits = 3),
    format.pval(test_inEI_noAR$p.value, eps = 0.001, digits = 3),
    format.pval(test_inEI_Cross$p.value, eps = 0.001, digits = 3),
    format.pval(test_inEI_Within$p.value, eps = 0.001, digits = 3),
    format.pval(test_outEI_Total$p.value, eps = 0.001, digits = 3),
    format.pval(test_outEI_noAR$p.value, eps = 0.001, digits = 3),
    format.pval(test_outEI_Cross$p.value, eps = 0.001, digits = 3),
    format.pval(test_outEI_Within$p.value, eps = 0.001, digits = 3)
  ),
  
  row.names = NULL
)
write_xlsx(Master_Stability_Table, "LB_CPH_T1T3_Complete_Stability_Matrix.xlsx")


# ====Bootstrap plots for edge weights and centrality stability====
#T1T2
p_edge_T1T2 <- plot(
  adjMat.T1T2.eboot,
  statistics = "edge",
  plot = "area",
  CIstyle = "quantiles",
  order = "sample",
  labels = TRUE,  # Show all edge labels
  legend = TRUE  # Show the legend
) +
  labs(
    title = "A: T1→T2 Cross-lagged Panel Network"
  ) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(
      hjust = 0, size = 32, face = "bold"
    ),
    strip.text = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "top"
  )

# Save the TIFF
ggsave(
  filename = "adjT1T2.edgeboot0.tiff",
  plot = p_edge_T1T2,
  device = "tiff",
  width = 20,
  height = 100,
  units = "in",
  dpi = 300,
  limitsize = FALSE,
  compression = "lzw",
  bg = "white"
)

#T2T3
p_edge_T2T3 <- plot(
  adjMat.T2T3.eboot,
  statistics = "edge",
  plot = "area",
  CIstyle = "quantiles",
  order = "sample",
  labels = TRUE,  # Show all edge labels
  legend = TRUE  # Show the legend
) +
  labs(
    title = "B: T2→T3 Cross-lagged Panel Network"
  ) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(
      hjust = 0, size = 32, face = "bold"
    ),
    strip.text = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "top"
  )

# Save the TIFF
ggsave(
  filename = "adjT2T3.edgeboot0.tiff",
  plot = p_edge_T2T3,
  device = "tiff",
  width = 20,
  height = 100,
  units = "in",
  dpi = 300,
  limitsize = FALSE,
  compression = "lzw",
  bg = "white"
)
# Plot case-dropping centrality stability
nodeboot_T1T2 <- plot(adjMat.T1T2.nboot, statistics=c("outExpectedInfluence", "inExpectedInfluence")) +
  ggtitle("A: T1→T2 Cross-lagged Panel Network") +
  theme(plot.title=element_text(hjust=0,face="bold",size=16),
        axis.text = element_text(size = 8),
        axis.title = element_text(size=8),
        text = element_text(family = "serif"))



nodeboot_T2T3 <- plot(adjMat.T2T3.nboot, statistics=c("outExpectedInfluence", "inExpectedInfluence")) +
  ggtitle("B: T2→T3 Cross-lagged Panel Network") +
  theme(plot.title=element_text(hjust=0,face="bold",size=16),
        axis.text = element_text(size = 8),
        axis.title = element_text(size=8),
        text = element_text(family = "serif"))


combined_nodeboot <- nodeboot_T1T2/nodeboot_T2T3
ggsave("corStability_Combined.tiff", plot = combined_nodeboot,
       width = 18, height = 20, units = "in", dpi = 300, compression = "lzw")
#edgeboot
combined_edgeboot <- p_edge_T1T2|p_edge_T2T3
ggsave("edgeboot_Combined.tiff", plot = combined_edgeboot,
       width = 45, height = 100, units = "in", dpi = 300, limitsize = FALSE,compression = "lzw")
