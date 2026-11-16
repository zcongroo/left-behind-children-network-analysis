# Cross-sectional network analyses at T1, T2, and T3.
# The variable suffixes t2, t3, and t4 refer to the source dataset waves
# and correspond to T1, T2, and T3 in the study, respectively.
# Supply the prepared SPSS wave files with CSN_T1_DATA_FILE, CSN_T2_DATA_FILE,
# and CSN_T3_DATA_FILE, or use the default filenames below.
# Reports, figures, and tables are written to the working directory.

#====Load required packages====
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!require(haven)) { install.packages("haven"); library(haven) } # SPSS data import
if (!require(mgm)) { install.packages("mgm"); library(mgm) } # Nodewise predictability
if (!require(qgraph)) { install.packages("qgraph"); library(qgraph) } # Network plots and node centrality
if (!require(bootnet)) { install.packages("bootnet"); library(bootnet) } # Network estimation and bootstrap analyses
if (!require(NetworkComparisonTest)) { install.packages("NetworkComparisonTest"); library(NetworkComparisonTest) } # Paired permutation tests
if (!require(networktools)) { install.packages("networktools"); library(networktools) } # Bridge centrality
if (!require(dplyr)) { install.packages("dplyr"); library(dplyr) } # Combine and summarize tables
if (!require(tidyr)) { install.packages("tidyr"); library(tidyr) } # Reshape bridge centrality tables
if (!require(ggplot2)) { install.packages("ggplot2"); library(ggplot2) } # Plot centrality and save figures
if (!require(patchwork)) { install.packages("patchwork"); library(patchwork) } # Combine plots across waves

#====Import the analysis data====
set.seed(100)
t1_file <- Sys.getenv("CSN_T1_DATA_FILE", unset = "analysis_data_T1.sav")
t2_file <- Sys.getenv("CSN_T2_DATA_FILE", unset = "analysis_data_T2.sav")
t3_file <- Sys.getenv("CSN_T3_DATA_FILE", unset = "analysis_data_T3.sav")
if (!all(file.exists(c(t1_file, t2_file, t3_file)))) {
  stop("Supply the three prepared SPSS wave files via CSN_T1_DATA_FILE, CSN_T2_DATA_FILE, and CSN_T3_DATA_FILE.")
}

# Select columns 3 through 41 to exclude the first two columns at import.
# The first 36 selected columns are item scores in the order listed below.
# The final three selected columns are Gender, Age, and SES.
framet2 <- as.data.frame(read_sav(t1_file, col_select = 3:41))
framet3 <- as.data.frame(read_sav(t2_file, col_select = 3:41))
framet4 <- as.data.frame(read_sav(t3_file, col_select = 3:41))

# Paired tests require one row per participant in the same order at all waves.
if (nrow(framet2) != nrow(framet3) || nrow(framet2) != nrow(framet4)) {
  stop("The three wave files must have the same participants in the same row order.")
}

LBCESDCPHT2 <- framet2[1:36]
LBCESDCPHt3 <- framet3[1:36]
LBCESDCPHt4 <- framet4[1:36]

# ====Define node labels and constructs====
# Read the analysis column names for each wave
labelnamesr2 = colnames(framet2)
labelnamesr3 = colnames(framet3)
labelnamesr4 = colnames(framet4)

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

# Set the construct sizes and assign nodes to constructs
# Set number of symptoms
k = 36
k1 = 2
k2 = 20
k3 = 6
k4 = 8
# Set number of covariate
Covar <- c("Gender","Age","SES")

groups <- c(rep("Smartphone Usage Time", k1), rep("CES-D", k2),
            rep("Positive Suicide Ideation",k3), rep("Negative Suicide Ideation", k4))



# ====Residualize the nodes on covariates and estimate MGM predictability====
adjusted_t2 <- as.data.frame(matrix(nrow = nrow(framet2), ncol = k))
colnames(adjusted_t2) <- labelnamesr2[1:36]
for (i in 1:(ncol(adjusted_t2))){
  formula <- paste(labelnamesr2[i], "~Gender+Age+SES")
  lm1 <- glm(formula = formula, family = "gaussian", data = framet2)
  y_res <- lm1$residuals
  adjusted_t2[,i] <- y_res
}
networkt2 <- as.matrix(adjusted_t2)
p <- ncol(networkt2)

set.seed(100)
fit_objt2 <- mgm(data = networkt2,
                 labels = labelnames,
                 type = rep("g",p),
                 level = rep(1,p),
                 lambdaSel = "CV",
                 ruleReg = "OR")

pred_objt2 <- predict(object = fit_objt2,
                      data = networkt2,
                      errorCon = "R2")
pred_objt2$errors

weights_t2 <- fit_objt2$pairwise$wadj
signs_t2 <- fit_objt2$pairwise$signs
signs_t2[is.na(signs_t2)] <- 1
networkt2.m <- weights_t2*signs_t2
#t3
adjusted_t3 <- as.data.frame(matrix(nrow = nrow(framet3), ncol = k))
colnames(adjusted_t3) <- labelnamesr3[1:36]
for (i in 1:(ncol(adjusted_t3))){
  formula <- paste(labelnamesr3[i], "~Gender+Age+SES")
  lm1 <- glm(formula = formula, family = "gaussian", data = framet3)
  y_res <- lm1$residuals
  adjusted_t3[,i] <- y_res
}
networkt3 <- as.matrix(adjusted_t3)
p <- ncol(networkt3)

set.seed(100)
fit_objt3 <- mgm(data = networkt3,
                 labels = labelnames,
                 type = rep("g",p),
                 level = rep(1,p),
                 lambdaSel = "CV",
                 ruleReg = "OR")

pred_objt3 <- predict(object = fit_objt3,
                      data = networkt3,
                      errorCon = "R2")
pred_objt3$errors

weights_t3 <- fit_objt3$pairwise$wadj
signs_t3 <- fit_objt3$pairwise$signs
signs_t3[is.na(signs_t3)] <- 1
networkt3.m <- weights_t3*signs_t3

#t4
adjusted_t4 <- as.data.frame(matrix(nrow = nrow(framet4), ncol = k))
colnames(adjusted_t4) <- labelnamesr4[1:36]
for (i in 1:(ncol(adjusted_t4))){
  formula <- paste(labelnamesr4[i], "~Gender+Age+SES")
  lm1 <- glm(formula = formula, family = "gaussian", data = framet4)
  y_res <- lm1$residuals
  adjusted_t4[,i] <- y_res
}
networkt4 <- as.matrix(adjusted_t4)
p <- ncol(networkt4)

set.seed(100)
fit_objt4 <- mgm(data = networkt4,
                 labels = labelnames,
                 type = rep("g",p),
                 level = rep(1,p),
                 lambdaSel = "CV",
                 ruleReg = "OR")

pred_objt4 <- predict(object = fit_objt4,
                      data = networkt4,
                      errorCon = "R2")
pred_objt4$errors

weights_t4 <- fit_objt4$pairwise$wadj
signs_t4 <- fit_objt4$pairwise$signs
signs_t4[is.na(signs_t4)] <- 1
networkt4.m <- weights_t4*signs_t4


# ====Estimate the EBICglasso networks using Spearman correlations====
set.seed(100)
LBCESDCPHt2.e<-estimateNetwork(data = networkt2, labels = labelnames, default = "EBICglasso",
                               tuning = 0.5, corMethod = "cor", corArgs = list(method = "spearman",
                                                                               use = "pairwise.complete.obs"))
sink("LBCESDCPHt1e.txt")
summary(LBCESDCPHt2.e)
sink()


set.seed(100)
LBCESDCPHt3.e<-estimateNetwork(data = networkt3, labels = labelnames, default = "EBICglasso",
                               tuning = 0.5, corMethod = "cor", corArgs = list(method = "spearman",
                                                                               use = "pairwise.complete.obs"))
sink("LBCESDCPHt2e.txt")
summary(LBCESDCPHt3.e)
sink()


set.seed(100)
LBCESDCPHt4.e<-estimateNetwork(data = networkt4, labels = labelnames, default = "EBICglasso",
                               tuning = 0.5, corMethod = "cor", corArgs = list(method = "spearman",
                                                                               use = "pairwise.complete.obs"))
sink("LBCESDCPHt3e.txt")
summary(LBCESDCPHt4.e)
sink()


# ====Bootstrap edge weights and test differences in edges and node centrality====
set.seed(100)
boott1<-bootnet(LBCESDCPHt2.e,nBoots = 1000, nCores=6,
                type="nonparametric",communities = groups, directed = F,
                statistics = c("edge","strength","expectedInfluence","closeness", "betweenness"))
gc()
write.csv(boott1$sampleTable, "bootT1.sampleTable.csv")
boott1.eboot <- summary(boott1,statistics = c("edge","strength","closeness","betweenness","expectedInfluence"))
write.csv(boott1.eboot,"sum_bootT1.csv", row.names = F)

tiff(filename = "edgebootT1.tiff", width = 20,height = 100,units = "in",res = 300,compression = "lzw")
plot(boott1,labels=T,order="sample")
dev.off()

tiff(filename = "strengthDifT1.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott1, "strength", plot = "difference",order="sample")
dev.off()

tiff(filename = "closenessDifT1.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott1, "Closeness", plot = "difference",order="sample")
dev.off()

tiff(filename = "betweenDifT1.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott1, "Betweenness", plot = "difference",order="sample")
dev.off()

tiff(filename = "EIDifT1.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott1, "ExpectedInfluence", plot = "difference",order = "sample")
dev.off()


# Run the case-dropping bootstrap for node centrality
set.seed(100)
bootNt1<-bootnet(LBCESDCPHt2.e,nBoots = 1000, nCores=6,type = "case", statistics = c("strength","closeness","betweenness","expectedInfluence"))
write.csv(bootNt1$sampleTable, "nbootT1.sampleTable.csv")
gc()


tiff(filename = "NStaT1.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(bootNt1,statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence"))
dev.off()

sink("NCST1.txt")
corStability(bootNt1)
sink()

#t2
set.seed(100)
boott2<-bootnet(LBCESDCPHt3.e,nBoots = 1000, nCores=6,
                type="nonparametric",communities = groups, directed = F,
                statistics = c("edge","strength","expectedInfluence","closeness", "betweenness"))

write.csv(boott2$sampleTable, "bootT2.sampleTable.csv")
boott2.eboot <- summary(boott2,statistics = c("edge","strength","closeness","betweenness","expectedInfluence"))
write.csv(boott2.eboot,"sum_bootT2.csv", row.names = F)

tiff(filename = "edgeboott2.tiff", width = 20,height = 100,units = "in",res = 300,compression = "lzw")
plot(boott2,labels=F,order="sample")
dev.off()

tiff(filename = "strengthDift2.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott2, "strength", plot = "difference",order="sample")
dev.off()

tiff(filename = "closenessDift2.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott2, "closeness", plot = "difference",order="sample")
dev.off()

tiff(filename = "betweenDift2.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott2, "betweenness", plot = "difference",order="sample")
dev.off()

tiff(filename = "EIDift2.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott2, "expectedInfluence", plot = "difference",order="sample")
dev.off()



set.seed(100)
bootNt2<-bootnet(LBCESDCPHt3.e,nBoots = 1000, nCores=6,type = "case", statistics = c("strength","closeness","betweenness","expectedInfluence"))
write.csv(bootNt2$sampleTable, "nbootT2.sampleTable.csv")


tiff(filename = "NStat2.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(bootNt2,statistics = c("strength","closeness","betweenness","expectedInfluence"))
dev.off()

sink("NCSt2.txt")
corStability(bootNt2)
sink()


#t3
set.seed(100)
boott3<-bootnet(LBCESDCPHt4.e,nBoots = 1000, nCores=6,
                type="nonparametric",communities = groups, directed = F,
                statistics = c("edge","strength","expectedInfluence","closeness", "betweenness"))

write.csv(boott3$sampleTable, "bootT3.sampleTable.csv")
boott3.eboot <- summary(boott3,statistics = c("edge","strength","closeness","betweenness","expectedInfluence"))
write.csv(boott3.eboot,"sum_bootT3.csv", row.names = F)


tiff(filename = "edgeboott3.tiff", width = 20,height = 100,units = "in",res = 300,compression = "lzw")
plot(boott3,labels=F,order="sample")
dev.off()

tiff(filename = "strengthDift3.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott3, "strength", plot = "difference",order="sample")
dev.off()

tiff(filename = "closenessDift3.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott3, "closeness", plot = "difference",order="sample")
dev.off()

tiff(filename = "betweenDift3.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott3, "betweenness", plot = "difference",order="sample")
dev.off()

tiff(filename = "EIDift3.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(boott3, "expectedInfluence", plot = "difference",order="sample")
dev.off()



set.seed(100)
bootNt3<-bootnet(LBCESDCPHt4.e,nBoots = 1000, nCores=6,type = "case", statistics = c("strength","closeness","betweenness","expectedInfluence"))
write.csv(bootNt3$sampleTable, "nbootT3.sampleTable.csv")


tiff(filename = "NStat3.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
plot(bootNt3,statistics = c("strength","closeness","betweenness","expectedInfluence"))
dev.off()

sink("NCSt3.txt")
corStability(bootNt3)
sink()



networkt2.q <- LBCESDCPHt2.e$graph
networkt3.q <- LBCESDCPHt3.e$graph
networkt4.q <- LBCESDCPHt4.e$graph

# ====Prepare the network plots====
legend_text <- paste(labelnames, longnames, sep = ": ")
# Set node colors by construct
colors_palette <- c("#FFFF33", "#a8e6cf", "#FC8D62", "#74ADD1")
node_colors <- c(
  rep(colors_palette[1], k1),
  rep(colors_palette[2], k2),
  rep(colors_palette[3], k3),
  rep(colors_palette[4], k4)
)

# Plot the qgraph networks
tiff(filename = "netT1-T3.tiff", width = 20,height = 8,units = "in",res = 300,compression = "lzw", family = "serif")
lay_mat <- matrix(c(1, 2, 3,
                    4, 4, 4), nrow = 2, byrow = TRUE)
layout(lay_mat, heights = c(0.75, 0.25))
par(family = "serif", font.main = 2)

LBCESDCPHt2.q<-qgraph(networkt2.q, labels=labelnames,groups = groups,layout = "spring",
                      vsize=6,cut=0,sampleSize=nrow(LBCESDCPHT2),legend = F, pie = pred_objt2$errors$R2,
                      border.color="black",color=c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                      nodeNames=longnames,legend.cex=.3,GLratio = 2,
                      title = expression(bold("T1 Network")), title.cex = 2)

LBCESDCPHt3.q<-qgraph(networkt3.q,labels=labelnames,groups = groups,layout = LBCESDCPHt2.q$layout,
                      vsize=6,cut=0,sampleSize=nrow(LBCESDCPHt3), legend = F, pie = pred_objt3$errors$R2,
                      border.color="black",color=c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                      nodeNames=longnames,legend.cex=.3,GLratio = 2,
                      title = expression(bold("T2 Network")), title.cex = 2)

LBCESDCPHt4.q<-qgraph(networkt4.q,labels=labelnames,groups = groups,layout = LBCESDCPHt2.q$layout,
                      vsize=6,cut=0,sampleSize=nrow(LBCESDCPHt4), legend = F, pie = pred_objt4$errors$R2,
                      border.color="black",color=c("#a8e6cf", "#74ADD1", "#FC8D62","#FFFF33"),
                      nodeNames=longnames,legend.cex=.3,GLratio = 2,
                      title = expression(bold("T3 Network")), title.cex = 2)
par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# Node legend title
text(x = 1, y = 1.4,  # Set the vertical position of the legend title
     labels = "Nodes",
     font = 2,  # Use bold text
     cex = 2)  # Set the title size
legend("center",
       legend = legend_text,  # Node legend text
       
       # Set the legend markers
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

# ====Calculate node centrality====
LBCESDCPHt2.c <- centrality(LBCESDCPHt2.q)
LBCESDCPHt2.cen <- centrality_auto(LBCESDCPHt2.q)
# Export T1 centrality to a text file
sink("LBCESDCPHt1cen.txt")
print(LBCESDCPHt2.cen)
sink()

LBCESDCPHt3.c <- centrality(LBCESDCPHt3.q)
LBCESDCPHt3.cen <- centrality_auto(LBCESDCPHt3.q)
# Export T2 centrality to a text file
sink("LBCESDCPHt2cen.txt")
print(LBCESDCPHt3.cen)
sink()

LBCESDCPHt4.c <- centrality(LBCESDCPHt4.q)
LBCESDCPHt4.cen <- centrality_auto(LBCESDCPHt4.q)
# Export T3 centrality to a text file
sink("LBCESDCPHt3cen.txt")
print(LBCESDCPHt4.cen)
sink()

# ====Plot node centrality====
tiff(filename = "centralityT1o.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt2.q,include = c("all"))
dev.off()

tiff(filename = "centralityZT1o.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt2.q, include =c("all"), scale = "z-scores")
dev.off()

tiff(filename = "centralityT2o.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt3.q,include = c("all"))
dev.off()

tiff(filename = "centralityZT2o.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt3.q, include =c("all"), scale = "z-scores")
dev.off()

tiff(filename = "centralityT3o.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt4.q,include = c("all"))
dev.off()

tiff(filename = "centralityZT3o.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt4.q, include =c("all"), scale = "z-scores")
dev.off()

#only S&EI
tiff(filename = "centralityT1.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt2.q,include = c("Strength","ExpectedInfluence"))
dev.off()

tiff(filename = "centralityZT1.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt2.q, include =c("Strength","ExpectedInfluence"), scale = "z-scores")
dev.off()

tiff(filename = "centralityT2.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt3.q,include = c("Strength","ExpectedInfluence"))
dev.off()

tiff(filename = "centralityZT2.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt3.q, include =c("Strength","ExpectedInfluence"), scale = "z-scores")
dev.off()

tiff(filename = "centralityT3.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt4.q,include = c("Strength","ExpectedInfluence"))
dev.off()

tiff(filename = "centralityZT3.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(LBCESDCPHt4.q, include =c("Strength","ExpectedInfluence"), scale = "z-scores")
dev.off()

tiff(filename = "centralityT1-T3Z.tiff", width = 15,height = 8,units = "in",res = 300,compression = "lzw")
centralityPlot(list(T1 = LBCESDCPHt2.q, T2 = LBCESDCPHt3.q, T3 = LBCESDCPHt4.q),include = c("all"), scale = "z-scores")
dev.off()

tiff(filename = "centralityT1-T3EIZ1.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
cenPlot <- centralityPlot(list(T1 = LBCESDCPHt2.q, T2 = LBCESDCPHt3.q, T3 = LBCESDCPHt4.q),include = c("Strength","ExpectedInfluence"), scale = "z-scores")
cenPlot <- cenPlot + labs(color = "Time", linetype = "Time", shape = "Time", x = "Centrality(Z-scored)", y = "Nodes") + theme(text = element_text(family = "serif"))
print(cenPlot)
dev.off()
# Plot raw centrality values
tiff(filename = "centralityT1-T3Origin.tiff", width = 8,height = 8,units = "in",res = 300,compression = "lzw")
cenPloto <- centralityPlot(list(T1 = LBCESDCPHt2.q, T2 = LBCESDCPHt3.q, T3 = LBCESDCPHt4.q),include = c("Strength","ExpectedInfluence"), scale = "raw")
cenPloto <- cenPloto + labs(color = "Time", linetype = "Time", shape = "Time", x = "Centrality(raw)", y = "Nodes") + theme(text = element_text(family = "serif"))
print(cenPloto)
dev.off()

# ====Export centrality tables====
cenTable_T1 <- centralityTable(LBCESDCPHt2.q, standardized = T)
write.csv(cenTable_T1, "Standardized_Centrality_T1.csv")

cenTable_T2 <- centralityTable(LBCESDCPHt3.q, standardized = T)
write.csv(cenTable_T2, "Standardized_Centrality_T2.csv")

cenTable_T3 <- centralityTable(LBCESDCPHt4.q, standardized = T)
write.csv(cenTable_T3, "Standardized_Centrality_T3.csv")
# ====Compare networks across waves====
LBCESDCPHt2.NCTe <- LBCESDCPHt2.e
LBCESDCPHt2.NCTe$labels <- labelnames
colnames(LBCESDCPHt2.NCTe$data) <- labelnames


LBCESDCPHt3.NCTe <- LBCESDCPHt3.e
LBCESDCPHt3.NCTe$labels <- labelnames
colnames(LBCESDCPHt3.NCTe$data) <- labelnames

LBCESDCPHt4.NCTe <- LBCESDCPHt4.e
LBCESDCPHt4.NCTe$labels <- labelnames
colnames(LBCESDCPHt4.NCTe$data) <- labelnames

# Compare the estimated networks using paired permutation tests
#t1t2
set.seed(100)
comT1T2 <- NCT(LBCESDCPHt2.NCTe, LBCESDCPHt3.NCTe, it = 1000, binary.data = F, paired = T, weighted = T,
               AND = F, abs = T, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
               centrality = c("strength", "expectedInfluence"),  progressbar = T,
               make.positive.definite = T, verbose = T)

sink("NCTt1t2.txt")
comT1T2
sink()

set.seed(100)
comT1T2GEI <- NCT(LBCESDCPHt2.NCTe, LBCESDCPHt3.NCTe, it = 1000, binary.data = F, paired = T, weighted = T,
                  AND = F, abs = F, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
                  centrality = c("strength", "expectedInfluence"),  progressbar = T,
                  make.positive.definite = T, verbose = T)

sink("NCTt1t2GEI.txt")
comT1T2GEI
sink()
#t2t3
set.seed(100)
comT2T3 <- NCT(LBCESDCPHt3.NCTe, LBCESDCPHt4.NCTe, it = 1000, binary.data = F, paired = T, weighted = T,
               AND = F, abs = T, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
               centrality = c("strength", "expectedInfluence"),  progressbar = T,
               make.positive.definite = T, verbose = T)

sink("NCTt2t3.txt")
comT2T3
sink()

set.seed(100)
comT2T3GEI <- NCT(LBCESDCPHt3.NCTe, LBCESDCPHt4.NCTe, it = 1000, binary.data = F, paired = T, weighted = T,
                  AND = F, abs = F, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
                  centrality = c("strength", "expectedInfluence"),  progressbar = T,
                  make.positive.definite = T, verbose = T)

sink("NCTt2t3EI.txt")
comT2T3GEI
sink()
#t1t3
set.seed(100)
comT1T3 <- NCT(LBCESDCPHt2.NCTe, LBCESDCPHt4.NCTe, it = 1000, binary.data = F, paired = T, weighted = T,
               AND = F, abs = T, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
               centrality = c("strength", "expectedInfluence"),  progressbar = T,
               make.positive.definite = T, verbose = T)

sink("NCTt1t3.txt")
comT1T3
sink()

set.seed(100)
comT1T3GEI <- NCT(LBCESDCPHt2.NCTe, LBCESDCPHt4.NCTe, it = 1000, binary.data = F, paired = T, weighted = T,
                  AND = F, abs = F, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
                  centrality = c("strength", "expectedInfluence"),  progressbar = T,
                  make.positive.definite = T, verbose = T)

sink("NCTt1t3GEI.txt")
comT1T3GEI
sink()


tiff(filename = "COMT1T2globle strength dif.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT1T2, what = "strength")
dev.off()

tiff(filename = "COMT1T2 network dif.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT1T2, what = "network")
dev.off()

# Create the folders used by the multipage NCT plots.
dir.create("T1-T2centrality_dif", showWarnings = FALSE)
dir.create("T1-T2edge_dif", showWarnings = FALSE)
dir.create("T2-T3centrality_dif", showWarnings = FALSE)
dir.create("T2-T3edge_dif", showWarnings = FALSE)

tiff(filename = "T1-T2centrality_dif/COMT1T2centrality_dif_%03d.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT1T2, what = "centrality")
dev.off()

tiff(filename = "T1-T2edge_dif/COMT1T2edge_dif_%03d.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT1T2, what = "edge")
dev.off()


tiff(filename = "COMT2T3globle strength dif.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT2T3, what = "strength")
dev.off()

tiff(filename = "COMT2T3 network dif.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT2T3, what = "network")
dev.off()

tiff(filename = "T2-T3centrality_dif/COMT2T3centrality_dif_%03d.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT2T3, what = "centrality")
dev.off()

tiff(filename = "T2-T3edge_dif/COMT2T3edge_dif_%03d.tiff", width = 10,height = 10,units = "in",res = 300,compression = "lzw")
plot(comT2T3, what = "edge")
dev.off()

# Supplementary paired network comparisons using the residualized item data
covmatrix_t2 <- adjusted_t2
covmatrix_t3 <- adjusted_t3
covmatrix_t4 <- adjusted_t4
colnames(covmatrix_t2) <- labelnames
colnames(covmatrix_t3) <- labelnames
colnames(covmatrix_t4) <- labelnames


set.seed(100)
comT1T2test <- NCT(covmatrix_t2, covmatrix_t3, it = 1000, binary.data = F, paired = T, weighted = T,
                   AND = F, abs = F, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
                   centrality = c("expectedInfluence"),  progressbar = T,
                   make.positive.definite = T, verbose = T)

sink("NCTt1t2test.txt")
comT1T2test
sink()

set.seed(100)
comT2T3test <- NCT(covmatrix_t3, covmatrix_t4, it = 1000, binary.data = F, paired = T, weighted = T,
                   AND = F, abs = F, test.edges = T, edges = "all", nodes = "all", test.centrality = T,
                   centrality = c("expectedInfluence"),  progressbar = T,
                   make.positive.definite = T, verbose = T)

sink("NCTt2t3test.txt")
comT2T3test
sink()

# ====Calculate bridge centrality====
#T2
LBCESDCPHT2.b <- bridge(LBCESDCPHt2.q, communities = groups, directed = F, nodes = labelnames)

tiff(filename = "bridge_all_T2.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
plot(LBCESDCPHT2.b, zscore = T, order = "value")
dev.off()

tiff(filename = "bridge_T2.tiff", width = 8,height = 10,units = "in",res = 300,compression = "lzw")
plot(LBCESDCPHT2.b, include = c("Bridge Strength", "Bridge Expected Influence (1-step)" ), zscore = T, order = "value")
dev.off()

#T3
LBCESDCPHT3.b <- bridge(LBCESDCPHt3.q, communities = groups, directed = F, nodes = labelnames)

tiff(filename = "bridge_all_T3.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
plot(LBCESDCPHT3.b, zscore = T, order = "value")
dev.off()

tiff(filename = "bridge_T3.tiff", width = 8,height = 10,units = "in",res = 300,compression = "lzw")
plot(LBCESDCPHT3.b, include = c("Bridge Strength", "Bridge Expected Influence (1-step)" ), zscore = T, order = "value")
dev.off()

#T4
LBCESDCPHT4.b <- bridge(LBCESDCPHt4.q, communities = groups, directed = F, nodes = labelnames)

tiff(filename = "bridge_all_T4.tiff", width = 18,height = 10,units = "in",res = 300,compression = "lzw")
plot(LBCESDCPHT4.b, zscore = T, order = "value")
dev.off()

tiff(filename = "bridge_T4.tiff", width = 8,height = 10,units = "in",res = 300,compression = "lzw")
plot(LBCESDCPHT4.b, include = c("Bridge Strength", "Bridge Expected Influence (1-step)" ), zscore = T, order = "value")
dev.off()
# Extract and standardize bridge centrality
get_bridge_table <- function(bridge_obj, time_label) {
  data.frame(
    Node = names(bridge_obj$`Bridge Expected Influence (1-step)`),
    Bridge_EI = as.numeric(scale(bridge_obj$`Bridge Expected Influence (1-step)`)),
    Bridge_Strength = as.numeric(scale(bridge_obj$`Bridge Strength`)),
    Time = time_label,
    stringsAsFactors = FALSE
  )
}


bt1 <- get_bridge_table(LBCESDCPHT2.b, "T1")
bt2 <- get_bridge_table(LBCESDCPHT3.b, "T2")
bt3 <- get_bridge_table(LBCESDCPHT4.b, "T3")


bridge_all_wide <- bind_rows(bt1, bt2, bt3) %>%
  pivot_wider(
    names_from = Time,
    values_from = c(Bridge_EI, Bridge_Strength),
    names_glue = "{Time}_{.value}"
  ) %>%
  mutate(Mean_Bridge_EI = (T1_Bridge_EI + T2_Bridge_EI + T3_Bridge_EI) / 3) %>%
  arrange(desc(Mean_Bridge_EI))


write.csv(bridge_all_wide, "Standardized_Bridge_Centrality_T1_T3.csv", row.names = FALSE)

# Plot standardized bridge centrality
custom_levels <- c(
  sprintf("C%02d", 1:20),  # C01 through C20
  sprintf("H%02d", 1:2),  # H01 and H02
  sprintf("SN%02d", 1:8),  # SN01 through SN08
  sprintf("SP%02d", 1:6)  # SP01 through SP06
)

extract_bridge <- function(bridge_obj, time_label, node_names) {
  data.frame(
    Node = factor(node_names, levels = custom_levels),
    `Bridge Strength` = as.numeric(scale(bridge_obj$`Bridge Strength`)),
    `Bridge Expected Influence (1-step)` = as.numeric(scale(bridge_obj$`Bridge Expected Influence (1-step)`)),
    Time = time_label,
    check.names = FALSE
  )
}


df_t1 <- extract_bridge(LBCESDCPHT2.b, "T1", labelnames)
df_t2 <- extract_bridge(LBCESDCPHT3.b, "T2", labelnames)
df_t3 <- extract_bridge(LBCESDCPHT4.b, "T3", labelnames)

df_all <- rbind(df_t1, df_t2, df_t3)
df_long <- pivot_longer(df_all,
                        cols = c("Bridge Strength", "Bridge Expected Influence (1-step)"),
                        names_to = "Measure",
                        values_to = "Value")

df_long$Measure <- factor(df_long$Measure,
                          levels = c("Bridge Strength", "Bridge Expected Influence (1-step)"))
df_long <- df_long %>% arrange(Measure, Time, Node)


tiff(filename = "bridge_T1-T3_combined_Final.tiff", width = 8, height = 8, units = "in", res = 300, compression = "lzw")

bridge_plot <- ggplot(df_long, aes(x = Value, y = Node, group = Time, color = Time)) +
  geom_path(linewidth = 0.5, linetype = "solid") +
  geom_point(size = 1.5, shape = 16) +
  facet_grid(. ~ Measure) +
  theme_bw() +
  labs(x = "Bridge Centrality (Z-scored)", y = "Nodes", color = "Time") +
  theme(text = element_text(family = "serif"),
        axis.text.y = element_text(size = 8),
        legend.position = "right")

print(bridge_plot)
dev.off()
# Plot raw bridge centrality values
bridge_measures <- c("Bridge Strength", "Bridge Expected Influence (1-step)")

extract_bridge_raw <- function(bridge_obj, time_label, node_names) {
  data.frame(
    Node = factor(node_names, levels = custom_levels),
    `Bridge Strength` = as.numeric(bridge_obj[["Bridge Strength"]]),
    `Bridge Expected Influence (1-step)` =
      as.numeric(bridge_obj[["Bridge Expected Influence (1-step)"]]),
    Time = time_label,
    check.names = FALSE
  )
}

df_bridge_raw <- rbind(
  extract_bridge_raw(LBCESDCPHT2.b, "T1", labelnames),
  extract_bridge_raw(LBCESDCPHT3.b, "T2", labelnames),
  extract_bridge_raw(LBCESDCPHT4.b, "T3", labelnames)
)

df_bridge_raw_long <- tidyr::pivot_longer(
  df_bridge_raw,
  cols = all_of(bridge_measures),
  names_to = "Measure",
  values_to = "Value"
)
df_bridge_raw_long$Measure <- factor(df_bridge_raw_long$Measure, levels = bridge_measures)
df_bridge_raw_long <- dplyr::arrange(df_bridge_raw_long, Measure, Time, Node)

bridge_plot_raw <- ggplot(df_bridge_raw_long, aes(x = Value, y = Node, group = Time, color = Time)) +
  geom_path(linewidth = 0.5, linetype = "solid") +
  geom_point(size = 1.5, shape = 16) +
  facet_grid(. ~ Measure) +
  theme_bw() +
  labs(x = "Bridge centrality (raw)", y = "Nodes", color = "Time") +
  theme(
    text = element_text(family = "serif"),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )

ggsave(
  filename = "bridge_T1-T3_raw.tiff",
  plot = bridge_plot_raw,
  device = "tiff",
  width = 8, height = 8, units = "in",
  dpi = 300, compression = "lzw"
)

# ====Assess bridge centrality stability====
set.seed(100)

bootNt1_bridge <- bootnet(LBCESDCPHt2.e, nBoots = 1000, nCores = 6,
                          type = "case", communities = groups,
                          statistics = c("bridgeStrength", "bridgeExpectedInfluence"))

tiff(filename = "Bridge_NStaT1.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
plot(bootNt1_bridge, statistics = c("bridgeStrength", "bridgeExpectedInfluence"))
dev.off()

sink("Bridge_NCST1.txt")
corStability(bootNt1_bridge)
sink()


#t2
set.seed(100)
bootNt2_bridge <- bootnet(LBCESDCPHt3.e, nBoots = 1000, nCores = 6,
                          type = "case", communities = groups,
                          statistics = c("bridgeStrength", "bridgeExpectedInfluence"))
tiff(filename = "Bridge_NStaT2.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
plot(bootNt2_bridge, statistics = c("bridgeStrength", "bridgeExpectedInfluence"))
dev.off()
sink("Bridge_NCST2.txt")
corStability(bootNt2_bridge)
sink()

#t3
set.seed(100)
bootNt3_bridge <- bootnet(LBCESDCPHt4.e, nBoots = 1000, nCores = 6,
                          type = "case", communities = groups,
                          statistics = c("bridgeStrength", "bridgeExpectedInfluence"))


tiff(filename = "Bridge_NStaT3.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
plot(bootNt3_bridge, statistics = c("bridgeStrength", "bridgeExpectedInfluence"))
dev.off()

sink("Bridge_NCST3.txt")
corStability(bootNt3_bridge)
sink()

# ====Test differences in bridge centrality====
set.seed(100)
boott1_bridge <- bootnet(LBCESDCPHt2.e, nBoots = 1000, nCores = 6,
                         type = "nonparametric", communities = groups,
                         statistics = c("bridgeStrength", "bridgeExpectedInfluence"))

tiff(filename = "Bridge_EIDifT1.tiff", width = 15, height = 8, units = "in", res = 300, compression = "lzw")
plot(boott1_bridge, "bridgeExpectedInfluence", plot = "difference", order="sample")
dev.off()

tiff(filename = "Bridge_SDifT1.tiff", width = 15, height = 8, units = "in", res = 300, compression = "lzw")
plot(boott1_bridge, "bridgeStrength", plot = "difference", order="sample")
dev.off()

#t2
set.seed(100)
boott2_bridge <- bootnet(LBCESDCPHt3.e, nBoots = 1000, nCores = 6,
                         type = "nonparametric", communities = groups,
                         statistics = c("bridgeStrength", "bridgeExpectedInfluence"))

tiff(filename = "Bridge_EIDifT2.tiff", width = 15, height = 8, units = "in", res = 300, compression = "lzw")
plot(boott2_bridge, "bridgeExpectedInfluence", plot = "difference", order="sample")
dev.off()

tiff(filename = "Bridge_SDifT2.tiff", width = 15, height = 8, units = "in", res = 300, compression = "lzw")
plot(boott2_bridge, "bridgeStrength", plot = "difference", order="sample")
dev.off()

#t3
set.seed(100)
boott3_bridge <- bootnet(LBCESDCPHt4.e, nBoots = 1000, nCores = 6,
                         type = "nonparametric", communities = groups,
                         statistics = c("bridgeStrength", "bridgeExpectedInfluence"))

tiff(filename = "Bridge_EIDifT3.tiff", width = 15, height = 8, units = "in", res = 300, compression = "lzw")
plot(boott3_bridge, "bridgeExpectedInfluence", plot = "difference", order="sample")
dev.off()

tiff(filename = "Bridge_SDifT3.tiff", width = 15, height = 8, units = "in", res = 300, compression = "lzw")
plot(boott3_bridge, "bridgeStrength", plot = "difference", order="sample")
dev.off()

# Plot pairwise differences between edge weights. This step can take substantial time.
tiff(filename = "edgeDifT1.tiff", width = 18,height = 15,units = "in",res = 300,compression = "lzw")
plot(boott1, "edge", plot = "difference",order = "sample", onlyNonZero = T)
dev.off()

tiff(filename = "edgeDift2.tiff", width = 18,height = 15,units = "in",res = 300,compression = "lzw")
plot(boott2, "edge", plot = "difference",order = "sample", onlyNonZero = T)
dev.off()

tiff(filename = "edgeDift3.tiff", width = 18,height = 15,units = "in",res = 300,compression = "lzw")
plot(boott3, "edge", plot = "difference",order = "sample", onlyNonZero = T)
dev.off()

# Combine edge-weight bootstrap plots
edgetheme <- theme(
  plot.title = element_text(
    family = "serif",       #
    face = "bold",  # Use bold text
    size = 32,  # Set the title size
    hjust = 0  # Align the title to the left
  )
)
T1edgeboot <- plot(boott1,labels=T,order="sample") + ggtitle("A: T1 Network") + edgetheme
T2edgeboot <- plot(boott2,labels=T,order="sample") + ggtitle("B: T2 Network") + edgetheme
T3edgeboot <- plot(boott3,labels=T,order="sample") + ggtitle("C: T3 Network") + edgetheme

combined_plot <- T1edgeboot|T2edgeboot|T3edgeboot
ggsave("Bootstrap_Edges_Combined.tiff", plot = combined_plot,
       width = 70, height = 100, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
# Combine case-dropping stability plots
nodetheme <- theme(
  plot.title = element_text(
    family = "serif",       #
    face = "bold",  # Use bold text
    size = 24,              #
    hjust = 0  # Align the title to the left
  )
)
T1nodeboot <- plot(bootNt1,statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence")) + ggtitle("A: T1 Network") + nodetheme
T2nodeboot <- plot(bootNt2,statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence")) + ggtitle("B: T2 Network") + nodetheme
T3nodeboot <- plot(bootNt3,statistics = c("Strength","Closeness","Betweenness","ExpectedInfluence")) + ggtitle("C: T3 Network") + nodetheme

combinednode_plot <- T1nodeboot/T2nodeboot/T3nodeboot
ggsave("Bootstrap_Nodes_Combined.tiff", plot = combinednode_plot,
       width = 15, height = 30, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
# Export pairwise centrality difference plots
ST1 <- plot(boott1, "strength", plot = "difference",order="sample")
ggsave("strengthDifT1.tiff", plot = ST1,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
ST2 <- plot(boott2, "strength", plot = "difference",order="sample")
ggsave("strengthDifT2.tiff", plot = ST2,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")

ST3 <- plot(boott3, "strength", plot = "difference",order="sample")
ggsave("strengthDifT3.tiff", plot = ST3,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")

EIT1 <- plot(boott1, "ExpectedInfluence", plot = "difference",order="sample")
ggsave("EIDifT1.tiff", plot = EIT1,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
EIT2 <- plot(boott2, "ExpectedInfluence", plot = "difference",order="sample")
ggsave("EIDifT2.tiff", plot = EIT2,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")

EIT3 <- plot(boott3, "ExpectedInfluence", plot = "difference",order="sample")
ggsave("EIDifT3.tiff", plot = EIT3,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")

BEIT1 <- plot(boott1_bridge, "bridgeExpectedInfluence", plot = "difference",order="sample")
ggsave("Bridge_EIDifT1.tiff", plot = BEIT1,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
BEIT2 <- plot(boott2_bridge, "bridgeExpectedInfluence", plot = "difference",order="sample")
ggsave("Bridge_EIDifT2.tiff", plot = BEIT2,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
BEIT3 <- plot(boott3_bridge, "bridgeExpectedInfluence", plot = "difference",order="sample")
ggsave("Bridge_EIDifT3.tiff", plot = BEIT3,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")

BST1 <- plot(boott1_bridge, "bridgeStrength", plot = "difference",order="sample")
ggsave("Bridge_SDifT1.tiff", plot = BST1,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
BST2 <- plot(boott2_bridge, "bridgeStrength", plot = "difference",order="sample")
ggsave("Bridge_SDifT2.tiff", plot = BST2,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
BST3 <- plot(boott3_bridge, "bridgeStrength", plot = "difference",order="sample")
ggsave("Bridge_SDifT3.tiff", plot = BST3,
       width = 15, height = 15, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
# Combine bridge centrality stability plots
Bnodetheme <- theme(
  plot.title = element_text(
    family = "serif",       #
    face = "bold",  # Use bold text
    size = 24,  # Set the title size
    hjust = 0               #
  )
)
T1Bnodeboot <- plot(bootNt1_bridge, statistics = c("bridgeStrength", "bridgeExpectedInfluence")) + ggtitle("A: T1 Network") + Bnodetheme
T2Bnodeboot <- plot(bootNt2_bridge, statistics = c("bridgeStrength", "bridgeExpectedInfluence")) + ggtitle("B: T2 Network") + Bnodetheme
T3Bnodeboot <- plot(bootNt3_bridge, statistics = c("bridgeStrength", "bridgeExpectedInfluence")) + ggtitle("C: T3 Network") + Bnodetheme

combinedBnode_plot <- T1Bnodeboot/T2Bnodeboot/T3Bnodeboot
ggsave("Bootstrap_BNodes_Combined.tiff", plot = combinedBnode_plot,
       width = 15, height = 30, units = "in", dpi = 300, limitsize = FALSE, compression = "lzw")
