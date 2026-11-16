# Correlation analysis at T1, T2, and T3.
# Supply the prepared SPSS analysis file by setting CLPN_DATA_FILE.
# If CLPN_DATA_FILE is unset, place analysis_data.sav in the working directory.
# The three TIFF heatmaps are written to the working directory.

#====Load required packages====
if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

if (!require(haven)) { install.packages("haven"); library(haven) } # SPSS data import
if (!require(qgraph)) { install.packages("qgraph"); library(qgraph) } # Automatic correlation estimation
if (!require(psych)) { install.packages("psych"); library(psych) } # Polychoric correlations and p-values
if (!require(corrplot)) { install.packages("corrplot"); library(corrplot) } # Correlation heatmaps

#====Import the analysis data====
data_file <- Sys.getenv("CLPN_DATA_FILE", unset = "analysis_data.sav")
if (!file.exists(data_file)) {
  stop("Set CLPN_DATA_FILE to the path of the prepared SPSS analysis file.")
}

# Read the 108 item columns in wave order from the prepared analysis file.
# In each wave: 2 smartphone usage items, 20 depression items, 6 positive
# suicidal ideation items, and 8 negative suicidal ideation items.
data1 <- read_sav(data_file, col_select = 3:110)
LBCPHT1 <- data1[, 1:36]
LBCPHT2 <- data1[, 37:72]
LBCPHT3 <- data1[, 73:108]

# Node labels in the order used at each wave.
labelnames = c("H01","H02",
               "C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20",
               "SP01","SP02","SP03","SP04","SP05","SP06",
               "SN01","SN02","SN03","SN04","SN05","SN06","SN07","SN08")

# ====Correlation analysis and heatmaps====

# Optional example for computing a correlation matrix: cor1 <- cor(...)
# Assign node labels
LBCPHT1.r <- LBCPHT1
LBCPHT2.r <- LBCPHT2
LBCPHT3.r <- LBCPHT3

colnames(LBCPHT1.r) <-labelnames
colnames(LBCPHT2.r) <-labelnames
colnames(LBCPHT3.r) <-labelnames

LBCPHT1.cor <- cor_auto(LBCPHT1.r)
LBCPHT2.cor <- cor_auto(LBCPHT2.r)
LBCPHT3.cor <- cor_auto(LBCPHT3.r)

# Calculate correlation coefficients and p-values
testRes1 = polychoric(LBCPHT1.r)
testRes1.m <- testRes1$rho
testRes1.p <- corr.p(testRes1.m,n = nrow(LBCPHT1.r))  # Calculate p-values

testRes2 = polychoric(LBCPHT2.r)
testRes2.m <- testRes2$rho
testRes2.p <- corr.p(testRes2.m,n = nrow(LBCPHT2.r))

testRes3 = polychoric(LBCPHT3.r)
testRes3.m <- testRes3$rho
testRes3.p <- corr.p(testRes3.m,n = nrow(LBCPHT3.r))


tiff(filename = "LBCPHT1cor.tiff",width = 20,height = 20,
     units = "in",res = 300,compression = "lzw")
# Plot the upper triangle
cor.plot <- corrplot(corr = LBCPHT1.cor,p.mat = testRes1.p$p,type = "upper",
                     tl.pos = "lt",tl.col = "black",  # Set label position and color
                     insig = "label_sig",sig.level = c(.001,.01, .05),  # Mark statistical significance at the specified thresholds
                     pch.cex = 1, pch.col = "black", order = "original")  # Set significance marker size and color
# Plot the lower triangle
cor.plot <- corrplot(corr = testRes1.m,type = "lower", add = T,
                     method = "number",
                     tl.pos = "n",tl.col = "black",tl.cex = 1.2,
                     col = "black",diag = F,cl.pos = "n",
                     pch.col = "black",number.cex = 0.8,
                     order = "original",number.digits = 3)
dev.off()

tiff(filename = "LBCPHT2cor.tiff",width = 20,height = 20,units = "in",res = 300,compression = "lzw")
cor.plot <- corrplot(corr = LBCPHT2.cor,p.mat = testRes2.p$p,type = "upper",
                     tl.pos = "lt",tl.col = "black",  # Set label position and color
                     insig = "label_sig",sig.level = c(.001,.01, .05),  # Mark statistical significance at the specified thresholds
                     pch.cex = 1, pch.col = "black", order = "original")  # Set significance marker size and color
cor.plot <- corrplot(corr = LBCPHT2.cor,type = "lower", add = T,
                     method = "number",
                     tl.pos = "n",tl.col = "black",tl.cex = 1.2,
                     col = "black",diag = F,cl.pos = "n",
                     pch.col = "black",number.cex = 0.8,
                     order = "original",number.digits = 3)
dev.off()

tiff(filename = "LBCPHT3cor.tiff",width = 20,height = 20,units = "in",res = 300,compression = "lzw")
cor.plot <- corrplot(corr = LBCPHT3.cor,p.mat = testRes3.p$p,type = "upper",
                     tl.pos = "lt",tl.col = "black",  # Set label position and color
                     insig = "label_sig",sig.level = c(.001,.01, .05),  # Mark statistical significance at the specified thresholds
                     pch.cex = 1, pch.col = "black", order = "original")  # Set significance marker size and color
cor.plot <- corrplot(corr = LBCPHT3.cor,type = "lower", add = T,
                     method = "number",
                     tl.pos = "n",tl.col = "black",tl.cex = 1.2,
                     col = "black",diag = F,cl.pos = "n",
                     pch.col = "black",number.cex = 0.8,
                     order = "original",number.digits = 3)
dev.off()
