# Tcell.SNE.densVM.server.R
# Chamith Fonseka
# 04 Jan 2017
#
# Objective: Perform SNE + densVM analysis on gated B cell synovial tissue CyTOF data
#
###############################################################################
require(Rtsne)

source("custom.densVM.Fn.R")

load("synData.Tcell.downsample.RData")


markers <- colnames(synData.Tcell.downsample)[1:35]

# Run SNE
perp <- 30 # default perplexity
print(paste("Starting SNE embedding at", format(Sys.time(), "%r"), "on", format(Sys.time(), "%D")))
sne <- Rtsne(X = as.matrix(synData.Tcell.downsample[,markers]), pca = F, perplexity = perp, theta = .5, check_duplicates = F, verbose = T)
print(paste("Finished SNE embedding at", format(Sys.time(), "%r"), "on", format(Sys.time(), "%D")))
synData.Tcell.downsample <- cbind(synData.Tcell.downsample, data.frame(SNE1 = sne$Y[,1], SNE2 = sne$Y[,2]))

# Run densVM
print(paste("Starting densVM clustering at", format(Sys.time(), "%r"), "on", format(Sys.time(), "%D")))
# Make xdata expression values used for SNE and ydata SNE axes
densVM.cluster <- densVM_cluster(xdata = as.matrix(synData.Tcell.downsample[,markers]),
                                 ydata = as.matrix(synData.Tcell.downsample[c("SNE1","SNE2")]))
print(paste("Finished densVM clustering at", format(Sys.time(), "%r"), "on", format(Sys.time(), "%D")))
synData.Tcell.downsample <- cbind(synData.Tcell.downsample, data.frame(SNE.cluster = densVM.cluster$clusters$cluster))

# Save file and densVM output as RData files
save(synData.Tcell.downsample, file = "synData.Tcell.downsample.SNE.RData")
save(densVM.cluster, file = "synData.Tcell.downsample.densVM.RData")

print("Finished")
