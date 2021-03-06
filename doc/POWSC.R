## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
suppressMessages(library(POWSC))

## ----quick_start, eval = FALSE------------------------------------------------
#  data("es_mef_sce")
#  sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
#  estParas = Est2Phase(sce)
#  # Users can customize how many cells they need as to change the number of n.
#  simData = Simulate2SCE(n=100, estParas1 = estParas, estParas2 = estParas)
#  de = runMAST(simData$sce)
#  estPower1 = Power_Disc(de, simData = simData)
#  estPower2 = Power_Cont(de, simData = simData)

## ----multiple-group, eval = F-------------------------------------------------
#  sim_size = 1000
#  cell_per = c(0.2, 0.3, 0.5)
#  data("GSE67835")
#  col = colData(sce)
#  exprs = assays(sce)$counts
#  (tb = table(colData(sce)$Patients, colData(sce)$cellTypes))
#  # use AB_S7 patient as example and take three cell types: astrocytes hybrid and neurons
#  estParas_set = NULL
#  celltypes = c("oligodendrocytes", "hybrid", "neurons")
#  for (cp in celltypes){
#      print(cp)
#      ix = intersect(grep(cp, col$cellTypes), grep("AB_S7", col$Patients))
#      tmp_mat = exprs[, ix]
#      tmp_paras = Est2Phase(tmp_mat)
#      estParas_set[[cp]] = tmp_paras
#  }
#  #########
#  #########  Simulation part
#  #########
#  sim = SimulateMultiSCEs(n = 6000, estParas_set = estParas_set, multiProb = cell_per)
#  
#  #########
#  #########  DE analysis part
#  #########
#  DE_rslt = NULL
#  for (comp in names(sim)){
#      tmp = runMAST(sim[[comp]]$sce)
#      DE_rslt[[comp]] = tmp
#  }
#  
#  #########
#  ######### Summarize the power result
#  #########
#  pow_rslt = pow1 = pow2 = pow1_marg = pow2_marg = NULL
#  TD = CD = NULL
#  for (comp in names(sim)){
#      tmp1 = Power_Disc(DE_rslt[[comp]], sim[[comp]])
#      tmp2 = Power_Cont(DE_rslt[[comp]], sim[[comp]])
#      TD = c(TD, tmp2$TD); CD = c(CD, tmp2$CD)
#      pow1_marg = c(pow1_marg, tmp1$power.marginal)
#      pow2_marg = c(pow2_marg, tmp2$power.marginal)
#      pow_rslt[[comp]] = list(pow1 = tmp1, pow2 = tmp2)
#      pow1 = rbind(pow1, tmp1$power)
#      pow2 = rbind(pow2, tmp2$power)
#  }
#  
#  #########
#  ######### Demonstrate the result by heatmap
#  #########
#  library(RColorBrewer); library(pheatmap)
#  breaksList = seq(0, 1, by = 0.01)
#  colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
#  dimnames(pow1) = list(names(sim), names(tmp1$CD))
#  dimnames(pow2) = list(names(sim), names(tmp2$CD))
#  pheatmap(pow1, display_numbers = T, color=colors, show_rownames = T,
#           cellwidth = 50, cellheight = 50, legend = T,
#           border_color = "grey96", na_col = "grey",
#           cluster_row = FALSE, cluster_cols = FALSE,
#           breaks = seq(0, 1, 0.01),
#           main = "")
#  

## -----------------------------------------------------------------------------
sessionInfo()

