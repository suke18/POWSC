#' Estimate characterized parameters for a given scRNA-seq data (SingleCellExperiment object or a count matrix).
#'
#' These parameters include four gene-wise parameters and two cell-wise parameters.
#'
#' @param sim_size a list of numbers
#' @param per_DE the precentage of the DE genes
#' @param est_Paras the template parameter estimated from one cell type
#' @param DE_Method is a string chosen from "MAST" or "SC2P".
#' @param Cell_Type is a string corresponding to the 1st scenario: same cell type comparison, and 2nd scenario: multiple cell types.
#' @param multi_Prob is the mixture cell proportions which sum up to 1. If not summing up to 1, then the package will internally do the normalization procedure.
#' @param alpha is the cutoff for the fdr which can be modified
#' @param disc_delta or the zero ratio change is the cutoff (=0.1) used to determined the high DE genes for Form II
#' @param cont_delta or the lfc is the cutoff (=0.5) used to determined the high DE genes for Form II
#' @return POWSC object
#' @examples
#' sim_size = c(100, 400, 1000) # A numeric vector
#' pow_rslt = runPOWSC(sim_size = sim_size, est_Paras = est_Paras,per_DE=0.05, DE_Method = "MAST", Cell_Type = "PW")
#' @export runPOWSC
######## POWSC ########
runPOWSC = function(sim_size = c(50, 100, 200, 800, 1000), per_DE = 0.05, est_Paras, DE_Method = c("MAST", "SC2P"),
                    Cell_Type = c("PW", "Multi"), multi_Prob = NULL,
                    alpha = 0.1, disc_delta = 0.1, cont_delta = 0.5){
    DE_Method = match.arg(DE_Method)
    Cell_Type = match.arg(Cell_Type)
    pow_rslt = NULL
    # Pair-wise comparison
    if (Cell_Type == "PW"){
        for (tmp_size in sim_size){
            sim_Data = Simulate2SCE(n = tmp_size, estParas1 = est_Paras, estParas2 = est_Paras)
            DE_rslt = runDE(sim_Data$sce, DE_Method = DE_Method)
            pow1_rslt = Power_Disc(DErslt = DE_rslt, simData = sim_Data, alpha = alpha, delta = disc_delta)
            pow2_rslt = Power_Cont(DErslt = DE_rslt, simData = sim_Data, alpha = alpha, delta = cont_delta)
            tmp_rslt = list(pow1 = pow1_rslt, pow2 = pow2_rslt)
            pow_rslt[[toString(tmp_size)]] = tmp_rslt
        }
    }else{ # Multiple Cell Types
        if (missing(multi_Prob)){
            stop("Please Estimate the Cell Proportions by SC3 or Seurat")
        }
        for (tmp_size in sim_size){
            pow1_rslt = pow2_rslt = NULL
            sim_Data = SimulateMultiSCEs(n = tmp_size, estParas_set = est_Paras, multiProb = multi_Prob)
            for (comp in names(sim_Data)){
                tmp_DE = runDE(sim_Data[[comp]]$sce, DE_Method = DE_Method)
                pow1_rslt[[comp]] = Power_Disc(DErslt = tmp_DE, simData = sim_Data[[comp]], alpha = alpha, delta = disc_delta)
                pow2_rslt[[comp]] = Power_Cont(DErslt = tmp_DE, simData = sim_Data[[comp]], alpha = alpha, delta = cont_delta)
            }
            tmp_rslt = list(pow1 = pow1_rslt, pow2 = pow2_rslt)
            pow_rslt[[toString(tmp_size)]] = tmp_rslt
        }
    }
    class(pow_rslt) = "POWSC"
    return(pow_rslt)
}

######## Plot ########
plot.POWSC = function(POWSCobj, Form = c("I", "II"), Cell_Type = c("PW", "Multi")){
    Cell_Type = match.arg(Cell_Type)
    Form = match.arg(Form)
    if (Cell_Type == "Multi"){
        colors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
        sim_size = names(POWSCobj)
        for (tmp_size in sim_size){
            tmp_pow_rslt = POWSCobj[[tmp_size]]
            pow1_mat = do.call(rbind, lapply(tmp_pow_rslt$pow1, function(x) x$power))
            pow2_mat = do.call(rbind, lapply(tmp_pow_rslt$pow2, function(x) x$power))
            colnames(pow1_mat) = names(POWSCobj[[1]][[1]][[1]][[1]])
            colnames(pow2_mat) = names(POWSCobj[[1]][[2]][[1]][[1]])
            if (Form == "I"){
                pow = pow1_mat
            }else{
                pow = pow2_mat
            }

            pheatmap(pow,display_numbers = TRUE, color=colors, show_rownames = FALSE,
                     show_colnames = TRUE, fontsize = 20, fontsize_col = 20,
                     cellwidth = 40, cellheight = 40, legend = FALSE,
                     border_color = "grey96", na_col = "grey",
                     cluster_row = FALSE, cluster_cols = FALSE,
                     breaks = seq(0, 1, 0.01),
                     main = paste0("Total Cell Number = ", tmp_size))
        }
    }
    else{
        pow1_mat = do.call(rbind, lapply(POWSCobj, function(x) x$pow1$power))
        pow2_mat = do.call(rbind, lapply(POWSCobj, function(x) x$pow2$power))
        nrep = rownames(pow1_mat)
        nm1 = names(POWSCobj[[1]][[1]][[1]])
        nm2 = names(POWSCobj[[1]][[2]][[1]])
        pow1 = data.frame(Strata = rep(nm1, each = nrow(pow1_mat)), Power = as.vector(pow1_mat), Reps = rep(nrep, ncol(pow1_mat)), stringsAsFactors = FALSE)
        pow2 = data.frame(Strata = rep(nm2, each = nrow(pow2_mat)), Power = as.vector(pow2_mat), Reps = rep(nrep, ncol(pow2_mat)), stringsAsFactors = FALSE)
        pow1$Strata = factor(pow1$Strata, levels = nm1); pow2$Strata = factor(pow2$Strata, levels = nm2)
        pow1$Reps = factor(pow1$Reps, levels = unique(pow1$Reps)); pow2$Reps = factor(pow2$Reps, levels = unique(pow2$Reps))
        if (Form == "I"){
            pow = pow1
            tit = ggtitle("Form I DE genes")
            tmpxlab = xlab("strata of zero ratios")
        }else{
            pow = pow2
            tit = ggtitle("Form II DE genes")
            tmpxlab = xlab("strata of average reads")
        }
        breaks = round(seq(0, 1, length = 6), 1)
        ggplot(pow, aes(x=Strata, y=Power, group=Reps, color=Reps)) +
            geom_line(aes(color=Reps, linetype = Reps), size = 0.6)+
            geom_point(aes(color=Reps, shape = Reps), size = 1.5) +
            scale_colour_brewer(palette = "Set1") +
            tmpxlab + tit + theme_classic()  +
            scale_y_continuous(breaks=breaks, limits = c(0, 1+0.1), labels=format(breaks, nsmall = 1))+
            theme(legend.position = c(0.5, 0.95), legend.direction = "horizontal") +
            theme(plot.title = element_text(size=16, face="bold"),
                  axis.text.x=element_text(colour="black", size = 14),
                  axis.text.y=element_text(colour="black", size = 14),
                  legend.title = element_text(size = 15),
                  legend.text = element_text(size = 15),
                  axis.title.y = element_text(size = 15))
    }
}

######## Summary Table ########
summary.POWSC = function(POWSCobj, Form = c("I", "II"), Cell_Type = c("PW", "Multi")){
    Cell_Type = match.arg(Cell_Type)
    Form = match.arg(Form)
    if (Cell_Type == "Multi"){
        pow_tab = NULL
        sim_size = names(POWSCobj)
        for (tmp_size in sim_size){
            tmp_pow_rslt = POWSCobj[[tmp_size]]
            pow1_mat = do.call(rbind, lapply(tmp_pow_rslt$pow1, function(x) x$power))
            pow2_mat = do.call(rbind, lapply(tmp_pow_rslt$pow2, function(x) x$power))
            colnames(pow1_mat) = names(POWSCobj[[1]][[1]][[1]][[1]])
            colnames(pow2_mat) = names(POWSCobj[[1]][[2]][[1]][[1]])
            if (Form == "I"){
                pow = pow1_mat
            }else{
                pow = pow2_mat
            }
            pow_tab[[tmp_size]] = pow
        }
        return(pow_tab)
    }else{
        pow1_mat = round(do.call(rbind, lapply(POWSCobj, function(x) x$pow1$power)), 4)
        pow2_mat = round(do.call(rbind, lapply(POWSCobj, function(x) x$pow2$power)), 4)
        nm1 = names(POWSCobj[[1]][[1]][[1]])
        nm2 = names(POWSCobj[[1]][[2]][[1]])
        colnames(pow1_mat) = nm1
        colnames(pow2_mat) = nm2
        if (Form == "I"){
            tab = pow1_mat
        }else{
            tab = pow2_mat
        }
        return(tab)
    }
}













