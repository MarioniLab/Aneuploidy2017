# PACKAGES ####


#load packages used for analysis
library(scran)
library(biomaRt)
library(reshape2)
library(edgeR)
library(RColorBrewer)
library(Matrix)
library(knitr)
library(BiocStyle)
library(parallel)
library(scploid)
library(ggplot2)
#configure a mart for human and mouse
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = ensembl)

mouse_ensembl = useMart("ensembl")
mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)


# UTILITY FUNCTIONS ####

# Multiple plot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' scale_colour_Publication is a set of colours that don't look too bad
#' use with ggplot like other scales
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour", "Publication",
                 manual_pal(values = c(
                   "#000000", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                   "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                   "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                   "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                   "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                   "#372101", "#FFB500", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                   "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                   "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                   "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                   "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                   "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                   "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                   "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")), ...)
}

#' rowMedians is a small wrapper for, unsurprisingly, getting row medians of matrices
#'
#' @param mat input matrix
#'
#' @return vector of row medians
#' @export
#'
#' @examples
rowMedians = function(mat){
  return(apply(mat, 1, median))
}


# SIMULATION FUNCTIONS ####

#' od_est gets an estimate of the overdispersion of data using
#' CV^2 - 1/mean. If this estimate is negative, 0 is used instead.
#'
#' @param vec Vector of data
#'
#' @return overdispersion estimate
#' @export
#'
#' @examples
od_est = function(vec){
  cv = sd(vec)/mean(vec)
  return(max( ( cv^2 - (1/mean(vec)) ), 0))
}

#' get_fit returns a loess model of the log(overdispersion) ~ log(mean)
#' relationship for all genes provided. It excludes genes with values of 0
#' in either parameter and those with very low overdispersion (<1e-5)
#'
#' @param counts_matrix is the counts matrix to make the model from
#'
#' @return loess model object
#' @export
#'
#' @examples
get_fit = function(counts_matrix){
  od = apply(counts_matrix, 1, od_est)
  mean = apply(counts_matrix, 1, mean)
  
  y = log10(od)
  x = log10(mean)
  
  #these genes really mess up the fit - remove them as they are not meaningful for our CPM>50 considerations
  remove = which(is.na(y) | is.na(x) | is.infinite(y) | is.infinite(x) | y < -5)
  y = y[-remove]
  x = x[-remove]
  
  fit = loess(y~x)
  return(fit)
}


#' subsample returns a simulated counts matrix based off of real data.
#'
#' @param counts counts matrix provided
#' @param downsample_frac the fraction of reads of the original matrix the simulated one will have
#' @param overdisperse_factor the factor by which overdispersion should be increased
#' compared to the loess model based on real data
#' @param const_mu if T, generates a NB model for each gene using mean counts fraction for all
#' cells provided, and overdispersion from the loess estimate * overdisperse_factor. If F,
#' generates a new NB model for each count in each cell using mean of the count fraction
#' for the gene in that cell and overdispersion from loess estimate * overdisperse_factor
#' @param fit if provided, overrides the calculates loess fit
#'
#' @return new simulated counts matrix
#' @export
#'
#' @examples
subsample = function(counts, downsample_frac = 1, overdisperse_factor = 1, const_mu = F, fit = NULL){
  if(is.null(fit))
    fit = get_fit(counts)
  
  lib.sizes = colSums(counts)
  
  phi = 10^predict(object = fit, newdata = log10(data.frame(x=rowMeans(counts))))
  #some genes cannot be fitted, as we excluded ultra-low expression genes from the fit
  #we therefore give them a log10(OD) (phi) that is predicted from the smallest gene that is predicted
  #this gives an approximately accurate value, and helps stop the production of NAs. This is not of huge consequence
  # for the simulations, as these genes will never be considered in our CPM>50 cutoff anyway
  phi[is.na(phi)] = 10^predict(object = fit, newdata = data.frame(x=log10(min(rowMeans(counts)[!is.na(phi)]))))
  
  new_matrix = matrix(NA, ncol = ncol(counts), nrow = nrow(counts),
                      dimnames = list(rownames(counts), colnames(counts)))
  
  lib_size = apply(counts, 2, sum)
  
  #if a constant mu, we take a mean count fraction across all cells for each gene and sample around this
  if(const_mu)
    means = rowMeans(counts)

  # old - slow
  # for(cell in 1:ncol(counts)){
  #   for(gene in 1:nrow(counts)){
  #     #if !constant mu, we take the vector of m_gi for this cell in particular
  #     if(!const_mu)
  #       mu = counts[,cell]
  #     #mean fraction is taken from the mean fraction across all cells
  #     #i.e. one NB distribution per gene
  #     new_matrix[gene, cell] = rnbinom(n = 1, mu = mu[gene] * downsample_frac,
  #                                        size = 1/(phi[gene]*overdisperse_factor))
  #   }
  # }
  
  if(const_mu){
    lib.scale = lib.sizes/mean(lib.sizes)
  } else {
    lib.scale = 1
  }

  
  for(gene in 1:nrow(counts)){
    
    if(!const_mu){
      #match the counts with the mean
      mu = counts[gene,]
    } else {
      #use the average count over all cells
      mu = means[gene]
    }
    new_matrix[gene,] = rnbinom(n = ncol(counts), mu = mu * downsample_frac * lib.scale,
                                size = 1/(phi[gene]*overdisperse_factor * lib.scale))
  }
  
  
  
  return(new_matrix)
}

#' subsample_simulate_aneuploidy also created a simulated matrix of aneuploid counts, but instead uses
#' fixed expression means of 0.5 or 1.5x the gene average given a data frame of aneuploidies
#'
#' @param counts counts matrix
#' @param downsample_frac fraction to multiply library sizes by
#' @param overdisperse_factor fraction to change NB dispersion by
#' @param gene_map map of genes (for setting of 1.5 or 0.5x expression mean by chromosome)
#' @param known_aneuploidy data frame of known aneuploidies (cell, chr, expr_down)
#'
#' @return simulated count matrix
#' @export
#'
#' @examples
subsample_simulate_aneuploidy = function(counts, downsample_frac = 1, overdisperse_factor = 1, gene_map = NULL, known_aneuploidy = NULL){
  if(is.null(gene_map) | is.null(known_aneuploidy))
    stop("Specify gene_map and known_aneuploidy")
  
  #get the fit, allocate memory, get library sizes
  fit = get_fit(counts)
  new_matrix = matrix(NA, ncol = ncol(counts), nrow = nrow(counts),
                      dimnames = list(rownames(counts), colnames(counts)))
  libs = apply(counts, 2, sum)
  
  
  #make a ploidy state matrix
  ploidy_mat = matrix(1, ncol = ncol(counts), nrow = nrow(counts),
                      dimnames = list(rownames(counts), colnames(counts)))
  for(i in 1:nrow(known_aneuploidy)){
    chr = as.numeric(as.character(known_aneuploidy$chr[i]))
    cell = as.character(known_aneuploidy$cell[i])
    factor_aneu = ifelse(known_aneuploidy$monosomy[i], 0.5, 1.5)
    
    target_genes = rownames(gene_map)[which(gene_map$chr==chr)]
    
    ploidy_mat[which(rownames(ploidy_mat)%in%target_genes), which(colnames(ploidy_mat)==cell)] = factor_aneu
  }  
  
  #estimate mu from the non-aneuploid cells
  normal_cells = counts[, !colnames(counts)%in%known_aneuploidy$cell ]
  mu = apply(normal_cells, 1, mean)
  #make a matrix of mus
  mu_matrix = matrix(mu, byrow = FALSE, ncol = ncol(counts), nrow = nrow(counts),
                     dimnames = list(rownames(counts), colnames(counts)))
  
  #correct mus to maintain library size differences
  lib_ratio = libs/mean(libs)
  mu_matrix = sweep(mu_matrix, 2, lib_ratio, "*")
  
  #adjust the mean expression by the ploidy matrix
  mu_matrix = mu_matrix * ploidy_mat
  


  #make the phi matrix, which will be the same dimensions and 
  #dimnames as the mu matrix
  phi_matrix = mu_matrix
  
  #calculate phis from the fit of the adjusted mus
  for(col in 1:ncol(phi_matrix)){
    #reassign columns of this matrix as we progress through it
    phi_matrix[, col] = 10^predict(object = fit, newdata = log10(data.frame(x = phi_matrix[, col])))
  }
  
  genes_with_na = apply(phi_matrix, 1, function(x) any(is.na(x)))
  #where phi is NA, we want instead to give it the value of the minimum considered mean expression in the fit
  phi_matrix[is.na(phi_matrix)] = 10^predict(object = fit, newdata = data.frame(x=log10(min(rowMeans(counts)[!genes_with_na]))))
  
  
  #sample from the matrix of mus, inferring phi from our fit for each mu separately
  for(cell in 1:ncol(counts)){
      
      new_matrix[, cell] = rnbinom(n = nrow(mu_matrix), mu = mu_matrix[, cell] * downsample_frac, 
                                       size = 1/(phi_matrix[, cell] * overdisperse_factor) )
      
  }
  return(new_matrix)
}


#' plot_test is a confusingly named function that simply plots real vs. simulated data
#' for comparisons of whether simulated results accurately represent the original
#' copy.
#'
#' @param real_counts is the counts matrix of the real data for comparison
#' @param sim_counts is the counts matrix of the simulated data for comparison
#' @param cpm_normalise is boolean T/F for whether to plot the CPM values (T) or
#' the raw counts (F)
#' @param limit_genes is boolean for whether to lower the number of genes
#' in the simulated data randomly to the number in the real data if the former
#' is larger.
#' @param plot if T, plots, else only returns the simulated reads, normalised by cpm 
#' and with randomly removed genes if chosen (limit_genes = T)
#'
#' @return the simulated CPM (with gene removal if requested)
#' @export
#'
#' @examples
plot_test = function(real_counts, sim_counts, cpm_normalise = T, limit_genes = T, plot = T){
  if(cpm_normalise){
    real_counts = sweep(real_counts, 2, colSums(real_counts), "/") * 1e6
    sim_counts = sweep(sim_counts, 2, colSums(sim_counts), "/") * 1e6
  }
  
  include = apply(real_counts, 1, median)>50
  real_mean = rowMeans(real_counts)[include]
  real_sd = apply(real_counts, 1, sd)[include]
  
  include = apply(sim_counts, 1, median)>50
  sim_mean = rowMeans(sim_counts)[include]
  sim_sd = apply(sim_counts, 1, sd)[include]
  
  all_mean = c(real_mean, sim_mean)
  all_sd = c(real_sd, sim_sd)
  
  if(plot){
    plot(x=real_mean, y = real_sd, log = "xy", pch = ".", col = "cornflowerblue",
         xlim = c(min(all_mean, na.rm = T), max(all_mean, na.rm = T)), ylim = c(min(all_sd, na.rm = T), max(all_sd, na.rm = T)),
         xlab = "Gene mean expression", ylab = "Gene expression S.D.")
    points(x=sim_mean, y = sim_sd, pch = ".", col = "coral")
    
    real_mod = lm(log10(real_sd) ~ log10(real_mean))
    abline(a=real_mod$coefficients[1], b=real_mod$coefficients[2], col = "blue")
    sim_mod = lm(log10(sim_sd) ~ log10(sim_mean))
    abline(a=sim_mod$coefficients[1], b=sim_mod$coefficients[2], col = "red")
    
    legend(x=2000, y = 100, legend = c("Real data", "Simulated data"), fill = c("cornflowerblue", "coral"))
  }
  if(length(sim_mean)>length(real_mean) & limit_genes){
  #we want to downsample the number of **ACTIVE** genes to match the real data
    
    active_genes = rownames(sim_counts[rowMedians(sim_counts)>50, ])
    
    active_to_keep = sample(active_genes, length(real_mean))
    active_to_bin = active_genes[!active_genes %in% active_to_keep]
    
    
    sim_counts = sim_counts[-which(rownames(sim_counts)%in%active_to_bin),]
  }
  
  return(sim_counts)
}

#' subsample_repeat_wrapper wraps the guts of subsample_repeat_analysis to allow
#' easy parallelisation
#'
#' @param raw_data_control raw counts (no aneu, using fixed mu)
#' @param raw_data_aneu raw counts (with aneu)
#' @param downsample downsample fraction
#' @param overdisperse overdispersion fraction
#' @param meta_in metadata for gt_test
#' @param gene_table_in gene table for gt_test
#' @param known_hits real aneuploidies for gt_test
#' @param nruns times to repeat analysis
#'
#' @return data frame with repeated metrics of success for each run
#' @export
#'
#' @examples
subsample_repeat_wrapper = function(raw_data_control = raw_fine,
                                    raw_data_aneu = raw_aneu,
                                    downsample = 1,
                                    overdisperse = 1,
                                    meta_in = emb8_meta,
                                    grouping = emb8_meta$treatment,
                                    gene_table_in = emb8_gene_table,
                                    known_hits = getKnownAneu(emb8)){
  
  #use overarching fit for these
  whole_fit = get_fit(cbind(raw_data_control, raw_data_aneu))
  control_aneu = subsample(counts = raw_data_aneu, downsample_frac = downsample, overdisperse_factor = overdisperse, fit = whole_fit)
  control_norm = subsample(counts = raw_data_control, downsample_frac = downsample, overdisperse_factor = overdisperse, const_mu = T, fit = whole_fit)
  control = cbind(control_aneu, control_norm)
  
  
  #normalise by cpm
  # control = plot_test(sim_counts = control, real_counts = cbind(raw_data_control, raw_data_aneu), plot = F)
  control = sweep(control, 2, colSums(control), "/") * 1e6
  
  
  new_pl = makeAneu(counts = control, 
                    genes = rownames(control), 
                    chrs = gene_table_in$chromosome_name, 
                    cellNames = colnames(control), 
                    cellGroups = as.character(meta_in$treatment))
  
  #for high overdispersion, many many genes are rejected by the algorithm.
  #this breaks the algorithm - there are not genes on all the chr
  #we therefore set a high threshold to ensure we can assess high overdispersion
  new_pl = setParam(new_pl, "extreme.gene.thresh", 1000)
  
  new_pl = setKnownAneu(new_pl, known_hits)
  new_pl = doAneu(new_pl)
  
  
  out = testPerformance(new_pl)
  return(out)
}

#' subsample_repeat_analysis repeats subsampling of the mouse 
#' embryo data according to subsample() and tests the prediction accuracies of the 
#' simulated results
#'
#' @param raw_data_control raw counts (no aneu, using fixed mu)
#' @param raw_data_aneu raw counts (with aneu)
#' @param downsample downsample fraction
#' @param overdisperse overdispersion fraction
#' @param meta_in metadata for gt_test
#' @param gene_table_in gene table for gt_test
#' @param known_hits real aneuploidies for gt_test
#' @param nruns times to repeat analysis
#'
#' @return data frame with repeated metrics of success for each run
#' @export
#'
#' @examples
subsample_repeat_analysis = function(raw_data_control = raw_fine,
                                     raw_data_aneu = raw_aneu,
                                     downsample = 1,
                                     overdisperse = 1,
                                     meta_in = emb8_meta,
                                     grouping = emb8_meta$treatment,
                                     gene_table_in = emb8_gene_table,
                                     known_hits = getKnownAneu(emb8),
                                     nruns = 1) {
  require(parallel)
  
  subsamples = mclapply(1:nruns, function(x) subsample_repeat_wrapper(raw_data_control = raw_data_control,
                                                                      raw_data_aneu = raw_data_aneu,
                                                                      downsample = downsample,
                                                                      overdisperse = overdisperse,
                                                                      meta_in = meta_in,
                                                                      grouping = grouping,
                                                                      gene_table_in = gene_table_in,
                                                                      known_hits = known_hits), 
                        mc.cores = 5)#keep at 5 to limit memory consumption
  res = do.call(rbind, subsamples)
  return(res)
}




subsample_repeat_fixed_wrapper = function(raw_counts = getCounts(emb8),
                                          downsample = 1,
                                          overdisperse = 1,
                                          meta_in = emb8_meta,
                                          gene_table_in = emb8_gene_table,
                                          known_hits = getKnownAneu(emb8)){
  # print(paste("Run", i, "of", nruns))
  control = subsample_simulate_aneuploidy(
    counts = raw_counts,
    downsample_frac = downsample,
    overdisperse_factor = overdisperse,
    gene_map = gene_table_in,
    known_aneuploidy = known_hits
  )
  
  #normalise by cpm
  control = sweep(control, 2, colSums(control), "/") * 1e6
  
  new_pl = makeAneu(counts = control, 
                    genes = rownames(control), 
                    chrs = gene_table_in$chromosome_name, 
                    cellNames = colnames(control), 
                    cellGroups = as.character(meta_in$treatment))
  
  new_pl = setKnownAneu(new_pl, known_hits)
  new_pl = doAneu(new_pl)
  
  
  out = testPerformance(new_pl)

  return(out)
}

#' subsample_repeat_analysis_fixed repeats subsampling with fixed 
#' effects(1.5/0.5) of the mouse 
#' embryo data and tests the prediction accuracies of the 
#' simulated results
#'
#' @param counts raw counts
#' @param downsample downsample fraction
#' @param overdisperse overdispersion fraction
#' @param meta_in metadata for gt_test
#' @param gene_table_in gene table for gt_test
#' @param known_hits real aneuploidies for gt_test
#' @param nruns times to repeat analysis
#'
#' @return data frame with repeated metrics of success for each run
#' @export
#'
#' @examples
subsample_repeat_analysis_fixed = function(raw_counts = getCounts(emb8),
                                           downsample = 1,
                                           overdisperse = 1,
                                           meta_in = emb8_meta,
                                           gene_table_in = emb8_gene_table,
                                           known_hits = getKnownAneu(emb8),
                                           nruns = 1) {
  
  sims = mclapply(1:nruns, function(x) subsample_repeat_fixed_wrapper(raw_counts = raw_counts,
                                                                      downsample = downsample,
                                                                      overdisperse = overdisperse,
                                                                      meta_in = meta_in,
                                                                      gene_table_in = gene_table_in,
                                                                      known_hits = known_hits),
                  mc.cores = nruns)
  return(do.call(rbind, sims))
  
  }

#' get_simulated_counts wraps subsample in order to lapply it and get out matrices of counts
#'
#' @param dispersion_factor od_factor for the resampling
#' @param aneu_counts counts for aneuploid cells; leave default here
#' @param norm_counts counts for normal ploidy cells; leave default here
#'
#' @return matrix of normalised gene counts.
#' @export
#'
#' @examples
get_simulated_counts = function(dispersion_factor, aneu_counts = raw_aneu, norm_counts = raw_fine, normalise = TRUE, downsample_factor = 1){
  whole_fit = get_fit(cbind(aneu_counts, norm_counts))
  aneu_sim = subsample(counts = aneu_counts, downsample_frac = downsample_factor, overdisperse_factor = dispersion_factor, fit = whole_fit)
  control_sim = subsample(counts = norm_counts, downsample_frac = downsample_factor, overdisperse_factor = dispersion_factor, const_mu = T, fit = whole_fit)
  sim = cbind(control_sim, aneu_sim)
  if(normalise){
    sim_norm = sweep(sim*1e6, 2, colSums(sim), "/")
  } else {
    sim_norm = sim
  }
  sim_norm = sim_norm[order(rownames(sim_norm), decreasing = FALSE),]
  return(sim_norm)
}

# ASE FUNCTIONS ####


#' my_appendList recursively descends two lits x and val until it finds a non-matching level
#' (matched by name) whereupon it will make a list of each object from x and val at this level
#' #many thanks to 42- at http://stackoverflow.com/questions/9519543/merge-two-lists-in-r
#'
#' @param x should be CAST matrix list
#' @param val should be C57 matrix list
#'
#' @return new list with new lowest level[[cast/c57]]
#' @export
#'
#' @examples
my_appendList <- function (x, val) 
{
  stopifnot(is.list(x), is.list(val))
  
  xnames <- names(x)
  for (v in names(val)) {
    x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
      my_appendList(x[[v]], val[[v]])
    else list(count1=x[[v]], count2=val[[v]])
  }
  x
}



#' score_ase gets the fraction of count1 expression of the total count1 + count2 expression
#' per chromosome per cell
#'
#' @param count1 counts matrix
#' @param count2 counts matrix
#' @param chr_map maps genes to chromosomes as in other functions
#'
#' @return long data frame of cell, chr, ratio, and developmental stage(for deng only)
#' @export
#'
#' @examples
score_ase = function(count1 = cast_counts,
                     count2 = c57_counts,
                     chr_map = deng_genes) {
  
  # if(any(rownames(chr_map) != rownames(count1)) | any(rownames(count1) != rownames(count2)))
  #   stop("Please order gene map and counts the same")
  
  keep_genes = rownames(count1)[rownames(count1) %in% rownames(count2) & rownames(count1) %in% chr_map$gene]
  count1 = count1[match(keep_genes, rownames(count1))[!is.na(match(keep_genes, rownames(count1)))], ]
  count2 = count2[match(keep_genes, rownames(count2))[!is.na(match(keep_genes, rownames(count2)))], ]
  chr_map = chr_map[match(keep_genes, chr_map$gene)[!is.na(match(keep_genes, chr_map$gene))], ]
  
  #split by chromosome 
  count1_split = split(as.data.frame(count1), chr_map$chr)
  count2_split = split(as.data.frame(count2), chr_map$chr)
  
  chr_sum_count1 = lapply(count1_split, colSums)
  chr_sum_count2 = lapply(count2_split, colSums)
  
  frac_list = lapply(1:length(count1_split), function(x) chr_sum_count1[[x]]/(chr_sum_count1[[x]] + chr_sum_count2[[x]]))
  names(frac_list) = names(chr_sum_count1)

  fracs = unlist(frac_list)
  #get chr, cell etc
  chrsplit = strsplit(names(fracs), split = ".", fixed = T)
  chrs = sapply(chrsplit, function(x) as.numeric(x[1]))
  #for deng, we need to add .txt; not for mesc!
  if(grepl(".txt", names(fracs)[1]))
    names = sapply(chrsplit, function(x) paste0(x[2], ".txt"))
  else
    names = sapply(chrsplit, function(x) x[2])
  namesplit = strsplit(names, "_")
  stage = sapply(namesplit, function(x) x[2])
  
  out_df = data.frame(ratio = fracs, stage = stage, chr = chrs, cell = names)
  
  return(out_df)
  
}

# DIFFERENTIAL EXPRESSION WRAPPERS ####

#' get_de_genes is a wrapper for edgeR differential expression analysis
#'
#' @param raw_counts raw counts matrix, will be normalised by edgeR
#' @param mm model.matrix for edgeR
#'
#' @return list of up (significantly upregulated genes), down (significantly downregulated genes),
#' all (all genes considered), and full_tab (full results table)
#' @export
#'
#' @examples
get_de_genes = function(raw_counts, mm){
  require(edgeR)
  #gene appears in at least 10% of cells 
  #edgeR_counts = raw_counts[rowSums(raw_counts>10)>(ncol(raw_counts)/10),]
  #mean count (CPM norm) > 10
  cpm = sweep(raw_counts, 2, colSums(raw_counts), "/")*1e6
  edgeR_counts = raw_counts[rowMeans(cpm)>10,]
  
  edge = DGEList(counts = edgeR_counts, genes = rownames(edgeR_counts))
  edge = calcNormFactors(edge)
  edge = estimateDisp(edge, design = mm)
  fit = glmFit(edge, mm)
  test = glmLRT(fit, coef = 2)
  
  res = test$table
  res$fdr = p.adjust(res$PValue, method = "fdr")
  res = res[order(res$PValue),]
  # return(res)
  
  #old method above, now we just return lists of "up", "down" and "all"
  up = rownames(res[res$fdr<0.1 & res$logFC>0,])
  down = rownames(res[res$fdr<0.1 & res$logFC<0,])
  all = rownames(res)
  
  return(list(up = up, down = down, all = all, full_tab = res))
}



#' do_target_genes is a wrapper to test for DE and plot differences between a set of genes
#' from the paper Vera-Rodriguez et al. 2015
#'
#' @param raw_counts raw counts matrix - will be CPM normalised inside function
#' @param design_matrix model.matrix for edgeR
#' @param ortholog_df data.frame of orthologous genes that is defined in the Rmd script
#' @param naming_scheme "ensembl" or "mgi" depenging on the dataset
#' @param aneu_cells character vector of the aneuploid cells, for plotting
#'
#' @return DE results data frame for the gene subset only, and plots the boxplot comparisons
#' @export
#'
#' @examples
do_target_genes = function(raw_counts, design_matrix, ortholog_df = orthologs, naming_scheme = c("ensembl", "mgi"), aneu_cells, title = "missing title") {
  #note that these are alphabetical - should match well if we sort
  hgnc_up = c("BUB1", "CASP2", "GAPDH", "GADD45A")
  hgnc_down = c("BUB3", "CDK7", "CTNNB1", "E2F1", "PTTG1", "TP53", "TSC2", "YBX2")
  ortholog_df = ortholog_df[order(ortholog_df$hgnc),]
  ortholog_df$ensembl = ortholog_df$mmusculus_homolog_ensembl_gene
  if(naming_scheme[1]%in%c("mgi", "ensembl")){
    up_genes = ortholog_df[ortholog_df$aneu_up, naming_scheme]
    up_labels = ortholog_df[ortholog_df$aneu_up, "hgnc"]
    down_genes = ortholog_df[!ortholog_df$aneu_up, naming_scheme]
    down_labels = ortholog_df[!ortholog_df$aneu_up, "hgnc"]
  } else {
    stop("select a naming scheme, hgnc or ensembl")
  }
  gene_selection = c(up_genes, down_genes)
  gene_counts = raw_counts[gene_selection, ]
  
  norm_counts = sweep(gene_counts, 2, colSums(raw_counts), "/")*1e6
  
  melted = melt(norm_counts)
  names(melted) = c("gene", "cell", "cpm")
  melted$gene = as.character(melted$gene)
  melted$cell = as.character(melted$cell)
  melted$aneu = melted$cell %in% aneu_cells
  melted$logcpm = log2(melted$cpm + 1)
  melted$hgnc = sapply(melted$gene, function(x)
    ortholog_df$hgnc[ortholog_df[,naming_scheme] == x])
  
  plot = ggplot(data = melted,
                mapping = aes(x = factor(hgnc, levels = c(up_labels, down_labels)), y = logcpm, fill = aneu)) + 
    geom_boxplot() + geom_vline(xintercept = 4.5) +
    annotate("text", x = 3, y = max(melted$logcpm)*1.1, label = "Aneuploidy higher") + 
    annotate("text", x = 6, y = max(melted$logcpm)*1.1, label = "Aneuploidy lower") + 
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    ylab("Log2 CPM") +
    ggtitle(title) +
    scale_fill_manual(labels = c("Normal Ploidy", "Aneuploid"), name = "", values = c("TRUE"="tomato3", "FALSE"="dodgerblue3"))
  
  
  #do the DE
  edge = DGEList(counts = raw_counts, genes = rownames(raw_counts))
  edge = calcNormFactors(edge)
  edge = estimateDisp(edge, design = design_matrix)
  fit = glmFit(edge, design_matrix)
  test = glmLRT(fit, coef = 2)
  
  key_genes = test$table[c(up_genes, down_genes),]
  key_genes$fdr = p.adjust(key_genes$PValue, method = "fdr")
  
  test$table$fdr = p.adjust(test$table$PValue, method = "fdr")
  
  sig_index = which(key_genes$fdr<0.1)
  
  if(length(sig_index)>0){
    plot = plot +
      annotate("text", x = sig_index, y = -1, label = "*", size = 12, col = "black") +
      geom_hline(yintercept = 0, col = "darkgrey")
  }
  plot(plot)
  
  return(key_genes)
}
