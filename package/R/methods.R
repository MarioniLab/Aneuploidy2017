
#' Aneuploidy Constructor function
#'
#' This function is the constructor for the ploidytest object
#' @param counts Raw counts matrix (i.e. integers)
#' @param genes Genes for the rows of the counts matrix
#' @param chrs Chromosome number for each of the genes
#' @param cellNames Names for each of the columns of the counts matrix
#' @param cellGroups Cell groups for use in analysis (i.e. different cell types)
#'
#' @return Ploidytest object
#' @export
#'
#' @examples
makeAneu = function(counts = NULL,
                    genes = NULL,
                    chrs = NULL,
                    cellNames = NULL,
                    cellGroups = NULL,
                    params = NULL){
  
  if(any(is.null(counts), is.null(genes), is.null(cellNames), is.null(chrs), is.null(cellGroups))){
    stop("Please specify all the data needed.")
  }
  
  if(!is.matrix(counts)){
    stop("Please supply counts as a matrix")
  }
  if(!is.vector(genes)){
    stop("Please supply genes as a vector")
  }
  if(!is.vector(chrs)){
    stop("Please supply chrs as a vector")
  }
  if(!is.vector(cellNames)){
    stop("Please supply cellNames as a vector")
  }
  if(!is.vector(cellGroups)){
    stop("Please supply cellGroups as a vector")
  }
  
  if(nrow(counts) != length(chrs) | nrow(counts) != length(genes)){
    stop("There aren't as many genes or chr as count rows")
  }
  
  if(ncol(counts) != length(cellNames) | ncol(counts) != length(cellGroups)){
    stop("cellName or cellGroup doesn't match the number of count columns")
  }
  
  #take default if unspecified
  if(is.null(params))
    params = list(p.thresh = 0.1, 
                  min.deviation = 0.2, 
                  center.cells = TRUE, 
                  min.median = 50, 
                  extreme.gene.thresh = Inf)
  
  cpm = sweep(counts*1e6,
              2,
              apply(counts, 2, sum),
              "/")
  hiexp = cpm[apply(cpm, 1, median)>50,]
  
  #from the data analysis in the paper - reversine model
  base_coefs = c(int = 0.1897, grad = 0.7479)
  
  means = log10(rowMeans(hiexp))
  sds = log10(apply(hiexp, 1, sd))
  fitted = base_coefs["int"] + means * base_coefs["grad"]
  score = sum(sds - fitted)/length(sds)
  
  metrics = list("ngenes" = nrow(hiexp),
                 "zeros" = apply(hiexp, 1, function(x) sum( x == 0 )/length(x)),
                 "residual" = score)

  obj = new("ploidytest", 
            counts = counts, 
            genes = as.character(genes), 
            chrs = as.character(chrs), 
            cellNames = as.character(cellNames), 
            cellGroups = as.character(cellGroups), 
            cpm = cpm,
            params = params,
            metrics = metrics)
  
  
  return(obj)
  
}

#' Separate cell analysis groups
#'
#' This function splits your cells into the cellGroups specified in the constructor.
#' @param ploidytest Ploidytest object
#'
#' @return List of ploidytest objects, names according to the cellGroups category represented.
#' @export
#'
#' @examples
splitCellsByGroup = function(ploidytest){
  groups = ploidytest@cellGroups
  
  list = list()
  for(i in unique(groups)){
    new_obj = makeAneu(counts = ploidytest@counts[,groups == i, drop = FALSE], 
                       genes = ploidytest@genes, 
                       chrs = ploidytest@chrs, 
                       cellNames = ploidytest@cellNames[groups == i], 
                       cellGroups = groups[groups == i],
                       params = ploidytest@params)
    list[[length(list)+1]] = new_obj
    names(list)[length(list)] = i
  }
  return(list)
}

#' Subset genes in ploidytest object
#'
#' @param ploidytest Ploidytest object
#' @param logical_keep Logical vector of genes to retain (TRUE) and drop (FALSE)
#'
#' @return Ploidytest object with (probably) fewer genes
#' @export
#'
#' @examples
subsetGenes = function(ploidytest, logical_keep){
  return(makeAneu(counts = ploidytest@counts[logical_keep, ], 
                  genes = ploidytest@genes[logical_keep], 
                  chrs = ploidytest@chrs[logical_keep], 
                  cellNames = ploidytest@cellNames, 
                  cellGroups = ploidytest@cellGroups,
                  params = ploidytest@params))
}

#' Exclude lowly expressed and DE genes (if specified)
#' 
#' This function excludes genes from analysis and is called during doAneu()
#' @param ploidytest Ploidytest object
#'
#' @return Ploidytest object with fewer genes
#' @export
#'
#' @examples
trimGenes = function(ploidytest){
  
  cpm = getCPM(ploidytest)
  
  #Remove low CPM genes
  keep_cpm = apply(cpm, 1, median) >= getParam(ploidytest, "min.median")
  ploidytest = subsetGenes(ploidytest, logical_keep = keep_cpm)

  cpm = getCPM(ploidytest)
  
  #Get a_gij
  a = sweep(cpm, 
            1, 
            apply(cpm, 1, median), 
            "/")
  
  
  amax = apply(a, 1, max)
  extreme_keep = amax < getParam(ploidytest, "extreme.gene.thresh")
  
  return(subsetGenes(ploidytest, logical_keep = extreme_keep))
}

#' Calculate aneuploidy scores for one cell group
#'
#' This is the guts of the aneuploidy calling method, which is called by doAneu().
#' You probably don't want to call this yourself
#' @param ploidytest Ploidytest object
#'
#' @return ploidytest object with aneuploidy results
#' @export
#'
#' @examples
calcAneu = function(ploidytest){

  #Remove low CPM/extreme genes
  ploidytest = trimGenes(ploidytest)
  cpm = getCPM(ploidytest)
  
  #Median normalise
  a = sweep(cpm, 
            1, 
            apply(cpm, 1, median), 
            "/")
  
  #Sum on chr
  spt = split(as.data.frame(a), getCHR(ploidytest))
  b = t(sapply(spt, colSums))
  
  #Normalise by number of genes
  r = sweep(b, 1, sapply(spt, nrow), "/")
  
  #Cell centering
  if(getParam(ploidytest, "center.cells")){
    s = sweep(r, 2, apply(r, 2, median), "-") + 1
  } else {
    s = r
  }

  
  #Convert to Z
  medians = apply(s, 1, median)
  mads = apply(s, 1, mad)
  z_score = sweep(sweep(s, 1, medians, "-"), 1, mads, "/")
  
  #Prepare output
  out = reshape2::melt(z_score)
  names(out) = c("chr", "cell", "z")
  out$score = reshape2::melt(s)$value
  #we do a one tailed test on the absolute value, then double it, as we want a two-tailed value.
  out$p =  sapply(abs(out$z), pt, df = ncol(z_score)-1, lower.tail = F) * 2

  #return calls
  return(out)
}

#' Perform Aneuploidy Assessment
#'
#' doAneu performs the aneuploidy assessment.
#' @param ploidytest Ploidytest object
#'
#' @return Ploidytest object with results
#' @export
#'
#' @examples
doAneu = function(ploidytest){
  #split by cell group
  spt = splitCellsByGroup(ploidytest)
  
  #run processing
  results = do.call(rbind, lapply(spt, calcAneu))
  
  #tidy up output
  results$p.adj = p.adjust(results$p, method = "fdr")
  results$monosomy = results$z < 0

  hits = results[results$p.adj < getParam(ploidytest, "p.thresh") & abs(results$score-1)>getParam(ploidytest, "min.deviation"), ]
  
  #merge results into object
  ploidytest@scores = results
  ploidytest@aneuploidies = hits
  
  #return result-populated object
  return(ploidytest)
}

#' Test the method's performance where aneuploidy is known
#'
#' @param ploidytest Ploidytest object
#'
#' @return Vector of performance metrics
#' @export
#'
#' @examples
testPerformance = function(ploidytest){
  
  if(nrow(getKnownAneu(ploidytest)) == 0)
    stop ("There are no known aneuploidies")
  
  if(nrow(getScores(ploidytest)) == 0)
    stop("You must run doAneu() first")

  hits = getHits(ploidytest)
  truth = getKnownAneu(ploidytest)
  
  hit_chr = paste(hits$cell, hits$chr, hits$monosomy)
  true_chr = paste(truth$cell, truth$chr, truth$monosomy)
  
  TP = sum(hit_chr %in% true_chr)
  FP = sum(! hit_chr %in% true_chr)
  FN = sum(! true_chr %in% hit_chr)
  TN = ncol(getCPM(ploidytest)) * length(unique(getCHR(ploidytest))) - (TP + FP + FN)

  return( c("sensitivity" =  TP / (TP + FN),
            "precision" = TP / (TP + FP),
            "fdr" = 1 - (TP / (TP + FP)),
            "specificity" = TN / (TN + FP),
            "accuracy" = (TP + TN) / (TP + TN + FP + FN),
            "f1" = 2 * TP / (2 * TP + FP + FN),
            "fpr" =  FP/(FP+TN) ) )
}


#' Plot a PCA of log-CPM counts
#'
#' @param ploidytest Ploidytest object
#' @param cols Colour labels e.g. cell line or embryo number
#'
#' @return invisible(0)
#' @export
#'
#' @examples
plotPCA = function(ploidytest, cols = NULL){
  
  ploidytest = trimGenes(ploidytest)
  pca = prcomp(t(log10(getCPM(ploidytest) + 1)))
  vars = pca$sdev^2
  
  if(is.null(cols))
    cols = getGroups(ploidytest)
  
  p = ggplot2::ggplot(as.data.frame(pca$x[,1:2]), ggplot2::aes(x=PC1, y=PC2, col = factor(cols))) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() + 
        ggplot2::labs(x = paste0("PC1, ", format(vars[1]/sum(vars)*100, digits = 3), "% variance"),
                 y = paste0("PC2, ", format(vars[2]/sum(vars)*100, digits = 3), "% variance")) +
        ggplot2::scale_color_brewer(palette = "Set1", name = "")
    
  print(p)
  invisible(0)
}

#' Present False Negative Results
#'
#' Function returns a nicely formatted data frame of false negative calls, when aneuploidies are known
#' @param ploidytest Ploidytest object
#'
#' @return data frame of false negative calls
#' @export
#'
#' @examples
presentFN = function(ploidytest){
  if(nrow(getScores(ploidytest)) == 0)
    stop("Run doAneu first")
  known = getKnownAneu(ploidytest)
  fn = getFN(ploidytest)
  out = fn[, c("cell", "chr", "score", "p.adj", "monosomy")]
  out$sig_p = out$p.adj<getParam(ploidytest, "p.thresh")
  out$sig_score = abs(out$score-1) > getParam(ploidytest, "min.deviation")
  po = paste(out$cell, out$chr)
  pk = paste(known$cell, known$chr)
  out$dir_correct = sapply(1:nrow(out), function(x) out$monosomy[x] == known$monosomy[match(po[x], pk)] )
  return(out)
}

#' Present Known Aneuploidy Scores
#'
#' Function returns a nicely formatted data frame of scores for true aneuploidies
#' @param ploidytest Ploidytest object
#'
#' @return data frame of true aneuploidy calls
#' @export
#'
#' @examples
presentKnown = function(ploidytest){
  if(nrow(getScores(ploidytest)) == 0)
    stop("Run doAneu first")
  known = getKnownAneu(ploidytest)
  scores = getScores(ploidytest)
  out = scores[paste(scores$cell, scores$chr) %in%
                 paste(known$cell, known$chr), 
               c("cell", "chr", "score", "p.adj", "monosomy")]
  names(out)[which(names(out) == "monosomy")] = "pred_monosomy"
  out$sig_p = out$p.adj<getParam(ploidytest, "p.thresh")
  out$sig_score = abs(out$score-1) > getParam(ploidytest, "min.deviation")
  po = paste(out$cell, out$chr)
  pk = paste(known$cell, known$chr)
  out$dir_correct = sapply(1:nrow(out), function(x) out$pred_monosomy[x] == known$monosomy[match(po[x], pk)] )
  return(out)
}

#' Randomly downsample cells
#'
#' @param ploidytest Ploidytest object
#' @param n Number of cells to sample
#'
#' @return Ploidytest object with fewer cells
#' @export
#'
#' @examples
sampleCells = function(ploidytest, n){
  if(n > ncol(scploid::getCounts(ploidytest)))
    stop("n larger than the number of cells")
  
  smp = sample(ncol(scploid::getCounts(ploidytest)), n, replace = FALSE)
  
  return(makeAneu(counts = scploid::getCounts(ploidytest)[, smp], 
                  genes = getGenes(ploidytest), 
                  chrs = getCHR(ploidytest), 
                  cellNames = getCellNames(ploidytest)[smp], 
                  cellGroups = getGroups(ploidytest)[smp],
                  params = ploidytest@params))
}

#' Provides broad-stroke data quality assessment
#'
#' assessMetrics looks at the residual score, the number of genes considered, and the fraction of 0's
#' in the considered genes and provides a simple quality assessment on each. This assumes default parameter choice.
#' @param ploidytest ploidytest object
#'
#' @return Vector of quality calls
#' @export
#'
#' @examples
assessMetrics = function(ploidytest){
  metrics = getMetrics(ploidytest)
  
  resid = ifelse(metrics$residual < 0.1, "Good quality",
                  ifelse(metrics$residual<0.3, "Exercise Caution", "Poor quality"))
  
  ngenes = ifelse(metrics$ngenes > 3000, "Good quality", "Poor quality")
  
  zeros = ifelse(quantile(metrics$zeros, 0.75)<0.05, "Good quality", "Poor quality")

  out = c(resid, ngenes, zeros)
  names(out) = c("residuals", "ngenes", "zeros")

  return(out)
  
}

