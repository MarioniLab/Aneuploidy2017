

#' Get counts matrix from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Counts matrix
#' @export
#'
#' @examples
getCounts = function(ploidytest){
  out = ploidytest@counts
  colnames(out) = ploidytest@cellNames
  rownames(out) = ploidytest@genes
  return(out)
}

#' Get CPM matrix from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return CPM matrix
#' @export
#'
#' @examples
getCPM = function(ploidytest){
  out = ploidytest@cpm
  colnames(out) = ploidytest@cellNames
  rownames(out) = ploidytest@genes
  return(out)
}

#' Get cell names from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Cell name vector
#' @export
#'
#' @examples
getCellNames = function(ploidytest){
  return(ploidytest@cellNames)
}


#' Get cell groups from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Cell group vector
#' @export
#'
#' @examples
getGroups = function(ploidytest){
  return(ploidytest@cellGroups)
}


#' Get chromosome locations for each gene from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Chromosome vector for each gene
#' @export
#'
#' @examples
getCHR = function(ploidytest){
  return(ploidytest@chrs)
}

#' Get gene names from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Gene name vector
#' @export
#'
#' @examples
getGenes = function(ploidytest){
  return(ploidytest@genes)
}

#' Get parameters for aneuploidy calling
#' 
#' @param ploidytest Ploidytest object
#' 
#' @return Parameter list. Contains:
#' \itemize{
#' \item "p.thresh": post FDR correction p-value for aneuploidy calling (default 0.1)
#' \item "min.deviation": minimum deviation of s_{ij} from 1 for aneuploidy calling (default 0.2)
#' \item "center.cells": logical value to allow or prevent within-cell s_{ij} centering (default TRUE)
#' \item "min.median": minimum CPM median value for gene retention for analysis (default 50)
#' \item "extreme.gene.thresh": maximum a_{gij} that one gene may show before removal as outlier (default Inf)
#' }
#' @export
#' 
#' @examples
getParams = function(ploidytest){
  return(ploidytest@params)
}

#' Get known aneuploidies from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Data frame of known aneuploidies
#' @export
#'
#' @examples
getKnownAneu = function(ploidytest){
  return(ploidytest@knownAneu)
}

#' Get prediction output from ploidytest object
#'
#' @param ploidytest Ploidytest object
#'
#' @return Data frame of prediction outputs for all chromosomes. Contains:
#' \itemize{
#' \item chr: chromosome
#' \item cell: cell
#' \item z: z-score for p-value calculation
#' \item score: s_{ij}
#' \item p: raw p-value
#' \item p.adj: p-value FDR-adjusted 
#' \item monosomy: TRUE if z<0 i.e. reduced expression - *not* an indicator of aneuploidy
#' }
#' @export
#'
#' @examples
getScores = function(ploidytest){
  if(nrow(ploidytest@scores) == 0)
    stop("Run doAneu first")
  
  return(ploidytest@scores)
}

#' Get prediction output from ploidytest object for predicted aneuploidies only
#'
#' @param ploidytest Ploidytest object
#'
#' @return Data frame of prediction outputs for aneuploid-called chromosomes. Contains:
#' chr: chromosome
#' cell: cell
#' z: z-score for p-value calculation
#' score: s_{ij}
#' p: raw p-value
#' p.adj: p-value FDR-adjusted 
#' monosomy: TRUE if z<0 i.e. reduced expression
#' @export
#'
#' @examples
getHits = function(ploidytest){
  if(nrow(getScores(ploidytest)) == 0)
    stop("Run doAneu first")
  return(ploidytest@aneuploidies)
}

#' Set parameters for aneuploidy calling
#' 
#' @param ploidytest Ploidytest object
#' @param param_name Name of parameter. May be any of:
#' \itemize{
#' \item "p.thresh": post FDR correction p-value for aneuploidy calling (default 0.1)
#' \item "min.deviation": minimum deviation of s_{ij} from 1 for aneuploidy calling (default 0.2)
#' \item "center.cells": logical value to allow or prevent within-cell s_{ij} centering (default TRUE)
#' \item "min.median": minimum CPM median value for gene retention for analysis (default 50)
#' \item "extreme.gene.thresh": maximum a_{gij} that one gene may show before removal as outlier (default Inf)
#' }
#' @param param_value New parameter value
#' @param print Prints the changes to the numbers of genes considered if you change extreme.gene.thresh if TRUE
#' 
#' @return Ploidytest object with revised parameter value
#' @export
#' 
#' @examples
setParam = function(ploidytest, param_name, param_value, print = FALSE){
  if(!param_name %in% names(getParams(ploidytest)))
    stop("parameter name doesn't match anything")
  
  if( param_name != "center.cells" & !is.numeric(param_value) )
    stop("Please insert a numeric parameter value")
  
  if( param_name == "center.cells" & !is.logical(param_value) )
    stop("Please use a logical parameter value")
  
  if(length(param_value)!=1)
    stop("Please use param_value with length 1")
  
  ploidytest@params[[as.character(param_name)]] = param_value
  
  if(param_name == "extreme.gene.thresh" & print){
    #first split, then perform calcs separately for each as in the real analysis
    spt = splitCellsByGroup(ploidytest)
    expr_genes = lapply(spt, function(x) 
      apply(getCPM(x), 1, median) > getParam(x, "min.median"))
    
    old_chr = lapply(1:length(spt), function(x) table(getCHR(spt[[x]])[expr_genes[[x]]]))
    old_chr = lapply(old_chr, function(x) as.vector(x[order(as.numeric(names(x)))]))
    
    new_chr = lapply(spt, function(x) table(getCHR(trimGenes(x))))
    new_chr = lapply(new_chr, function(x) as.vector(x[order(as.numeric(names(x)))]))
    
    tabs = lapply(1:length(spt), function(x) matrix(c(old_chr[[x]], old_chr[[x]] - new_chr[[x]]),
                                                    nrow = 2, byrow = TRUE,
                                                    dimnames = list(
                                                      paste(names(spt)[x],
                                                            c("genes considered",
                                                              "extreme genes dropped")),
                                                      1:length(old_chr[[x]]))
                                                    ))
    
    tab = do.call(rbind, tabs)
    print(tab)
  }
  

  return(ploidytest)
}

#' Get parameter for aneuploidy calling
#' 
#' @param ploidytest Ploidytest object
#' @param param_name Name of parameter (see ?getParams or ?setParams)
#' 
#' @return named parameter value
#' @export
#' 
#' @examples
getParam = function(ploidytest, param_name){
  if(!param_name %in% names(getParams(ploidytest)))
    stop("parameter name doesn't match anything")
  
  return(ploidytest@params[[as.character(param_name)]])
}

#' Add known aneuploidy data to ploidytest object
#'
#' @param ploidytest Ploidytest object
#' @param aneu_df Data frame of aneuploidies. Requires columns:
#' \itemize{
#' \item "chr": Chromosome
#' \item "cell": Cell name
#' \item "monosomy": TRUE if monosomy, FALSE if trisomy
#' }
#'
#' @return Ploidytest object with known aneuploidies
#' @export
#'
#' @examples
setKnownAneu = function(ploidytest, aneu_df){
  if(!is.data.frame(aneu_df))
    stop("Please supply a dataframe")
  
  if(any(! c("cell", "chr", "monosomy") %in% names(aneu_df)))
    stop("data frame needs cell, chr, and monosomy columns")
  
  if(any(! aneu_df$cell %in% getCellNames(ploidytest)))
    stop("Cells are in the aneu_df data frame that are not in cellNames")
  
  if(any(! aneu_df$chr %in% getCHR(ploidytest)))
    stop("Chromosoms are in the aneu_df data frame that are not in chr")
  
  ploidytest@knownAneu = aneu_df
  return(ploidytest)
}

#' Get mean-variance relationship of genes
#'
#' Function gets the mean and variance of expression for each gene's CPM that would be used in aneuploidy
#' assessment
#' @param ploidytest Ploidytest object
#'
#' @return Data frame of mean (gene mean CPM) and sd (gene CPM sd)
#' @export
#'
#' @examples
getMeanVar = function(ploidytest){
  counts = getCPM(trimGenes(ploidytest))
  return(data.frame(mean = rowMeans(counts), sd = apply(counts, 1, sd)))
}

#' Get call scores for False positive aneuploidy calls
#'
#' @param ploidytest Ploidytest object with known aneuploidies
#'
#' @return Score data frame for false positive calls only
#' @export
#'
#' @examples
getFP = function(ploidytest){
  if(nrow(getScores(ploidytest)) == 0)
    print("You must run doAneu first")
  
  hits = getHits(ploidytest)
  known = getKnownAneu(ploidytest)
  
  return(hits[!paste(hits$cell, hits$chr, hits$monosomy) %in% 
         paste(known$cell, known$chr, known$monosomy), ])
}

#' Get call scores for False negative aneuploidy calls
#'
#' @param ploidytest Ploidytest object with known aneuploidies
#'
#' @return Score data frame for false negative calls only
#' @export
#'
#' @examples
getFN = function(ploidytest){
  if(nrow(getScores(ploidytest)) == 0)
    print("You must run doAneu first")
  
  hits = getHits(ploidytest)
  known = getKnownAneu(ploidytest)
  score = getScores(ploidytest)
  
  pk = paste(known$cell, known$chr, known$monosomy)
  ph = paste(hits$cell, hits$chr, hits$monosomy)
  #which are known, but not in hits
  fn = pk[! pk %in% ph]
  
  return(score[paste(score$cell, score$chr, score$monosomy) %in% fn, ])
}


#' Get call scores for all true aneuploidies
#'
#' @param ploidytest Ploidytest object with known aneuploidies
#'
#' @return Score data frame for true aneuploidies
#' @export
#'
#' @examples
getTrueScores = function(ploidytest){
  if(nrow(getScores(ploidytest)) == 0)
    print("You must run doAneu first")
  
  scores = getScores(ploidytest)
  known = getKnownAneu(ploidytest)
  
  return(scores[!paste(scores$cell, scores$chr, scores$monosomy) %in% 
                paste(known$cell, known$chr, known$monosomy), ])
}

#' Get matrix of sij
#'
#' @param ploidytest Ploidytest object
#'
#' @return Sij matrix, rows chromosomes, columns cells
#' @export
#'
#' @examples
getSijMat = function(ploidytest){
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
  
  return(s)
}

#' Get data metrics to estimate performance
#'
#' @param ploidytest Ploidytest object
#'
#' @return list of data metrics:
#' \itemize{
#' \item "ngenes": number of genes that qualify for aneuploidy assessment (median CPM>50 only)
#' \item "zeros": fraction of 0 counts for each gene that qualifies for aneuploidy assessment (median CPM>50 only) 
#' \item "residual": residual score vs. G&T-seq 8-cell stage embryos (see supplemental document of paper)
#' }
#' @export
#'
#' @examples
getMetrics = function(ploidytest){
  return(ploidytest@metrics)
}

#' Get a gene-CHR map
#'
#' @param ploidytest Ploidytest object
#'
#' @return data frame with gene and chr columns
#' @export
#'
#' @examples
getGeneTable = function(ploidytest){
  return(data.frame(gene = getGenes(ploidytest), chr = getCHR(ploidytest), row.names = getGenes(ploidytest)))
}

#' Gets the maximum value of a_{gij} that each considered gene contributes
#'
#' @param ploidytest ploidytest object
#'
#' @return Vector of all genes considered for analysis in any batch, named with the gene ID, containing the maximum
#' a_{gij} score observed
#' @export
#'
#' @examples
getMaxA = function(ploidytest){
  
  spt = splitCellsByGroup(ploidytest)
  #exclude specified genes
  spt = lapply(spt, trimGenes)
  as = lapply(spt, function(x) sweep(getCPM(x), 
                                     1, 
                                     apply(getCPM(x), 1, median), 
                                     "/"))
  maxas = lapply(as, function(x) apply(x, 1, max))
  
  dfs =lapply(maxas, function(x) data.frame(gene = names(x), maxa = x))
  df = do.call(rbind, dfs)
  df = df[order(df$maxa, decreasing = TRUE),]
  df = df[!duplicated(df$gene), ]
  
  out = df$maxa
  names(out) = df$gene
  
  if(any(out>20))
    warning("At least one gene scores a maximum a_{gij} value of 20. Consider removing it.")

  return(out)
}

