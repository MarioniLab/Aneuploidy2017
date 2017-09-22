# This script turns counts data into R objects that are saved in .RData files
# These files are read into the Rmd file for analysis

# LOAD AND PROCESS G&T DATA ####

raw_counts = read.table(
  paste0(folder_location, "raw_data/gt_8cell_counts.txt"),
  header = T,
  row.names = 1
)

raw_control_counts = read.table(
  paste0(folder_location, "raw_data/gt_8cell_control_counts.txt"),
  header = T,
  row.names = 1
)

raw_meta = read.table(
  paste0(folder_location, "raw_data/gt_8cell_meta.csv"),
  header = T,
  sep = ","
)

known_hits = read.csv(paste0(folder_location, "raw_data/known_gt_8cell_aneu.csv"), row.names = 1)

# METADATA
#remove trailing whitespace
raw_meta = raw_meta[-which(raw_meta$Sample.name == ""),]

#select interesting columns
meta_columns = c("Sample.name", "Embryo", "QC.Pass.fail", "Treatment")
meta = raw_meta[, meta_columns]

names(meta) = c("cell", "embryo", "qc", "treatment")

#tidy up some of the entries to be more concise
meta$embryo = sub("Embryo ", "", meta$embryo)
meta$cell = sub("8_cell_", "", meta$cell)

meta = meta[meta$qc == "Pass",]

# CHR MAP

gene_table = getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  mart = mouse_ensembl,
  values = as.character(rownames(raw_counts)),
  filters = "ensembl_gene_id"
)
rownames(gene_table) = gene_table$ensembl_gene_id
#remove sex chromosomes
gene_table = gene_table[gene_table$chromosome_name %in% c(1:19),]

# COUNTS

all_raw_counts = cbind(raw_counts, raw_control_counts)

#match up the counts to the metadata and gene table
split_cell = strsplit(colnames(all_raw_counts), ".", fixed = TRUE)
colnames(all_raw_counts) = sapply(split_cell, function(x) paste0(substr(x[[1]], nchar(x[[1]]), nchar(x[[1]])),
                                                             substr(x[[2]], nchar(x[[2]]), nchar(x[[2]]))))

all_raw_counts = all_raw_counts[gene_table$ensembl_gene_id, meta$cell]


emb8 = makeAneu(counts = as.matrix(all_raw_counts), 
                genes = gene_table$ensembl_gene_id, 
                chrs = gene_table$chromosome_name, 
                cellNames = meta$cell, 
                cellGroups = as.character(meta$treatment))

emb8 = setKnownAneu(emb8, known_hits)
emb8_meta = meta
emb8_gene_table = gene_table

save(
  emb8,
  emb8_meta,
  emb8_gene_table,
  file = paste0(folder_location, "proc_data/emb8_data.RData")
)

# LOAD AND PROCESS HCC38 CELL LINE DATA ####

hcc = read.table(paste0(folder_location, "raw_data/hcc38.mtx"), header = T, row.names = 1)

#keep only the cells that were present in the G&T figure in the paper
cells_38 = paste0("Cell_", c(38, 28, 77, 73, 86, 61, 62, 74, 29, 42, 66, 50, 14, 39, 65, 40, 49, 75, 5, 16, 41, 30, 51))
cells_bl = paste0("Cell_", c(43, 44, 69, 33, 80, 57, 31, 93, 58, 46, 72, 55, 59, 48, 21, 32, 95, 83, 35, 60, 94, 91, 23, 
                            24, 82, 11, 36, 84, 71, 45, 67, 68, 22, 70, 12, 56, 79, 19, 81))

hcc = hcc[,colnames(hcc)%in%c(cells_38, cells_bl)]

#keep only the human genes
hcc = hcc[grepl("ENSG", rownames(hcc)),]
#remove the genes with 0 reads
hcc = hcc[apply(hcc, 1, sum)>0, ]


#make meta
meta = data.frame(row.names = c(cells_38, cells_bl), 
                  cell = c(cells_38, cells_bl), 
                  line = c(rep("HCC38", length(cells_38)), rep("HCC38-BL", length(cells_bl)))
)

#make the map for gene location
gene_table = getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  mart = ensembl,
  values = as.character(rownames(hcc)),
  filters = 'ensembl_gene_id'
)

#keep only autosomal genes + order to meta
gene_table = gene_table[gene_table$chromosome_name%in%1:22,]
rownames(gene_table) = gene_table$ensembl_gene_id
hcc = hcc[rownames(gene_table), meta$cell]


#note the aneuploidies from the paper
hcc_known_aneuploidies = data.frame(
  cell = paste0("Cell_", c(68, 22, 70, 12, 82, 56, 79)),
  chr = c(rep(11, 4), rep(16, 3)),
  monosomy = c(rep(F, 4), F, rep(T, 2))
)





hcc = makeAneu(counts = as.matrix(hcc), 
                genes = gene_table$ensembl_gene_id, 
                chrs = gene_table$chromosome_name, 
                cellNames = as.character(meta$cell), 
                cellGroups = as.character(meta$line))

lst = splitCellsByGroup(hcc)
hcc38 = lst[["HCC38"]]
hccbl = lst[["HCC38-BL"]]

hcc = setKnownAneu(hcc, hcc_known_aneuploidies)
hccbl = setKnownAneu(hccbl, hcc_known_aneuploidies[hcc_known_aneuploidies$cell %in% cells_bl, ])
hcc_meta = meta

#get the copy number info for HCC38 (relative to the numbered cell, here)
rel_5 = read.table(paste0(folder_location, "raw_data/hcc38_rel_cell5.txt"), header = T)

get_abs_cn = function(cn_mat, ref_col = "Cell_5"){
  #change the names to match the meta
  names(cn_mat) = sub("HCC38", "Cell", names(cn_mat))
  #get the columns to do weighted mean
  cell_columns = which(grepl("Cell", names(cn_mat)))
  #only autosomes
  cn_mat = cn_mat[cn_mat$chr!="X",]
  #add length
  cn_mat$length = cn_mat$end - cn_mat$start
  #make absolute cn calls
  #number is X = Cell_ref - Cell_N
  #hence absolute number is Cell_N = Cell_ref - X
  ref_col = which(colnames(cn_mat)==ref_col)
  other_cells = cell_columns[cell_columns!=ref_col]
  cn_mat[, other_cells] = cn_mat[, ref_col] - cn_mat[, other_cells]
  
  #split to chr
  split_cn_mat = split(cn_mat, cn_mat$chr)
  
  #get relative length for each chr
  split_cn_mat = lapply(split_cn_mat, function(x) {
    x$rel_length = x$length/sum(x$length)
    return(x)}
  )
  
  #do the weighted mean
  means = lapply(split_cn_mat, function(x){
    apply(x[,cell_columns], 2, weighted.mean, w = x$rel_length)
  })
  
  unlisted = unlist(means)
  chrs = sapply(strsplit(names(unlisted), split = ".", fixed = T), function(x) x[1])
  cells = sapply(strsplit(names(unlisted), split = ".", fixed = T), function(x) x[2])
  
  #make and fill matrix
  out_mat = matrix(0, nrow = 22, ncol = length(cell_columns), dimnames = list(c(1:22), colnames(cn_mat)[cell_columns]))
  
  for(i in 1:length(unlisted)){
    out_mat[which(rownames(out_mat) == chrs[i]), which(colnames(out_mat) == cells[i])] = unlisted[i]
  }
  
  
  return(out_mat)
  
}

cn_base_5 = get_abs_cn(rel_5, ref_col = "Cell_5")

chr_medians = apply(cn_base_5, 1, median)
deviance = sweep(cn_base_5, MARGIN = 1, FUN = "-", STATS = chr_medians)

hits = which(abs(deviance)>=1)

cn_change = data.frame(cell = colnames(deviance)[col(deviance)[hits]],
                       chr = rownames(deviance)[row(deviance)[hits]],
                       stringsAsFactors = F)

cn_change$monosomy = deviance[hits]<0

hcc38 = setKnownAneu(hcc38, cn_change)

save(hcc,
     hcc_meta,
     hccbl,
     hcc38,
     cn_base_5,
     file = paste0(folder_location, "proc_data/hcc_data.RData"))

# LOAD AND PROCESS TRIS21 CELL LINE DATA ####

counts = read.table(paste0(folder_location, "raw_data/tris13_counts.mtx"), row.names = 1, header = T)

tris = counts[, ( grepl("Control", colnames(counts)) | 
                        grepl("T21", colnames(counts)) ) ]

#make meta
meta = data.frame(cell = as.character(colnames(tris)),
                       tris = ifelse(grepl("T21", colnames(tris)), "T21", "Diploid"),
                       row.names = colnames(tris))



#keep only the human genes
tris = tris[grepl("ENSG", rownames(tris)),]
#remove the genes with 0 reads
tris = tris[apply(tris, 1, sum)>0, ]

#make the map for gene location
gene_table = getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  mart = ensembl,
  values = as.character(rownames(tris)),
  filters = 'ensembl_gene_id'
)
#keep only autosomal genes
gene_table = gene_table[gene_table$chromosome_name%in%1:22,]
gene_table$chr = gene_table$chromosome_name
rownames(gene_table) = gene_table$ensembl_gene_id
tris = tris[rownames(gene_table), meta$cell]

t21 = makeAneu(counts = as.matrix(tris), 
               genes = gene_table$ensembl_gene_id, 
               chrs = gene_table$chromosome_name, 
               cellNames = as.character(meta$cell), 
               cellGroups = as.character(meta$tris))

# COPY NUMBER
#these csv files were available from the figure in the online version of the paper
norm_cn = read.table(paste0(folder_location, "raw_data/tris_copynumbers_natmeth/tris_norm_cn.csv"), sep = ",", header = T)
cn = read.table(paste0(folder_location, "raw_data/tris_copynumbers_natmeth/tris_tris_cn.csv"), sep = ",", header = T)

get_cn_table = function(cn_df){
  
  
  cn_df = cn_df[,-2]
  cn_df = cn_df[cn_df$Chr!="X",]
  cn_df[,1] = as.numeric(as.character(cn_df[,1]))
  
  split_chr = split(cn_df, cn_df$Chr)
  
  # #remove trisomy 21 character if it is present - we will do delta copy number later
  # if(any(grepl("Trisomy", names(split_chr[[21]]))))
  #   split_chr[[21]] = split_chr[[21]] - 1
  
  split_chr = lapply(split_chr, function(x) x[,-1])
  
  means = lapply(split_chr, colMeans)
  
  unlist = unlist(means) - 2
  #take those with a copy number that rounds to CN of 1 or >3
  candidates = unlist[abs(unlist)>0.5]
  
  split_names = strsplit(names(candidates), ".", fixed = T)
  chrs = sapply(split_names, function(x) x[1])
  cells = sapply(split_names, function(x) x[2])
  
  return( data.frame(cell = cells, chr = chrs, monosomy = candidates<0) )
  
  
}

known_aneuploidies = rbind(get_cn_table(norm_cn), get_cn_table(cn))
known_aneuploidies$cell = as.character(known_aneuploidies$cell)
known_aneuploidies$chr = as.numeric(as.character(known_aneuploidies$chr))

#make names match properly
known_aneuploidies$cell = gsub("_", "_Cell_", known_aneuploidies$cell)
known_aneuploidies$cell = gsub("Trisomy", "T", known_aneuploidies$cell)
#exclude a cell that isn't in the counts
known_aneuploidies = known_aneuploidies[known_aneuploidies$cell %in% getCellNames(t21),]

t21 = setKnownAneu(t21, known_aneuploidies)
tris_meta = meta
tris_meta$tris = tris_meta$cell %in% getKnownAneu(t21)$cell
save(t21,
     tris_meta,
     file = paste0(folder_location, "proc_data/tris_data.RData"))


# LOAD AND PROCESS MESC CELL DATA ####

mesc_counts = read.table(
  paste0(folder_location, "raw_data/mesc.csv"),
  header = T,
  row.names = 1
)
mesc_counts = as.matrix(mesc_counts)

colnames(mesc_counts) = sub(pattern = ".counts", replacement = "", x = colnames(mesc_counts))

#METADATA
#We need to extract the information from cell names:
breakup = strsplit(x = colnames(mesc_counts), split = "_")
treatment = sapply(breakup, function(x) x[3])
batch = sapply(breakup, function(x) x[4])
cell_number = sapply(breakup, function(x) x[5])

#get information into a metadata data frame
mesc_meta = data.frame(
  cell = colnames(mesc_counts),
  treatment = treatment,
  batch = batch,
  embryo = paste(treatment, batch),
  row.names = colnames(mesc_counts),
  cellnum = cell_number
)

# COUNTS

#remove rows of 0's
mesc_counts = mesc_counts[rowSums(mesc_counts) > 0,]
#make gene table
mesc_gene_table = getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  mart = mouse_ensembl,
  values = as.character(rownames(mesc_counts)),
  filters = "ensembl_gene_id"
)
#keep only the genes on autosomal chromosomes
mesc_gene_table = mesc_gene_table[mesc_gene_table$chromosome_name %in% c(1:19),]
rownames(mesc_gene_table) = mesc_gene_table$ensembl_gene_id
# rownames(mesc_gene_table) = mesc_gene_table$ensembl_gene_id
# #add a chr column for compatability with functions
# mesc_gene_table$chr = mesc_gene_table$chromosome_name
# #keep only genes that we have found on autosomes
# mesc_counts = mesc_counts[rownames(mesc_counts)%in%rownames(mesc_gene_table), ]


mesc = makeAneu(counts = as.matrix(mesc_counts[mesc_gene_table$ensembl_gene_id, as.character(mesc_meta$cell)]), 
                genes = mesc_gene_table$ensembl_gene_id, 
                chrs = mesc_gene_table$chromosome_name, 
                cellNames = as.character(mesc_meta$cell), 
                cellGroups = as.character(paste(mesc_meta$treatment, mesc_meta$batch, sep = ", B")))

# ALLELE

#ASE data is in a snp-by-snp format
#These are a few functions to collapse down the SNP level reads into gene level ones

sum_split = function(split_element){
  out = vector(length = ncol(split_element))
  
  for(col in 1:ncol(split_element)){
    if(is.numeric(split_element[,col])){
      out[col] = sum(as.numeric(split_element[,col]))
    } else {
      out[col] = as.character(split_element[1, col])
    }
  }
  return(out)
}

collapse_duplicates = function(df){
  df$gene_stable_id = as.character(df$gene_stable_id)
  rows = length(unique(df$gene_stable_id))
  
  out = matrix(NA, nrow = rows, ncol = ncol(df))
  
  split = split(df, df$gene_stable_id)
  
  summed_list = lapply(split, sum_split)
  
  summed = data.frame(matrix(unlist(summed_list), nrow = length(summed_list), byrow=T),stringsAsFactors=FALSE)
  names(summed) = names(df)
  
  
  for(col in 11:ncol(summed)){
    summed[,col] = as.numeric(summed[,col])
  }
  
  return(summed)
  
}

process_mesc_allele = function(table, filename){
  #remove duplicate SNPS
  table = table[!duplicated(paste(table$chr, table$pos)),]
  table = table[table$qual=="PASS",]
  
  #separate the alleles
  constants = table[,1:10]
  counts = table[,11:ncol(table)]
  
  bl6 = counts[,grepl("BL6", colnames(counts))]
  s129 = counts[,grepl("S129", colnames(counts))]
  
  bl6 = cbind(constants, bl6)
  s129 = cbind(constants, s129)
  
  #collapse down to gene level
  bl6_collapsed = collapse_duplicates(bl6)
  s129_collapsed = collapse_duplicates(s129)
  
  #remove the chaff
  rownames(bl6_collapsed) = bl6_collapsed$gene_stable_id
  bl6_collapsed = bl6_collapsed[,11:ncol(bl6_collapsed)]
  rownames(s129_collapsed) = s129_collapsed$gene_stable_id
  s129_collapsed = s129_collapsed[,11:ncol(s129_collapsed)]
  
  
  #get cell names
  name = basename(filename)
  split = strsplit(name, "_")[[1]]
  treatment = split[3]
  batch = substr(split[4], 1, 1)
  
  cellnames = names(bl6_collapsed)
  split = strsplit(cellnames, ".", fixed = T)
  cellnum = lapply(split, function(x) x[1])
  cellnum = sub("X", "", cellnum)
  
  full_cellnames = paste0("ola_mES_", treatment, "_", batch, "_", cellnum)
  
  names(bl6_collapsed) = full_cellnames
  names(s129_collapsed) = full_cellnames
  
  return(list(bl6 = bl6_collapsed, s129 = s129_collapsed))
  
}

allele_files = dir(paste0(folder_location, "raw_data/mesc_alleles/"), full.names = T)
allele_files = allele_files[grepl(".transcript", allele_files)]


for(file in allele_files){
  # print(paste("File", which(allele_files==file), "of", length(allele_files)))
  in_file = read.table(file, header = T)
  process = process_mesc_allele(table = in_file, filename = file)
  if(file == allele_files[1]){
    #make objects if doing first run
    bl6_mat = process$bl6
    s129_mat = process$s129
  } else {
    bl6_mat = cbind(bl6_mat, process$bl6)
    s129_mat = cbind(s129_mat, process$s129)
  }
}

#standardise to our mESC count matrix i.e. remove XY chr genes, and
#remove QC fail cells
mesc_bl6 = bl6_mat[rownames(mesc_counts), colnames(mesc_counts)]
mesc_s129 = s129_mat[rownames(mesc_counts), colnames(mesc_counts)]


save(mesc, mesc_meta, mesc_bl6, mesc_s129, mesc_gene_table, file = paste0(folder_location, "/proc_data/mesc_data.RData"))



# LOAD AND PROCESS DENG MOUSE DATA ####

#Each file here contains info for one cell
deng_folder = paste0(folder_location, "raw_data/deng/")
#text file list
files = dir(paste0(deng_folder, "raw"), full.names = T, pattern = "*.txt")
#we read in a dummy file to get the genes we need to generate the chromosome map
temp = read.table(files[1])
#make gene table
deng_gene_table = getBM(attributes = c("mgi_symbol", "chromosome_name"), 
                   filters = "mgi_symbol", 
                   values = temp[,1], mart = 
                     mouse_ensembl)
deng_gene_table$chr = deng_gene_table$chromosome_name
#keep only autosome
deng_gene_table = deng_gene_table[deng_gene_table$chromosome_name %in% 1:19,]
deng_gene_table = deng_gene_table[!duplicated(deng_gene_table$mgi_symbol),]
rownames(deng_gene_table) = deng_gene_table$mgi_symbol

#all the genes are the same for all cells i.e. they have used full database
#so we can prepare the chromosome name vector from one entry just to add on to each other cell
#this vector maps each count to a chromosome for easier summation
chr_loc = sapply(temp[,1], function(x) ifelse(x %in% deng_gene_table$mgi_symbol,
                                              deng_gene_table[min(which(deng_gene_table$mgi_symbol==x)), "chromosome_name"],
                                              NA))

#go through files, reading them in and storing in a big list for lapply later
deng_list = list()
for(file in files){
  temp = read.table(file)
  #relabel
  names(temp) = c("gene", "refseq", "rpkm", "count", "cast", "c57")
  #add ratio cast/c57
  temp$ratio = temp$cast/temp$c57
  #add chromosome names
  temp$chr = chr_loc
  
  deng_list[[length(deng_list)+1]] = temp
  names(deng_list)[[length(deng_list)]] = basename(file)
}


#make meta - extract info from cell names
namesplit = strsplit(names(deng_list), split = "_")
stage = sapply(namesplit, function(x) x[2])
embcell = sapply(namesplit, function(x) x[3])
embsplit = strsplit(embcell, split = "-")
emb = sapply(embsplit, function(x) x[1])
cell = sapply(embsplit, function(x) x[2])

deng_meta = data.frame(
  cell = names(deng_list),
  stage = stage,
  emb_num = emb,
  cell_num = cell,
  embryo = paste(stage, emb),
  row.names = names(deng_list)
)


#drop the cells that are 100% one set of alleles: (non-informative about aneuploidy)
#c57 and also the early 2 cell and zygotes
remove = c("C57twocell", "early2cell", paste0("zy", 1:4), "fibroblast", "BXC")
deng_meta = deng_meta[!deng_meta$stage%in%remove,]
deng_meta$stage = factor(deng_meta$stage)


#get counts matrix
deng_counts = sapply(deng_list, function(x) x$count)
rownames(deng_counts) = deng_list[[1]][,1]

#keep only genes present in the counts 
deng_gene_table = deng_gene_table[deng_gene_table$mgi_symbol%in%rownames(deng_counts),]

deng_counts = deng_counts[deng_gene_table$mgi_symbol, rownames(deng_meta)]

#make deng object
deng = makeAneu(counts = deng_counts, genes = rownames(deng_counts), chrs = deng_gene_table$chromosome_name, cellNames = colnames(deng_counts), cellGroups = as.character(deng_meta$stage))

#make counts matrices with exclusively the reads from the different alleles
cast_counts = sapply(deng_list, function(x) x$cast)
rownames(cast_counts) = deng_list[[1]][,1]
cast_counts = cast_counts[rownames(deng_counts), rownames(deng_meta)]

c57_counts = sapply(deng_list, function(x) x$c57)
rownames(c57_counts) = deng_list[[1]][,1]
c57_counts = c57_counts[rownames(deng_counts), rownames(deng_meta)]



#save the lot
save(deng, deng_gene_table, deng_meta, cast_counts, c57_counts, file = paste0(folder_location, "proc_data/deng_data.RData"))


# LOAD AND PROCESS SUPPLEMENTAL DATA ####
# i.e. those for the method quality assessment.

#first the zeisel umis
umi = readRDS(paste0(folder_location, "raw_data/3prime/brain_data.rds"))
umi_celltypes = umi@phenoData@data$level1class
#take oligodendrocytes
odc_counts = umi@assayData$counts[, umi_celltypes=="oligodendrocytes"]
#normalise
odc = makeAneu(counts = odc_counts, 
               genes = rownames(odc_counts), 
               chrs = rep(NA, nrow(odc_counts)), 
               cellNames = paste0("cell_", 1:ncol(odc_counts)), 
               cellGroups = rep(1, ncol(odc_counts)))

#now the 10X data
ten = readMM(paste0(folder_location, "raw_data/3prime/10x_counts.mtx"))
ten_barcodes = read.table(paste0(folder_location, "raw_data/3prime/10x_barcodes.tsv"))
ten_genes = read.table(paste0(folder_location, "raw_data/3prime/10x_genes.tsv"))
ten_clusters = read.table(paste0(folder_location, "raw_data/3prime/10x_clusters.csv"), sep = ",", header = T)

#table(ten_clusters$Barcode == ten_barcodes$V1) #tick
#subset cluster 2
ten = ten[,which(ten_clusters$Cluster==2)]
colnames(ten) = ten_barcodes[which(ten_clusters$Cluster==2),]
rownames(ten) = ten_genes$V1
#detected genes
ten = ten[rowSums(ten)>0,]
#normalise
ten = makeAneu(counts = as.matrix(ten), 
               genes = rownames(ten), 
               chrs = rep(NA, nrow(ten)), 
               cellNames = colnames(ten), 
               cellGroups = rep(1, ncol(ten)))


#Then the high quality data from the Fluidigm UMI approach (Tung et al)
umi_reads = read.table(paste0(folder_location, "raw_data/3prime/gilad_reads.txt"), header = T, row.names = 1)
umi_reads = makeAneu(counts = as.matrix(umi_reads), 
                     genes = rownames(umi_reads), 
                     chrs = rep(NA, nrow(umi_reads)), 
                     cellNames = colnames(umi_reads), 
                     cellGroups = rep(1, ncol(umi_reads)))

umi_mols = read.table(paste0(folder_location, "raw_data/3prime/gilad_molecules.txt"), header = T, row.names = 1)
umi_mols = makeAneu(counts = as.matrix(umi_mols), 
                     genes = rownames(umi_mols), 
                     chrs = rep(NA, nrow(umi_mols)), 
                     cellNames = colnames(umi_mols), 
                     cellGroups = rep(1, ncol(umi_mols)))

# Cel-seq2 data
cel = read.table(paste0(folder_location, "raw_data/cel_seq.txt"), header = TRUE, row.names = 1, sep = "\t")
sc = strsplit(names(cel), ".", fixed = TRUE)
is.cell = sapply(sapply(sc, function(x) x[1]), function(y) substr(y, start =2, stop = 2))
# Remove doublets etc. i.e. is.cell == 0
cel = cel[, as.logical(as.numeric(is.cell))]

cel = makeAneu(counts = as.matrix(cel), 
               genes = rownames(cel),
               chrs = rep(1, nrow(cel)), 
               cellNames = colnames(cel), 
               cellGroups = rep(1, ncol(cel)))


# LOAD AND PROCESS SCIALDONE DATA ####
as_counts = read.table(paste0(folder_location, "raw_data/scialdone_counts.txt"), header = TRUE, row.names= 1)
as_meta = read.table(paste0(folder_location, "raw_data/scialdone_meta.txt"), header = TRUE)

#take turquoise cluster: epiblast
as_counts = as_counts[, colnames(as_counts) %in% as_meta$cellName[as_meta$cluster=="turquoise"]]


sci_gene_table = getBM(
  attributes = c("ensembl_gene_id", "chromosome_name"),
  mart = mouse_ensembl,
  values = as.character(rownames(as_counts)),
  filters = "ensembl_gene_id"
)

sci_gene_table = sci_gene_table[sci_gene_table$chromosome_name %in% 1:19,]
as_counts = as_counts[rownames(as_counts) %in% sci_gene_table$ensembl_gene_id, ]

#make object
scialdone = makeAneu(counts = as.matrix(as_counts),
                     genes = rownames(as_counts), 
                     chrs = sci_gene_table$chromosome_name[match(rownames(as_counts), 
                                                                 sci_gene_table$ensembl_gene_id)], 
                     cellNames = colnames(as_counts), 
                     cellGroups = rep(1, ncol(as_counts)))

save(odc, ten, umi_reads, umi_mols, cel, scialdone, file = paste0(folder_location, "proc_data/performance.Rdata"))
