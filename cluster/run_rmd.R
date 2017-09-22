library(rmarkdown)

fout = render("/nfs/research2/marioni/jonny/Aneuploidy2017/PaperSupplementary.Rmd", run_pandoc = TRUE, clean = FALSE)
print(paste("Output at", fout))
