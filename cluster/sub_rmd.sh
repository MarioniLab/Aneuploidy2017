#INITIALISE FOLDERS
my_folder=/nfs/research2/marioni/jonny
out_folder=${my_folder}/clust_out
err_folder=${my_folder}/clust_err

#SELECT SCRIPT
script_name=aneu_rmd
script_to_run=${my_folder}/${script_name}

#CHOOSE PARAMETERS
#RAM in megabytes
memory=100000
r_command="rusage[mem=${memory}]"
#num_processors
nproc=12



bsub -e ${err_folder}/${script_name} \
-o ${out_folder}/${script_name} \
-M $memory -R $r_command -n $nproc -J ${script_name} \
"Rscript /nfs/research2/marioni/jonny/Aneuploidy2017/cluster/run_rmd.R"
