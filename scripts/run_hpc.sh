#$ -cwd -V
#$ -N benchdeconv
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -t 1-10
#$ -N ArrayJob


#scripthere
./scripts/run.R --scdata "~/project/data/scRNA_wu" --scmeta "/localdisk/home/s2600569/project/data/scRNA_wu/metadata.csv" --outdir ./data/results/run_${SGE_TASK_ID}






#$ -o output.o
#$ -e errors.e


#work in cwd (-V = use global vars)
#give name
#give runtime
#give required mem # remove if not needed, but may wait more for scheduling




#give output outfile 
#give error outfile
