#$ -cwd -V
#$ -N benchdeconv
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -t 1-10
#$ -N ArrayJob


#scripthere
conda activate new_conda

./scripts/run.R --scdata "data/scRNA_wu" \
--scmeta "data/scRNA_wu/metadata.csv" \
--outdir /data/results/run_${SGE_TASK_ID} \
--seed ${SGE_TASK_ID} \
--test 1 \
--grain-lvl "celltype_major" \
--gene-column 1 \
--synth-dataset "artificial_regional_rare_celltype_diverse"


#$ -o output.o
#$ -e errors.e


#work in cwd (-V = use global vars)
#give name
#give runtime
#give required mem # remove if not needed, but may wait more for scheduling

#give output outfile 
#give error outfile
