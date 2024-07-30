#$ -cwd -V
#$ -N benchdeconv
#$ -l h_rt=2:00:00
#$ -l h_vmem=16G
#$ -l rl9=true
#$ -pe sharedmem 8
#$ -t 1-2
#$ -N ArrayJob
#$ -o output.txt
#$ -e errors.txt

#scripthere
. /etc/profile.d/modules.sh
module load anaconda
conda activate new_conda

Rscript --vanilla ./scripts/run.R --scdata "data/scRNA_wu" \
--scmeta "data/scRNA_wu/metadata.csv" \
--outdir /exports/eddie3_homes_local/s2600569/benchdeconv/data/results/all_runs/run_${SGE_TASK_ID} \
--seed ${SGE_TASK_ID} \
--test 1 \
--grain_lvl "celltype_major" \
--gene_column 1 \
--synth_dataset "artificial_regional_rare_celltype_diverse"


#work in cwd (-V = use global vars)
#give name
#give runtime
#give required mem # remove if not needed, but may wait more for scheduling

#give output outfile 
#give error outfile
