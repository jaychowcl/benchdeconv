#$ -cwd -V
#$ -N benchdeconv
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G
#$ -l rl9=true
#$ -pe sharedmem 8
#$ -t 1-50
#$ -N ArrayJob
#$ -o ./log/output.txt
#$ -e ./log/errors.txt

#scripthere
. /etc/profile.d/modules.sh
module load anaconda
conda activate new_conda

# Determine downsize value based on task ID
if [ ${SGE_TASK_ID} -le 10 ]; then
  DOWNSIZE=0
elif [ ${SGE_TASK_ID} -le 20 ]; then
  DOWNSIZE=500
elif [ ${SGE_TASK_ID} -le 30 ]; then
  DOWNSIZE=1000
elif [ ${SGE_TASK_ID} -le 40 ]; then
  DOWNSIZE=2500
elif [ ${SGE_TASK_ID} -le 50 ]; then
  DOWNSIZE=5000
else
  echo "Error: Task ID out of range"
  exit 1
fi

Rscript --vanilla ./scripts/run.R --scdata "data/scRNA_wu" \
--scmeta "data/scRNA_wu/metadata.csv" \
--outdir /exports/eddie3_homes_local/s2600569/benchdeconv/data/results/all_runs/run_${SGE_TASK_ID} \
--seed ${SGE_TASK_ID} \
--downsize ${DOWNSIZE} \
--grain_lvl "celltype_major" \
--gene_column 1 \
--synth_dataset "artificial_regional_rare_celltype_diverse" \
--select_celltype "T-cells"


#work in cwd (-V = use global vars)
#give name
#give runtime
#give required mem # remove if not needed, but may wait more for scheduling

#give output outfile 
#give error outfile
