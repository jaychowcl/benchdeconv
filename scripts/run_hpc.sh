#$ -cwd -V
#$ -N benchdeconv
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G
#$ -l rl9=true
#$ -pe sharedmem 8
#$ -t 1-180
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
  mintest=0
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 20 ]; then
  DOWNSIZE=500
  mintest=0
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 30 ]; then
  DOWNSIZE=1000
  mintest=0
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 40 ]; then
  DOWNSIZE=2500
  mintest=0
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1.csv"  
elif [ ${SGE_TASK_ID} -le 50 ]; then
  DOWNSIZE=5000
  mintest=0
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 60 ]; then
  DOWNSIZE=5000
  mintest=0
  mintest_id="T-cells"
  subtype="HER2+"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 70 ]; then
  DOWNSIZE=5000
  mintest=0
  mintest_id="T-cells"
  subtype="ER+"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 80 ]; then
  DOWNSIZE=5000
  mintest=0
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1.csv"
elif [ ${SGE_TASK_ID} -le 85 ]; then
  DOWNSIZE=5000
  mintest=0.05
  mintest_id="B-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 90 ]; then
  DOWNSIZE=5000
  mintest=0.1
  mintest_id="B-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 95 ]; then
  DOWNSIZE=5000
  mintest=0.2
  mintest_id="B-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 100 ]; then
  DOWNSIZE=5000
  mintest=0.5
  mintest_id="B-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 105 ]; then
  DOWNSIZE=5000
  mintest=0.8
  mintest_id="B-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 110 ]; then
  DOWNSIZE=5000
  mintest=0.05
  mintest_id="Plasmablasts"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 115 ]; then
  DOWNSIZE=5000
  mintest=0.1
  mintest_id="Plasmablasts"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 120 ]; then
  DOWNSIZE=5000
  mintest=0.2
  mintest_id="Plasmablasts"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 125 ]; then
  DOWNSIZE=5000
  mintest=0.5
  mintest_id="Plasmablasts"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 130 ]; then
  DOWNSIZE=5000
  mintest=0.8
  mintest_id="Plasmablasts"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 135 ]; then
  DOWNSIZE=5000
  mintest=0.05
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 140 ]; then
  DOWNSIZE=5000
  mintest=0.1
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 145 ]; then
  DOWNSIZE=5000
  mintest=0.2
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 150 ]; then
  DOWNSIZE=5000
  mintest=0.5
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 155 ]; then
  DOWNSIZE=5000
  mintest=0.8
  mintest_id="T-cells"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
  elif [ ${SGE_TASK_ID} -le 160 ]; then
  DOWNSIZE=5000
  mintest=0.05
  mintest_id="Myeloid"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 165 ]; then
  DOWNSIZE=5000
  mintest=0.1
  mintest_id="Myeloid"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 170 ]; then
  DOWNSIZE=5000
  mintest=0.2
  mintest_id="Myeloid"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 175 ]; then
  DOWNSIZE=5000
  mintest=0.5
  mintest_id="Myeloid"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
elif [ ${SGE_TASK_ID} -le 180 ]; then
  DOWNSIZE=5000
  mintest=0.8
  mintest_id="Myeloid"
  subtype="TNBC"
  coords="./data/spot_coords/out1_mintest.csv"
else
  echo "Error: Task ID out of range"
  exit 1
fi

Rscript --vanilla ./scripts/run.R --scdata "data/scRNA_wu" \
--scmeta "data/scRNA_wu/metadata.csv" \
--outdir /exports/eddie/scratch/s2600569/benchdeconv_results/all_runs/run_${SGE_TASK_ID} \
--seed ${SGE_TASK_ID} \
--downsize ${DOWNSIZE} \
--grain_lvl "celltype_major" \
--gene_column 1 \
--synth_dataset "artificial_diverse_overlap" \
--select_celltype "T-cells" \
--n_cells_max 15 \
--min_cell_id_test ${mintest} \
--select_celltype_min_id ${mintest_id} \
--subtype ${subtype} \
--coords ${coords} \
--coords_total "./data/spot_coords/Spatial-Projection.csv" 


#work in cwd (-V = use global vars)
#give name
#give runtime
#give required mem # remove if not needed, but may wait more for scheduling

#give output outfile 
#give error outfile
