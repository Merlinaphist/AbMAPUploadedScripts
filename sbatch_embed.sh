#!/bin/bash
#SBATCH --job-name=CMV_DENGUE_EMBED  # Job name
#SBATCH --output=%j_output.log  # Output file (with job ID)
#SBATCH --error=%j_error.log    # Error file (with job ID)
#SBATCH --ntasks=1                     # Number of tasks (usually set to 1 for serial jobs)
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=16G                      # Memory per node (e.g., 16GB)
#SBATCH --partition=singhlab-gpu       # Partition name
#SBATCH --gres=gpu:1                   # Request 1 GPU

echo `date`": program starts"

cd /hpc/home/jm688/projects/AbMAP_JM/AbMAPUploadedScripts
source /hpc/group/singhlab/tools/conda/miniconda3/etc/profile.d/conda.sh
conda activate abmap_jm

echo `date`": define variables"
FASTA_PATH=/hpc/home/jm688/projects/AbMAP_JM/diseases_heavy/FASTA
export IgBLASTwd=/hpc/home/jm688/tools/ncbi-igblast-1.22.0
# export IGDATA=/hpc/home/jm688/tools/ncbi-igblast-1.22.0

for n_sample in 1000 10000
do
echo `date`": "$n_sample" sample sequences"
python sample.py -n $n_sample

echo `date`": "$n_sample" Run BLAST+ for alignment"
HEADERS="qseqid\tqstart\tqend\tsseqid\tsstart\tsend\tevalue\tpident\tqframe\tsframe"
for DISEASE in "DEN" "CMV"
do
# for REGION in "V" "J"
# do
blastx -query $FASTA_PATH/${DISEASE}_${n_sample}.fasta \
        -db $IgBLASTwd/database/IGHV/IGHV_AA.fasta \
        -num_alignments 1 -max_hsps 1 \
        -outfmt "6 qseqid qstart qend sseqid sstart send evalue pident qframe sframe" \
        -out output.tsv
echo -e $HEADERS | cat - output.tsv > $FASTA_PATH/${DISEASE}_${n_sample}_BLASTX_IGHV.tsv
rm output.tsv

blastn -query $FASTA_PATH/${DISEASE}_${n_sample}.fasta \
        -task "blastn" \
        -db $IgBLASTwd/database/IGHJ/Edited_IGHJ.fasta \
        -num_alignments 1 -max_hsps 1 \
        -outfmt "6 qseqid qstart qend sseqid sstart send evalue pident qframe sframe" \
        -out output.tsv
echo -e $HEADERS | cat - output.tsv > $FASTA_PATH/${DISEASE}_${n_sample}_BLASTN_IGHJ.tsv
rm output.tsv

# done
done

echo `date`": "$n_sample" impute"
python impute.py -n $n_sample

echo `date`": "$n_sample" Re-run IgBLAST for Translation and Quality Check"
for DISEASE in "DEN" "CMV"
do
$IgBLASTwd/bin/igblastn -germline_db_V $IgBLASTwd/database/IGHV/Edited_IGHV.fasta \
            -germline_db_J $IgBLASTwd/database/IGHJ/Edited_IGHJ.fasta \
            -germline_db_D $IgBLASTwd/database/IGHD/Edited_IGHD.fasta \
            -c_region_db $IgBLASTwd/database/c_genes/ncbi_human_c_genes \
            -auxiliary_data $IgBLASTwd/optional_file/human_gl.aux \
            -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 \
            -query $FASTA_PATH/${DISEASE}_${n_sample}_imputed.fasta \
            -outfmt 19 -out $FASTA_PATH/${DISEASE}_${n_sample}_imputed_IgBLAST.tsv \
            -show_translation
done

python post_quality.py -n $n_sample

echo `date`": "$n_sample" AbMAP for embedding"
python embed.py -n $n_sample

done

echo `date`": UMAP for visualization"
python plot_embed.py

echo `date`": program ends"
