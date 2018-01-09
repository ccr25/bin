#!/bin/bash

#SBATCH --mem=80000
#SBATCH --cpus-per-task=6
#SBATCH --time=36:00:00
module load java/1.8.0_20

filename=${1}
reference=${2}
filebase=${filename%%.*}


jid1=$(srun java -jar /packages/gatk/3.3-0/GenomeAnalysisTK.jar -T RealignerTargetCreator -o ${filebase}.intervals -I $filename -R $reference)

jid2=$(srun --dependency=afterok:$jid1 java -jar /packages/gatk/3.3-0/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filename -targetIntervals ${filebase}.intervals -o ${filebase}_realigned.bam)

jid3=$(srun --dependency=afterok:$jid2 java -jar /packages/gatk/3.3-0/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I ${filebase}_realigned.bam -o ${filebase}.vcf -hets 0.01 -out_mode EMIT_ALL_CONFIDENT_SITES -ploidy 1)

jid4=$(srun --dependency=afterok:$jid3 java -jar /packages/gatk/3.3-0/GenomeAnalysisTK.jar -T VariantFiltration -S silent -R $reference -o ${filebase}_filtered.vcf --variant ${filebase}.vcf --filterExpression "QD = 2.0 || FS_filter = 60.0 || MQ_filter = 30.0 || MQ_Rank_Sum_filter = -12.5 || Read_Pos_Rank_Sum_filter = -8" --filterName ${filebase}_filterName)
