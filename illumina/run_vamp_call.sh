

ls ../fastq_aspa_vamp/*_R1_001.fastq.gz  |  parallel -P 4 ./call_zerotol_paired.sh {};
Rscript merge_counts.r *_counts.txt

rm *.txt *.un1 *.un2 *.join *.out aspa_toxicity_call1.tgz R1.fastq.gz R2.fastq.gz
mv raw.rda raw_vamp.rda
