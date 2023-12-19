# construct ref
# 1.download fastq and gtf
cd ref
aria2c https://ftp.ensembl.org/pub/release-109/fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa.gz
gunzip Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa.gz
aira2c https://ftp.ensembl.org/pub/release-109/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.109.gtf.gz
gunzip Pan_troglodytes.Pan_tro_3.0.109.gtf.gz
# 2.filter gtf
cellranger mkgtf \
  Pan_troglodytes.Pan_tro_3.0.109 \
  Pan_troglodytes.Pan_tro_3.0.109.filtered.gtf \
  --attribute=gene_biotype:protein_coding
cellranger mkref \
  --genome=Chimpanzee_genome \
  --fasta=Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa \
  --genes=Pan_troglodytes.Pan_tro_3.0.109.filtered.gtf 

# download raw data(fastq)
cd ../fastq/
aria2c -c -i download_urls.txt

# run cellranger
cd ..
for i in $(ls fastq/)
do
  mkdir "fastq2/"${i%*_S01_L003-001.fastq.tar} 
  tar -xvf "fastq/"$i -C "fastq2/"${i%*_S01_L003-001.fastq.tar}
  rm "fastq/"$i
  mv "fastq2/"${i%*_S01_L003-001.fastq.tar}"/"${i%*.fastq.tar}/* "fastq2/"${i%*_S01_L003-001.fastq.tar}
  rm -rf "fastq2/"${i%*_S01_L003-001.fastq.tar}"/"${i%*.fastq.tar}
  rename S01 S1 "fastq2/"${i%*_S01_L003-001.fastq.tar}/*
# 3.cellranger pipeline
  cellranger count --id="results_"${i%*_S01_L003-001.fastq.tar}   --transcriptome=ref/Chimpanzee_genome/   --fastqs=fastq2   --sample=${i%*_S01_L003-001.fastq.tar}   --localcores=20 --nosecondary
# 4.extract all samples
  mkdir "results/"${i%*_S01_L003-001.fastq.tar}
  mv "results_"${i%*_S01_L003-001.fastq.tar}"/outs/filtered_feature_bc_matrix.h5" "results/"${i%*_S01_L003-001.fastq.tar}"/filtered_feature_bc_matrix.h5"
  mv "results_"${i%*_S01_L003-001.fastq.tar}"/outs/filtered_feature_bc_matrix" "results/"${i%*_S01_L003-001.fastq.tar}"/filtered_feature_bc_matrix"
  mv "results_"${i%*_S01_L003-001.fastq.tar}"/outs/molecule_info.h5" "results/"${i%*_S01_L003-001.fastq.tar}"/molecule_info.h5"
  rm -rf "results_"${i%*_S01_L003-001.fastq.tar}
done

# construct cellranger aggragate table
ls results/ >> cellranger_aggr_table.csv
paste -d ',' cellranger_aggr_table.csv cellranger_aggr_table.csv >> cellranger_aggr_table2.csv
for i in {1..26};do echo "molecule_info.h5";done >> test.csv
paste -d '/' cellranger_aggr_table2.csv test.csv >cellranger_aggr_table.csv
# run cellranger aggr
cellranger aggr --id=aggr_chimpanzee --csv cellranger_aggr_table.csv --nosecondary  --localcores 20