# construct ref
# 1.download fastq and gtf
cd ref
aria2c https://ftp.ensembl.org/pub/release-109/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
gunzip Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
aira2c https://ftp.ensembl.org/pub/release-109/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.109.gtf.gz
gunzip Macaca_mulatta.Mmul_10.109.gtf.gz
# 2.filter gtf
cellranger mkgtf \
  Macaca_mulatta.Mmul_10.109.gtf \
  Macaca_mulatta.Mmul_10.109.filtered.gtf \
  --attribute=gene_biotype:protein_coding
cellranger mkref \
  --genome=Macaca_mulatta_genome \
  --fasta=Macaca_mulatta.Mmul_10.dna.toplevel.fa \
  --genes=Macaca_mulatta.Mmul_10.109.filtered.gtf 

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
# 3.cellranger流程
  cellranger count --id="results_"${i%*_S01_L003-001.fastq.tar}   --transcriptome=ref/Macaca_mulatta_genome/   --fastqs=fastq2   --sample=${i%*_S01_L003-001.fastq.tar}   --localcores=20 --nosecondary
# 4.将所有outs提取出来  
  mkdir "results/"${i%*_S01_L003-001.fastq.tar}
  mv "results_"${i%*_S01_L003-001.fastq.tar}"/outs/filtered_feature_bc_matrix.h5" "results/"${i%*_S01_L003-001.fastq.tar}"/filtered_feature_bc_matrix.h5"
  mv "results_"${i%*_S01_L003-001.fastq.tar}"/outs/filtered_feature_bc_matrix" "results/"${i%*_S01_L003-001.fastq.tar}"/filtered_feature_bc_matrix"
  rm -rf "results_"${i%*_S01_L003-001.fastq.tar}
done