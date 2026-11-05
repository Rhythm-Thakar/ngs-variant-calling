#!/usr/bin/env bash
# E. coli variant calling: SRR2584863 vs REL606 (GCA_000017985.1)
set -euo pipefail
mkdir -p ref raw qc trimmed aln var annot snpeff/data/rel606
REF=ref/GCA_000017985.1_ASM1798v1_genomic.fna

# 1. Download reference and reads
wget -P ref https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip ref/*.fna.gz
(cd raw && prefetch SRR2584863 && fasterq-dump SRR2584863 --split-files)

# 2. QC and trimming
fastqc raw/SRR2584863_1.fastq raw/SRR2584863_2.fastq -o qc -t 2
fastp -i raw/SRR2584863_1.fastq -I raw/SRR2584863_2.fastq \
  -o trimmed/R1.fq.gz -O trimmed/R2.fq.gz \
  --detect_adapter_for_pe --trim_poly_x \
  --cut_right --cut_right_mean_quality 20 \
  --length_required 50 --thread 2 \
  -h trimmed/fastp.html -j trimmed/fastp.json

# 3. Align, sort, index
bwa-mem2 index $REF
samtools faidx $REF
bwa-mem2 mem -t 2 -R '@RG\tID:S1\tSM:S1\tPL:ILLUMINA' $REF \
  trimmed/R1.fq.gz trimmed/R2.fq.gz | samtools sort -@ 2 -o aln/S1.sorted.bam -
samtools index aln/S1.sorted.bam
samtools flagstat aln/S1.sorted.bam
samtools coverage aln/S1.sorted.bam

# 4. Call and filter variants (haploid)
bcftools mpileup -f $REF -Ou aln/S1.sorted.bam \
  | bcftools call -mv --ploidy 1 -Ov -o var/S1.raw.vcf
bcftools filter -i 'QUAL>=30 && DP>=10 && MQ>=40' var/S1.raw.vcf -Ov -o var/S1.final.vcf

# 5. Annotate with a local snpEff database built from the GenBank file
wget -O snpeff/data/rel606/genes.gbk.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.gbff.gz
gunzip snpeff/data/rel606/genes.gbk.gz
cp /etc/snpeff/snpEff.config snpeff/
sed -i "s|^data.dir.*|data.dir = $PWD/snpeff/data/|" snpeff/snpEff.config
echo "rel606.genome : E. coli REL606" >> snpeff/snpEff.config
echo "rel606.CP000819.1.codonTable : Bacterial_and_Plant_Plastid" >> snpeff/snpEff.config
snpEff build -genbank -noCheckCds -c snpeff/snpEff.config rel606
snpEff -Xmx2g -c snpeff/snpEff.config -stats annot/final.html rel606 \
  var/S1.final.vcf > annot/S1.final.annotated.vcf
# remove a stray wrapped header line so bcftools can read the file
awk '/^#CHROM/{h=1} h || /^#/' annot/S1.final.annotated.vcf > annot/S1.final.clean.vcf

# 6. Variant table for the Python notebook
bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/MQ\t%INFO/DP4\t%INFO/ANN\n' annot/S1.final.clean.vcf \
| awk -F'\t' 'BEGIN{OFS="\t"; print "pos","change","effect","impact","gene","DP","MQ","alt_frac"}
{split($7,d,","); tot=d[1]+d[2]+d[3]+d[4]; split($8,b,","); split(b[1],c,"|");
 r=(length($2)>12)?"long":$2; a=(length($3)>12)?"long":$3;
 print $1, r">"a, c[2], c[3], c[4], $5, $6, sprintf("%.2f",(d[3]+d[4])/tot)}' > annot/variant_table.tsv
