# E. coli variant calling (NGS)

End-to-end short-read variant calling on an E. coli sample (SRA: SRR2584863) against the REL606 reference (GCA_000017985.1), with manual IGV review and a Python summary.

## Pipeline
FastQC -> fastp -> bwa-mem2 -> samtools -> bcftools -> snpEff -> IGV -> pandas/matplotlib

Full commands: `scripts/pipeline.sh`. Tool versions I could confirm: samtools 1.13, bcftools 1.13, bwa-mem2 2.2.1, snpEff 5.1.

## Results
- 1,553,259 read pairs; 1,286,105 pairs kept after fastp (83%), Q30 96.0% (R1) and 92.6% (R2)
- 100% of reads mapped, 98.16% properly paired, 99.85% of the genome covered, mean depth about 74x
- 30 raw variants -> 27 after QUAL/depth filter -> 22 after adding MQ>=40
- Manual IGV review of 11 of the 22 calls: 7 real, 2 likely real, 1 probable (tandem-repeat insertion, exact size uncertain), 1 artifact (coverage cliff at `yieO`, removed). The other 11 were not inspected.

Files: `results/final_calls_with_verdicts.csv`, `results/depth_vs_altfrac.png`, `results/S1.final.clean.vcf`, `notebooks/`, `qc/`.

## What I learned
- Chromosome names must match across genome, BAM and VCF (`CP000819` vs `CP000819.1`).
- Repeats and homopolymers produce doubtful calls; mapping quality and IGV inspection catch many.
- Raw calls need filtering and visual checks before they are trusted.
