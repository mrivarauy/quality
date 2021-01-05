#!/usr/bin/env bash

## EDITAR ESTA SECCIÓN

CONDA=/media3/anaconda3/etc/profile.d/conda.sh
ASSEMBLY=flye.fasta
PREFIX=r10
NPROC=28
READS_NANO=/media3/lucia/zymo_mock_community/Zymo-GridION-EVEN-3Peaks-R103-merged.fq.gz
MEDAKA_MODEL=r103_min_high_g345
MP_PARAMS=/media3/lucia/zymo_mock_community/marginPolish/params/allParams.np.microbial.r103g324.json
HELEN_MODEL=/media3/lucia/zymo_mock_community/helen-models/HELEN_r103_guppy_microbial.pkl
HMP_DIR=/media3/lucia/zymo_mock_community/hmp/
HOMO_MASH=/media3/lucia/zymo_mock_community/hmp/bacteria.msh
HOMO_MODEL=/media3/lucia/zymo_mock_community/hmp/R10.3.pkl
REF_DIR_ILL=/media3/lucia/zymo_mock_community/illumina_ref_asm/bac/
METAQUAST_OUTDIR=ill
## flags, cambiar a false lo que no quiera correr
raconflag=true
medakaflag=true
marginflag=false
homopolishflag=false
metaquastflag=true

# -----------------------------

source ${CONDA}

##polishing con racon
if [ $raconflag == true ]; then
    minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_NANO} > ${PREFIX}-map1.sam
    racon -t ${NPROC} ${READS_NANO} ${PREFIX}-map1.sam ${ASSEMBLY} > ${PREFIX}-racon_r1.fasta
    minimap2 -ax map-ont -t ${NPROC} ${PREFIX}-racon_r1.fasta ${READS_NANO} > ${PREFIX}-map2.sam
    racon -t ${NPROC} ${READS_NANO} ${PREFIX}-map2.sam ${PREFIX}-racon_r1.fasta > ${PREFIX}-racon_r2.fasta
    minimap2 -ax map-ont -t ${NPROC} ${PREFIX}-racon_r2.fasta ${READS_NANO} > ${PREFIX}-map3.sam
    racon -t ${NPROC} ${READS_NANO} ${PREFIX}-map3.sam ${PREFIX}-racon_r2.fasta > ${PREFIX}-racon_r3.fasta
    minimap2 -ax map-ont -t ${NPROC} ${PREFIX}-racon_r3.fasta ${READS_NANO} > ${PREFIX}-map4.sam
    racon -t ${NPROC} ${READS_NANO} ${PREFIX}-map4.sam ${PREFIX}-racon_r3.fasta > ${PREFIX}-racon_r4.fasta
fi

## polishing con medaka
if [ $medakaflag == true ]; then
    medaka_consensus -i ${READS_NANO} -d ${ASSEMBLY} -o ${PREFIX}-medaka -m ${MEDAKA_MODEL} -t ${NPROC}
    mv ${PREFIX}-medaka/consensus.fasta ${PREFIX}-medaka.fasta
    medaka_consensus -i ${READS_NANO} -d ${PREFIX}-racon_r1.fasta -o ${PREFIX}-r1-medaka -m ${MEDAKA_MODEL} -t ${NPROC}
    mv ${PREFIX}-r1-medaka/consensus.fasta ${PREFIX}-r1-medaka.fasta
    medaka_consensus -i ${READS_NANO} -d ${PREFIX}-racon_r2.fasta -o ${PREFIX}-r2-medaka -m ${MEDAKA_MODEL} -t ${NPROC}
    mv ${PREFIX}-r2-medaka/consensus.fasta ${PREFIX}-r2-medaka.fasta
    medaka_consensus -i ${READS_NANO} -d ${PREFIX}-racon_r4.fasta -o ${PREFIX}-r4-medaka -m ${MEDAKA_MODEL} -t ${NPROC}
    mv ${PREFIX}-r4-medaka/consensus.fasta ${PREFIX}-r4-medaka.fasta
fi

##polishing con marginPolish
if [ $marginflag == true ]; then
    samtools view -@ ${NPROC} -T ${ASSEMBLY} -F 2308 -Sb ${PREFIX}-map1.sam | samtools sort -@ ${NPROC} - -o ${ASSEMBLY}.bam
    samtools index ${ASSEMBLY}.bam
    mkdir ${PREFIX}-margin
    marginPolish ${ASSEMBLY}.bam ${ASSEMBLY} ${MP_PARAMS} -t ${NPROC} -o ${PREFIX}-margin -f
    helen polish -i ${PREFIX}-margin/ -m ${HELEN_MODEL} -b 256 -w 4 -t 8 -o ./ -p ${PREFIX}-margin-helen
    mv ${PREFIX}-margin/output.fa ${PREFIX}-margin.fasta
    mv ${PREFIX}-margin-helen.fa ${PREFIX}-margin-helen.fasta
fi

##polishing con homopolish
if [ $homopolishflag == true ]; then
    conda activate ${HMP_DIR}
    homopolish polish -a ${PREFIX}-racon_r1.fasta -s ${HOMO_MASH} -m ${HOMO_MODEL} -o ${PREFIX}-r1-homopolish -t ${NPROC}
    homopolish polish -a ${PREFIX}-racon_r2.fasta -s ${HOMO_MASH} -m ${HOMO_MODEL} -o ${PREFIX}-r2-homopolish -t ${NPROC}
    homopolish polish -a ${PREFIX}-r1-medaka.fasta -s ${HOMO_MASH} -m ${HOMO_MODEL} -o ${PREFIX}-r1-medaka-homopolish -t ${NPROC}
    homopolish polish -a ${PREFIX}-r2-medaka.fasta -s ${HOMO_MASH} -m ${HOMO_MODEL} -o ${PREFIX}-r2-medaka-homopolish -t ${NPROC}
    homopolish polish -a ${PREFIX}-margin.fasta -s ${HOMO_MASH} -m ${HOMO_MODEL} -o ${PREFIX}-margin-homopolish -t ${NPROC}
    homopolish polish -a ${PREFIX}-margin-helen.fasta -s ${HOMO_MASH} -m ${HOMO_MODEL} -o ${PREFIX}-margin-helen-homopolish -t ${NPROC}
    mv *homopolish*/*.fasta .
fi

##evaluacion con metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --fragmented --min-identity 90 --min-contig 5000 --threads ${NPROC} -r ${REF_DIR_ILL} -o metaquast-${METAQUAST_OUTDIR} *.fasta
fi
