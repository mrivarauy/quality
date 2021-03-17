#!/usr/bin/env bash

## EDITAR ESTA SECCIÓN

PREFIX=r10-3
ASSEMBLY=raw-${PREFIX}.fasta
NPROC=32
READS_NANO=/media3/lucia/zymo_mock_community/Zymo-GridION-EVEN-3Peaks-R103-merged.fq.gz
MEDAKA_MODEL=r103_min_high_g345
MP_PARAMS=/media3/lucia/zymo_mock_community/marginPolish/params/allParams.np.microbial.r103g324.json
HELEN_MODEL=/media3/lucia/zymo_mock_community/helen-models/HELEN_r103_guppy_microbial.pkl
REF_DIR_ILL=/media3/lucia/zymo_mock_community/illumina_ref_asm/bac/
METAQUAST_OUTDIR=r10
## flags, cambiar a false lo que no quiera correr
flyeflag=true
raconflag=true
medakaflag=true
marginflag=true
helenflag=true
metaquastflag=true
buscoflag=false

# -----------------------------

##assembly con flye
if [ $flyeflag == true ]; then
	flye --nano-raw ${READS_NANO} --meta --iterations 0 --threads ${NPROC} --out-dir assembly-${PREFIX}
	infoseq --only --length --name assembly-${PREFIX}/assembly.fasta | awk '{print $2, $1}' | 
		sort -nr | head -n8 | cut -d' ' -f2 > lista_contigs_mayores
	fastaUtils.pl -u assembly-${PREFIX}/assembly.fasta | 
		grep --no-group-separator -A 1 -w -f lista_contigs_mayores > raw-${PREFIX}.fasta
	ASSEMBLY=raw-${PREFIX}.fasta
fi

##polishing con racon
if [ $raconflag == true ]; then
    minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_NANO} > map1-${PREFIX}.sam
    racon -t ${NPROC} ${READS_NANO} map1-${PREFIX}.sam ${ASSEMBLY} > racon_r1-${PREFIX}.fasta
    minimap2 -ax map-ont -t ${NPROC} racon_r1-${PREFIX}.fasta ${READS_NANO} > map2-${PREFIX}.sam
    racon -t ${NPROC} ${READS_NANO} map2-${PREFIX}.sam racon_r1-${PREFIX}.fasta > racon_r2-${PREFIX}.fasta
fi

## polishing con medaka
if [ $medakaflag == true ]; then
    medaka_consensus -i ${READS_NANO} -d racon_r1-${PREFIX}.fasta -o r1-medaka-${PREFIX} -m ${MEDAKA_MODEL} -t ${NPROC}
    mv r1-medaka-${PREFIX}/consensus.fasta r1-medaka-${PREFIX}.fasta && rm -rf r1-medaka-${PREFIX}
    medaka_consensus -i ${READS_NANO} -d racon_r2-${PREFIX}.fasta -o r2-medaka-${PREFIX} -m ${MEDAKA_MODEL} -t ${NPROC}
    mv r2-medaka-${PREFIX}/consensus.fasta r2-medaka-${PREFIX}.fasta && rm -rf r2-medaka-${PREFIX}
fi

##polishing con marginPolish
if [ $marginflag == true ]; then
	if [ -f map1-${PREFIX}.sam ]; then
	    samtools view -T ${ASSEMBLY} -F 2308 -Sb map1-${PREFIX}.sam | samtools sort -@ ${NPROC} - -o ${ASSEMBLY}.bam
	else
		minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_NANO} | samtools view -T ${ASSEMBLY} -F 2308 -b - |
			samtools sort -@ ${NPROC} -o ${ASSEMBLY}.bam
	fi
	samtools index -@ ${NPROC} ${ASSEMBLY}.bam
	mkdir MP-${PREFIX}
	marginPolish ${ASSEMBLY}.bam ${ASSEMBLY} ${MP_PARAMS} -t ${NPROC} -o MP-${PREFIX} -f
	mv MP-${PREFIX}/output.fa MP-${PREFIX}.fasta
fi

##polishing con helen (requiere polishing con marginPolish)
if [ $helenflag == true ]; then
	cp ${HELEN_MODEL} $(pwd)
	model=${HELEN_MODEL##*/}
	docker run --rm -it --ipc=host -v $(pwd):/helen/ kishwars/helen:latest \
		helen polish -i /helen/MP-${PREFIX}/ -m /helen/${model} -b 256 -w 4 -t 8 -o /helen/ -p MP-helen-${PREFIX}
	mv MP-helen-${PREFIX}.fa MP-helen-${PREFIX}.fasta
fi

##evaluacion con metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o metaquast-${METAQUAST_OUTDIR} *.fasta
fi
