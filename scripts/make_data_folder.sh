#!/bin/bash

# ** Not a literal script **
# The directories have to be created manually for now. Use this for reference.

# Make archive
tar -zcvf chromatin_information_manuscript_data.tar.gz data



# Scripts
cp /home/albanus/manuscripts/atac-seq_information/micro_analyses/vplots_chip_sameN/downsample2.R \
    scripts/
cp /home/albanus/manuscripts/atac-seq_information/micro_analyses/vplots_chip_sameN/downsample.R \
    scripts/
cp /home/albanus/manuscripts/atac-seq_information/micro_analyses/vplots_chip_sameN/config.yml \
    rules/same_n.yml
cp /home/albanus/manuscripts/atac-seq_information/micro_analyses/vplots_chip_sameN/Snakefile \
    rules/same_n.smk
cp ~/scripts/get_ROC_all.R bin/
cp ~/utils/midpointBed/midpointBed /bin/
cp ~/utils/midpointBed/midpointBed bin/
cp /home/albanus/analyses/2019/atacseq_peaks_macs2/bin/calculate_F1_noPIQ.R \
    bin/calculate_F1.R

# Data
cp /lab/work/albanus/2019_encode_chipseq_redo/associations_top1_2019_hepg2.txt data/chipseq_associations_hepg2.txt
cp /lab/work/albanus/2019_encode_chipseq_redo/associations_top1_2019_gm12878.txt data/chipseq_associations_gm12878.txt
cp /lab/work/albanus/2019_hepg2_atac/work/tf-binding/output/hepg2_{1,3}.fvice.chip.out data/processed/tf-binding/hepg2/
cp /lab/work/albanus/2018_redo_all/work/f-vices_1KG/output/lab_gm12878__information_track_clusters.txt \
    data/processed/fvices/
cp /lab/work/albanus/k-mer_restimes/bed_files_sorted_noBL/6-mers/*.bed \
    data/6-mers_scans/ && \
    gzip data/6-mers_scans/*.bed 
cp /lab/work/albanus/2018_redo_all/work/tf-binding/output/lab_gm12878.new_fvice.out \
    data/processed/tf-binding/gm12878/chip_seq_fvices.lab_gm12878.txt
cp /lab/work/albanus/2018_redo_all/work/tf-binding/output/buenrostro_rep1.new_fvice.out \
    data/processed/tf-binding/gm12878/chip_seq_fvices.buenrostro_rep1.txt
cp ~/data/motifs/motifInformationTotalBits.dat data/motifInformationTotalBits.dat
cp /lab/work/albanus/2018_redo_all/work/f-vices_1KG/output/protein_domains_lambert_bonf0.05.RData \
    data/protein_domains/
cp /lab/work/albanus/2019_eqtl_enrichments/work/f-vices_1KG/gregor_eqtl_redo_same_n/output/*  \
    data/processed/eqtl/

# NGS plot
indir="/lab/work/albanus/2018_gm12878_mnase-seq/work/ngsplot_per_motif/\
lab_gm12878"
outdir="data/processed/mnase/lab_gm12878"
mkdir -p ${outdir}
for f in `ls ${indir}/SRR452483__*/avgprof.RData`
do
    base=`echo ${f} | sed 's/.*\/\(SRR\)/\1/g; s/\//./g'`
    cp ${f} ${outdir}/${base}
done

indir="/lab/work/albanus/2018_gm12878_mnase-seq/work/ctcf_cohesin/ngsplot/\
lab_gm12878/ext_500"
outdir="data/processed/mnase/ctcf_cohesin"
mkdir -p ${outdir}
for f in `ls ${indir}/SRR452483.CTCF*/avgprof.RData`
do
    base=`echo ${f} | sed 's/.*\/\(SRR\)/\1/g; s/\//./g'`
    cp ${f} ${outdir}/${base}
done

# GREGOR neiborhood files
outdir="data/processed/eqtl/neighborhood"
mkdir -p ${output}
glob_ls="/lab/work/albanus/2019_eqtl_enrichments/work/f-vices_1KG/\
gregor_eqtl_redo_same_n/*/output*/*.bed/neighborhoodFile.txt"
for f in `ls ${glob_ls} | grep -v "dbd"`
do
    base=`echo ${f} | sed 's/.*gregor_eqtl_redo_same_n\///g'`
    base=`echo ${base} | sed 's/\//__/g'`
    cp ${f} ${outdir}/${base}
done

# BMO predictions
samples="lab_gm12878 buenrostro_rep1 abcu196_4 alpha3 beta3 cd4_87 \
cd4_88 hepg2_1 hepg2_3"
for i in ${samples}
do
    indir="/lab/work/albanus/2018_redo_all/work/f-vices_1KG/bmo/${i}/bound"
    outdir="data/processed/bmo_predictions/${i}"
    mkdir -p ${outdir}
    for f in `ls ${indir}/*.bed`
    do
        base=`basename ${f}`
        out="${outdir}/${base}.gz"
        gzip -c ${f} > ${out}
    done
done

## V-plots
samples="lab_gm12878 buenrostro_rep1 abcu196_4 alpha3 beta3 cd4_87 \
cd4_88 hepg2_1 hepg2_3"
for i in ${samples}
do
    indir="/lab/work/albanus/2019_vplots_new/work/f-vices_1KG/vplots/${i}"
    outdir="data/processed/chromatin_information/${i}"
    mkdir -p ${outdir}
    for f in `ls ${indir}/*.posinfo.gz`
    do
        base=`basename ${f}`
        out="${outdir}/${base}"
        cp ${f} ${out}
    done
done

indir="/lab/work/albanus/2018_cohesin_and_tfs/work/vplots_same_n_2019/\
lab_gm12878"
outdir="data/processed/ctcf_cohesin/lab_gm12878"
mkdir -p ${outdir}
for f in `ls ${indir}/*.posinfo.gz`
do
    base=`basename ${f}`
    out="${outdir}/${base}"
    cp ${f} ${out}
done

indir="/lab/work/albanus/2018_cohesin_and_tfs/work/vplots_same_n_2019/\
lab_gm12878"
outdir="data/processed/ctcf_cohesin/lab_gm12878"
mkdir -p ${outdir}
for f in `ls ${indir}/*.posinfo.gz`
do
    base=`basename ${f}`
    out="${outdir}/${base}"
    cp ${f} ${out}
done

samples="lab_gm12878_sonicated_3 lab_gm12878_subsampled_3"
for i in ${samples}
do
    indir="/lab/work/albanus/2018_cohesin_and_tfs/work/vplots_same_n_manual_redo2/${i}"
    outdir="data/processed/ctcf_cohesin/${i}"
    mkdir -p ${outdir}
    for f in `ls ${indir}/*.posinfo.gz`
    do
        base=`basename ${f}`
        out="${outdir}/${base}"
        cp ${f} ${out}
    done
done
mv data/processed/ctcf_cohesin/lab_gm12878_subsampled_3 \
    data/processed/ctcf_cohesin/lab_gm12878_downsampled

# Asymmetry
# indir="/lab/work/albanus/2018_redo_all/work/f-vices_1KG/vplot_asymmetry/lab_gm12878"
# outdir="data/processed/asymmetry/permutations"
# mkdir -p ${outdir}
# while read motif
# do
#     cp ${indir}/${motif}.everything.RData ${outdir}
# done < data/motif_list.trimmed.txt

outdir="data/processed/asymmetry"
cp /lab/work/albanus/2018_redo_all/work/f-vices_1KG/vplot_asymmetry/*.RData ${outdir}
cp /lab/work/albanus/2018_redo_all/work/f-vices_1KG/vplot_asymmetry/sides_prox_dista_lab_gm12878.RData ${outdir}
cp /lab/work/albanus/2018_redo_all/work/f-vices_1KG/vplot_asymmetry_cage10/info_prox_distal_lab_gm12878.RData ${outdir}

# Motif pwms
outdir="data/pwm"
mkdir ${outdir}
for f in `ls /lab/data/motifs/pwm/ENCODE2013/hg18_matrixes/*.mat | grep -v "_disc"`
do
    cp ${f} ${outdir}
done
for f in `ls /lab/data/motifs/pwm/JASPAR2014/hg18_matrixes/*.mat`
do
    cp ${f} ${outdir}
done


