#/bin/bash

echo "Downsampling signals..."

# Downsample V-plots to make them directly comparable
Rscript bin/downsample_vsignals.R

echo "Generating V-plots..."

indir="work/vplots_sonicated"
outdir="${indir}"
motif="CTCF__CTCF_known2__ENCFF963PJY-CTCF__ENCFF002CPK"
for i in lab_gm12878_sonicated_3 lab_gm12878_subsampled_3
do
    for type in plus minus
    do
        mkdir -p ${outdir}/${i}
        f="${indir}/${i}/${motif}.${type}.vsignal.gz"
        out="${outdir}/${i}/${motif}.${type}"
        handle="${motif}.${type}"
        params="--size 2 --xlim 250 --ylim 1.0 --ylim2 500 --alpha 0.1 --split -n ${handle}"
        Rscript ../../bin/makeVplots.R -f ${f} -o ${out} ${params} &
    done
done

echo "Done!"
echo "Results are be at ${outdir}"