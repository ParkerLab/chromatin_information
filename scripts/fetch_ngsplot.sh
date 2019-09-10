indir="/lab/work/albanus/2018_gm12878_mnase-seq/work/ngsplot_per_motif/\
lab_gm12878"
outdir="data/processed/mnase/lab_gm12878"
mkdir -p ${outdir}
for f in `ls ${indir}/SRR452483__*/avgprof.RData`
do
    base=`echo ${f} | sed 's/.*\/\(SRR\)/\1/g; s/\//./g'`
    cp ${f} ${outdir}/${base}
done

