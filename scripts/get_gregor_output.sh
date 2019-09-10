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