from os.path import join
import math


associations = "associations_hepg2_top1.txt"

MOTIF_DIR = config["motif_dir"]
PREV = join("work", "tf-binding")
BMO_DIR = join(PREV, 'bmo', "{sample}", "all")
PIQ_DIR = join(PREV, 'piq')
TP_dir = "work/tf-binding/true_positives"

# Target directories (relative to root)
results = join("work", "piq-centered")
PIQ_DIR2 = join(results, "piq", "{sample}")
BMO_DIR2 = join(results, 'bmo', "{sample}")
F1_dir = join(results, "f1", "{sample}")
ROC_dir = join(results, "roc", "{sample}")


def get_motifs():
    """List all motifs"""
    with open(config["chip_assoc"], "r") as stream:
        header = stream.readline()
        motifs = []
        for line in stream:
            line = line.strip().split("\t")
            motif = line[7]
            motifs.append(motif)
        motifs = list(set(motifs))
        motifs.sort()
        return motifs


def get_samples():
    """List all samples"""
    for sample in sorted(config['samples'].keys()):
        yield sample


def bam_from_sample(sample):
    """Return the bam file associated with a sample"""
    bamfile = config['samples'][sample]['bamfile']
    return bamfile


def subsampled_bam_from_sample(sample):
    """Return subsampling depths and fractions for a sample"""
    depth = config['samples'][sample]['depth'] / 10**6
    sub_depths = [5, 10]
    while sub_depths[-1] + 10 < depth:
        sub_depths.append(sub_depths[-1] + 10)
    sub_depths.append(math.floor(depth))
    bam_names = [os.path.join(
        BAM_SUB_DIR, "{}.{}M.bam"
    ).format(sample, x) for x in sub_depths]
    return bam_names


def subsampled_depths_from_sample(sample):
    """Return subsampling depths and fractions for a sample"""
    depth = config['samples'][sample]['depth'] / 10**6
    sub_depths = [5, 10]
    while sub_depths[-1] + 10 < depth:
        sub_depths.append(sub_depths[-1] + 10)
    sub_depths.append(math.floor(depth))
    # return sub_depths
    return [max(sub_depths)]  # only highest depth


def subsampled_bam_names_from_depth(sample, sub_depth):
    """Return subsampling depths and fractions for a sample"""
    bam_name = os.path.join(
        BAM_SUB_DIR, "{}.{}M.bam".format(sample, int(sub_depth))
    )
    return bam_name


def subsampled_bam_fraction_from_depth(sample, sub_depth):
    """Return subsampling depths and fractions for a sample"""
    depth = config['samples'][sample]['depth'] / 10**6
    fraction = sub_depth / depth
    return fraction


def get_piq_handle(motif):
    handle = motif_handle_dic[motif]
    return handle


def get_motif_index(motif):
    index = motif_dic[motif]
    return index


def get_motif_from_handle(handle):
    """Split ChIP-seq handle experiment__tf__motif"""
    handle = handle.split("__")
    return handle[2]


def get_dir(method):
    """Method's output directory"""
    return method


def get_column(method):
    """Column for F1 calculation"""
    return "7"


def get_roc_column(method):
    """Column for ROC-PR calculation"""
    if method == "piq":
        return "5"
    else:
        return "7"


def get_threshold(method):
    """What is the method's cutoff threshold and direction (lt, gt)"""
    if method == "bmo":
        return "1.30103"
    else:
        return "0"


# Subsample dictionary - will be used to determine the output files for
# each sample
subsample_depths = {
    x: subsampled_depths_from_sample(x) for x in get_samples()
}
# subsample_depths = {"buenrostro_rep1" : [80]} # debug

# Number of motifs and motif-number dictionary - used for PIQ
# e.g. 9-CTCFknown2, 20-NFKBknown10
all_motifs = [x for x in get_motifs()]
all_motifs.sort()
nmotifs = len(all_motifs)
motif_dic = {}
motif_handle_dic = {}
for i in range(0, nmotifs):
    motif = all_motifs[i]
    handle = "{}-{}".format(i + 1, all_motifs[i].replace("_", ""))
    motif_dic[motif] = str(i + 1)
    motif_handle_dic[motif] = handle

# Get al the motif - ChIP associations
chip_handles = []
with open(associations, "r") as stream:
    header = stream.readline()
    for line in stream:
        line = line.strip().split("\t")
        handle = line[8]
        chip_handles.append(handle)

methods = ["bmo", "piq"]


rule all:
    input:
        [expand(join(PIQ_DIR2, "{motif}.{depth}.bed"),
                sample=sample, depth=depths, motif=get_motifs()
         ) for sample, depths in subsample_depths.items()],
        [expand(join(BMO_DIR2, "{motif}.{depth}.bed"),
                sample=sample, depth=depths, motif=get_motifs()
         ) for sample, depths in subsample_depths.items()],
        

rule all_metrics:
    input:
        [expand(join(F1_dir, "{sample}.{method}.{handle}__{depth}.out"),
               sample=sample, depth=depths, method=methods, 
               handle=chip_handles
         )for sample, depths in subsample_depths.items()],
        [expand(join(ROC_dir, "{sample}.{method}.{handle}__{depth}.prauc.out"),
               sample=sample, depth=depths, method=methods, 
               handle=chip_handles
         )for sample, depths in subsample_depths.items()],


rule aggregate_f1:
    input:
        [expand(join(F1_dir, "{sample}.{method}.{handle}__{depth}.out"),
               sample=sample, depth=depths, method=methods, 
               handle=chip_handles
         )for sample, depths in subsample_depths.items()],
    output:
        join(results, "f1_data.txt")
    shell:
        """
        echo -e "precision\\trecall\\tf1\\tmethod\\texp\\ttrue_positives" \
            > {output}
        cat {input} >> {output}
        """

rule aggregate_prauc:
    input:
        [expand(join(ROC_dir, "{sample}.{method}.{handle}__{depth}.prauc.out"),
               sample=sample, depth=depths, method=methods, 
               handle=chip_handles
         )for sample, depths in subsample_depths.items()],
    output:
        join(results, "prauc_data.txt")
    shell:
        """
        echo -e "motif\\tmethod\\tprauc\\ttrue_positives" > {output}
        cat {input} >> {output}
        """


rule piq:
    input:
        piq = join(PIQ_DIR, "processed_output_scores", "{sample}", 
                   "{depth}M.{motif}.bed"),
    output:
        binary = join(PIQ_DIR2, "{motif}.{depth}.bed"),
    params:
        handle = lambda wildcards: get_piq_handle(wildcards.motif),
        handle_tmp = "binary_{sample}__{depth}__{motif}",
        f_fw = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out",
                "{}-calls.csv".format(get_piq_handle(wildcards.motif))
            ),
        f_rc = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out",
                "{}.RC-calls.csv".format(get_piq_handle(wildcards.motif))
            ),
        f_size_fw = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out",
                "{}-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
        f_size_rc = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out",
                "{}.RC-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
    shell:
        """
        size_fw=`head -n 25 {params.f_size_fw} | awk '{{print $3 - $2}}' | \
            sort -nr | head -1`
        size_rc=`head -n 25 {params.f_size_rc} | awk '{{print $3 - $2}}' | \
            sort -nr | head -1`
        tmp_fw="{params.handle_tmp}.fw.bed.tmp"
        tmp_rc="{params.handle_tmp}.rc.bed.tmp"
        tmp_scores="{params.handle_tmp}.scores.bed.tmp"

        ionice -c2 -n7 intersectBed -u -wa -a {input.piq} \
            -b {params.f_size_fw} {params.f_size_rc} | 
            grep -vP \"\\t-1\\t\" > ${{tmp_scores}} && \

        tail -n+2 {params.f_fw} | sed 's/\"//g; s/,/\\t/g' | \
            awk \"{{OFS=\\\"\\t\\\"; print \\$2,\\$3,\\$3 + ${{size_fw}}}}\" \
            > ${{tmp_fw}} && \
        tail -n+2 {params.f_rc} | sed 's/\"//g; s/,/\\t/g' | \
            awk \"{{OFS=\\\"\\t\\\"; print \\$2,\\$3,\\$3 + ${{size_rc}}}}\" \
            > ${{tmp_rc}} && \
        ionice -c2 -n7 cut -f 1-6 ${{tmp_scores}} | \
            intersectBed -c -a stdin -b ${{tmp_fw}} ${{tmp_rc}} | \
            sed 's/\\t[1-9]$/\\t1/g' > {output} && \
        rm ${{tmp_fw}} ${{tmp_rc}} ${{tmp_scores}}
        """

rule get_bmo:
    input:
        bmo = join(BMO_DIR, "{motif}.{depth}.all.bed"),
        piq = rules.piq.output,
    output:
        join(BMO_DIR2, "{motif}.{depth}.bed"),
    shell:
        """
        ionice -c2 -n7 intersectBed -f 1.0 -F 1.0 -a {input.bmo} \
            -b {input.piq} | cut -f 1-6,13 > {output}
        """

rule get_f1:
    input:
        motif = lambda wildcards: join(
            results, "{}".format(get_dir(wildcards.method)),
            "{}".format(wildcards.sample),
            "{}.{}.bed".format(
                get_motif_from_handle(wildcards.handle),
                wildcards.depth
            )
        ),
        tp = join(TP_dir, "{handle}.bed")
    output:
        join(F1_dir, "{sample}.{method}.{handle}__{depth}.out"),
    params:
        main = "Rscript bin/calculate_F1_individual.R",
        outdir = F1_dir,
        col = lambda wildcards: get_column(wildcards.method),
        handle2 = "{sample}.{method}.{handle}__{depth}",
        method = "{method}",
        thresh = lambda wildcards: get_threshold(wildcards.method),
    shell:
        """
        ionice -c2 -n7 {params.main} -f {input.motif} -o {params.outdir} \
            -c {params.col} -n {params.handle2} -m {params.method} \
            -t {input.tp} -r {params.thresh}
        """

rule get_roc_prauc:
    input:
        motif = lambda wildcards: join(
            results, "{}".format(get_dir(wildcards.method)),
            "{}".format(wildcards.sample),
            "{}.{}.bed".format(
                get_motif_from_handle(wildcards.handle),
                wildcards.depth
            )
        ),
        tp = join(TP_dir, "{handle}.bed")
    output:
        join(ROC_dir, "{sample}.{method}.{handle}__{depth}.prauc.out"),
    params:
        main = "Rscript ~albanus/scripts/get_ROC_single.R",
        outdir = ROC_dir,
        col = lambda wildcards: get_roc_column(wildcards.method),
        handle2 = "{sample}.{method}.{handle}__{depth}",
        method = "{method}",
    shell:
        """
        ionice -c2 -n7 {params.main} -f {input.motif} -c {params.col} \
            -n {params.handle2} -t {input.tp} -m {params.method} \
            -o {params.outdir}
        """