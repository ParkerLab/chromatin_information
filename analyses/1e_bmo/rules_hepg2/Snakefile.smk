import itertools
import os
import sys
import functools
import math

from os.path import join

###
# GENERAL
#

# Data

MOTIF_DIR = config["motif_dir"]
CHIP_dir = config["chip_dir"]

# Path generator function
prefix_results = functools.partial(join, config['results'])

# Target directories (relative to root)
BAM_SUB_DIR = prefix_results("bam_subsampled")
MACS2_SUB_DIR = prefix_results("macs2_subsampled")
TP_DIR = prefix_results('true_positives')
BMO_DIR = prefix_results('bmo')
CENTI_DIR = prefix_results('centipede')
CENTI_SB_DIR = prefix_results('centipede_singlebin')
D2TF_DIR = prefix_results('dnase2tf')
HINT_DIR = prefix_results('hint')
PIQ_DIR = prefix_results('piq')
MS_CENTI_DIR = prefix_results("mscentipede", "{sample}")
PEAKS_DIR = prefix_results('motifs_in_atac-seq_peaks')
PEAKS2_DIR = prefix_results('motifs_in_atac-seq_peaks_scores')
CONCAT_DIR = prefix_results('concatenated_files')
COMPARE_DIR = prefix_results('comparisons')
TP_DIR = prefix_results("true_positives")
NOOV_TP_DIR = prefix_results("true_positives_noov")
NOOV_DIR = prefix_results("bmo_bound_noov", "{sample}")
VPLOT_dir = prefix_results("vplots_chip-seq")
VPLOT_BMO_dir = prefix_results("vplots_bmo")

# Logs
LOG_DIR = prefix_results('logs')
VERSION_DIR = prefix_results('versions')


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


def get_handles():
    """List all motif-ChIP handles"""
    with open(config["chip_assoc"], "r") as stream:
        header = stream.readline()
        handles = []
        for line in stream:
            line = line.strip().split("\t")
            handle = line[8]
            handles.append(handle)
        motifs = list(set(handles))
        return handles


def get_exp_from_handle(handle):
    return handle.split("__")[0]


def get_motif_from_handle(handle):
    return handle.split("__")[2]


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
    return sub_depths
    # return [max(sub_depths)]  # only highest depth


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


def get_chip_handles():
    with open(config["chip_assoc"]) as stream:
        header = stream.readline()
        for line in stream:
            line = line.strip().split("\t")
            handle = line[8]
            yield(handle)


def get_motif_from_handle(handle):
    handle = handle.split("__")
    motif = handle[2]
    f = join(config["motif_dir"], "{}.bed.gz".format(motif))
    return f


def get_chip_from_handle(handle):
    handle = handle.split("__")
    chip = handle[0]
    f = join(config["chip_dir"], "{}.bed.gz".format(chip))
    return f


# Subsample dictionary - will be used to determine the output files for
# each sample
subsample_depths = {x: subsampled_depths_from_sample(x) for x in get_samples()}
# subsample_depths = {"hepg2_3" : [34]} # debug

# Number of motifs and motif-number dictionary - used for PIQ
# e.g. 9-CTCFknown2, 20-NFKBknown10
motifs = get_motifs()
motifs = list(set(motifs)) # remove duplicates
motifs.sort()
nmotifs = len(motifs)
motif_dic = {}
motif_handle_dic = {}
for i in range(0, nmotifs):
    motif = motifs[i]
    handle = "{}-{}".format(i + 1, motifs[i].replace("_", "").replace(".", ""))
    motif_dic[motif] = str(i + 1)
    motif_handle_dic[motif] = handle



rule all:
    input:
        # Peak calls
        [expand(
            os.path.join(
                MACS2_SUB_DIR,
                "{sample}.{depth}_peaks.broadPeak.fdr0.05.noblacklist"
            ),
            sample=sample,
            depth=depths
        ) for sample, depths in subsample_depths.items()],
        # HINT
        [expand(
            os.path.join(
                HINT_DIR,
                "{sample}/motif_intersections/{motif}.{depth}.bed"
            ),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        # DNase2TF
        [expand(
            os.path.join(
                D2TF_DIR,
                "{sample}/motif_intersections/{motif}.{depth}.bed"
            ),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        # CENTIPEDE
        [expand(
            os.path.join(
                CENTI_DIR, "{sample}/posteriors/{motif}.{depth}.bed.gz"
            ),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        # BMO
        [expand(
            os.path.join(
                BMO_DIR, "{sample}", "bound", "{motif}.{depth}.bound.bed"
            ),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        [expand(
            os.path.join(
                CENTI_SB_DIR, "{sample}/posteriors/{motif}.{depth}.bed.gz"
            ),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        # PIQ - binarized files
        [expand(
            os.path.join(
                PIQ_DIR,
                "processed_output_binary", "{sample}", "{depth}M.{motif}.bed"
            ),
            sample=sample,
            depth=depths,
            motif=motifs
        ) for sample, depths in subsample_depths.items()],
        # PIQ - scores files
        [expand(
            os.path.join(
                PIQ_DIR,
                "processed_output_scores", "{sample}", "{depth}M.{motif}.bed"
            ),
            sample=sample,
            depth=depths,
            motif=motifs
        ) for sample, depths in subsample_depths.items()],
        # Motifs in peaks
        [expand(
            os.path.join(PEAKS_DIR, "{sample}", "{motif}.{depth}.bed"),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        [expand(
            os.path.join(PEAKS2_DIR, "{sample}", "{motif}.{depth}.bed"),
            sample=sample,
            depth=depths,
            motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        # True positives
        expand(
            prefix_results("true_positives", "{chip_handle}.bed"), 
            chip_handle=get_chip_handles()
        ),
        

rule all_evaluate:
    input:
        [expand(
            os.path.join(CONCAT_DIR, "{sample}/{motif}.{depth}.bed"),
            sample=sample,
            depth=depths, motif=get_motifs()
        ) for sample, depths in subsample_depths.items()],
        expand(
            prefix_results("f1", "{sample}", "get_f1.ok"), 
            sample=get_samples()
        ),
        expand(
            prefix_results("roc", "{sample}", "get_roc.ok"), 
            sample=get_samples()
        )

rule all_roc_f1:
    input:
        expand(
            prefix_results("roc", "{sample}", "prauc_data.txt"), 
            sample=get_samples()
        ),
        expand(
            prefix_results("roc", "{sample}", "roc_data.txt"), 
            sample=get_samples()
        ),
        expand(
            prefix_results("f1", "{sample}", "f1_data.txt"), 
            sample=get_samples()
        )

rule all_vplots:
    input:
        expand(
            join(VPLOT_dir, "{sample}", "{chip_handle}.png"), 
            sample=get_samples(), chip_handle=get_chip_handles()
        ),
        expand(
            join(VPLOT_BMO_dir, "{sample}", "{motif}.png"), 
            sample=get_samples(), motif=get_motifs()
        ),


rule all_fvices:
    input:
        expand(
            prefix_results("output", "{sample}.fvice.chip.out"), 
            sample=get_samples(),
        ),
        expand(
            prefix_results("output", "{sample}.fvice.bmo.out"), 
            sample=get_samples(),
        ),

# Call peaks on downsampled data
include: "tf-binding_subsample.smk"

# Run TF-binding prediction methods
include: "tf-binding_hint.smk"
include: "tf-binding_dnase2tf.smk"
include: "tf-binding_centipede.smk"
include: "tf-binding_bmo.smk"
include: "tf-binding_piq.smk"
include: "tf-binding_mscentipede.smk"
include: "tf-binding_centipede_singlebin.smk"

# Compare methods
include: "tf-binding_truepositives.smk"
include: "tf-binding_evaluate.smk"

# Generate V-plots
include: "tf-binding_chip-seq_fvice.smk"
include: "tf-binding_bmo_fvice.smk"