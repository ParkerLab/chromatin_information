import os

rule subsampled_bam_to_bed:
    input:
        os.path.join(BAM_SUB_DIR, "{sample}.{depth}M.bam"),
    output:
        os.path.join(BED_SUB_DIR, "{sample}.{depth}M.bed")
    shell:
        """
        ionice -c2 -n7 bamToBed -i {input} > {output}
        """


rule call_subsampled_peaks:
    input:
        rules.subsampled_bam_to_bed.output,
    output:
        os.path.join(MACS2_SUB_DIR, "{sample}.{depth}_peaks.broadPeak")
    params:
        name = '{sample}.{depth}',
        genome_size = 'hs',
        stderr_location = os.path.join(MACS2_SUB_DIR, "{sample}.{depth}.log"),
        outdir = MACS2_SUB_DIR
    resources:
        io_lmit = 1
    shell:
        """
        ionice -c2 -n7 macs2 callpeak -t {input} --outdir {params.outdir} \
            -f BED -n {params.name} -g {params.genome_size} --nomodel \
            --shift -100 --seed 762873 --extsize 200 -B --broad \
            --keep-dup all &> {params.stderr_location}
        """

rule blacklist_filter:
    input:
        rules.call_subsampled_peaks.output
    output:
        os.path.join(
            MACS2_SUB_DIR,
            "{sample}.{depth}_peaks.broadPeak.fdr0.05.noblacklist"
        )
    params:
        BL1 = os.path.join(
            "/lab/data/reference/human/hg19/annot",
            "wgEncodeDukeMapabilityRegionsExcludable.bed.gz"
        ),
        BL2 = os.path.join(
            "/lab/data/reference/human/hg19/annot",
            "wgEncodeDacMapabilityConsensusExcludable.bed.gz"
        ),
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 intersectBed -v -a {input} -b {params.BL1} \
            {params.BL2} | awk '$9>1.30103{{print}}' > {output}
        """

rule motifs_in_peaks:
    input:
        motif = lambda wildcards: "{}/{}.bed.gz".format(
            config['motif_dir'], wildcards.motif),
        peaks = rules.blacklist_filter.output
    output:
        os.path.join(PEAKS_DIR, "{sample}/{motif}.{depth}.bed")
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 intersectBed -f 1.0 -loj -a {input.motif} \
            -b {input.peaks} | cut -f 1-7 | \
            awk '{{OFS="\\t"; if($7 == ".")\
            {{print $1,$2,$3,$4,$5,$6,0}}else{{print $1,$2,$3,$4,$5,$6,1}}}}'\
            > {output}
        """
