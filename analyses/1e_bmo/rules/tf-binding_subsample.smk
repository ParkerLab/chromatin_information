# Rule has to be anonymous or it will clash between samples
for sample, depths in subsample_depths.items():
    rule:
        input:
            bam_from_sample(sample)
        output:
            expand(
                join(BAM_SUB_DIR, "{sample}.{depth}M.bam"),
                sample=sample, depth=depths
            ),
            expand(
                join(BAM_SUB_DIR, "{sample}.{depth}M.bam.bai"),
                sample=sample, depth=depths
            )
        resources:
            io_limit = 2
        run:
            for depth in depths:
                if depth < depths[-1]:
                    fraction = subsampled_bam_fraction_from_depth(
                        sample, depth)
                    out = subsampled_bam_names_from_depth(sample, depth)
                    if not os.path.isfile(out):
                        shell(
                            """
                            ionice -c2 -n7 samtools view -b \
                                -s {fraction} {input} > {out} && \
                            samtools index {out}
                            """)
                else:
                    out = subsampled_bam_names_from_depth(sample, depth)
                    ln = os.path.abspath(input[0])
                    shell("ln -sf {ln} {out} && ln -sf {ln}.bai {out}.bai")

rule call_subsampled_peaks:
    input:
        join(BAM_SUB_DIR, "{sample}.{depth}M.bam")
    output:
        join(MACS2_SUB_DIR, "{sample}.{depth}_peaks.broadPeak")
    params:
        name = '{sample}.{depth}',
        genome_size = 'hs',
        stderr_location = join(MACS2_SUB_DIR, "{sample}.{depth}.log"),
        outdir = "work/tf-binding/subsampled_peaks"
    resources:
        io_lmit = 1
    shell:
        """
        ionice -c2 -n7 macs2 callpeak -t {input} --outdir {params.outdir} \
            -f BAM -n {params.name} -g {params.genome_size} --nomodel \
            --shift -100 --seed 762873 --extsize 200 -B --broad \
            --keep-dup all &> {params.stderr_location}
        """

rule blacklist_filter:
    input:
        rules.call_subsampled_peaks.output
    output:
        join(MACS2_SUB_DIR,
             "{sample}.{depth}_peaks.broadPeak.fdr0.05.noblacklist")
    params:
        BL1="../../data/annot/wgEncodeDukeMapabilityRegionsExcludable.bed.gz",
        BL2="../../data/annot/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 intersectBed -v -a {input} -b {params.BL1} \
            {params.BL2} | awk '$9>1.30103{{print}}' > {output}
        """

rule motifs_in_peaks:
    input:
        motif = join(MOTIF_DIR, "{motif}.bed.gz"),
        peaks = rules.blacklist_filter.output,
    output:
        join(PEAKS2_DIR, "{sample}/{motif}.{depth}.bed")
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

rule motifs_in_peaks2:
    input:
        motif = join(MOTIF_DIR, "{motif}.bed.gz"),
        peaks = rules.blacklist_filter.output
    output:
        join(PEAKS2_DIR, "{sample}/{motif}.{depth}.bed")
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 intersectBed -f 1.0 -loj -a {input.motif} \
            -b {input.peaks} | cut -f 1-6,14 | \
            sed 's/\\.$/0/g' > {output}
        """