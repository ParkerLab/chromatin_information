# BMO
RS_DIR = prefix_results("raw_signals", "{sample}")
NB_DIR = prefix_results("negative_binomials", "{sample}")

rule measure_raw_signal:
    input:
        motif = os.path.join(MOTIF_DIR, "{motif}.bed.gz"),
        bam = os.path.join(BAM_SUB_DIR, "{sample}.{depth}M.bam"),
    output:
        os.path.join(RS_DIR, "all_regions", "{motif}.{depth}.bed"),
    threads: 10
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 ../../bin/measureRawSignal.py -p {threads} \
            -b {input.bam} -m {input.motif} > {output}
        """

rule motifs_outside_peaks:
    input:
        motif = rules.measure_raw_signal.output,
        peaks = rules.blacklist_filter.output
    output:
        os.path.join(RS_DIR, "outside_peaks", "{motif}.{depth}.bed")
    shell:
        """
        ionice -c2 -n7 intersectBed -v -a {input.motif} \
            -b {input.peaks} > {output}
        """

rule count_overlapping_motifs:
    input:
        rules.measure_raw_signal.output
    output:
        os.path.join(RS_DIR, "co_occurring", "{motif}.{depth}.bed")
    params: 100
    shell:
        """
        ionice -c2 -n7 ../../bin/count_co-occuring_motifs.sh {input} \
            {params} > {output}
        """

rule fit_nbinoms:
    input:
        raw_signal = rules.measure_raw_signal.output,
        outside_peaks = rules.motifs_outside_peaks.output
    output:
       os.path.join(NB_DIR, "{motif}.{depth}.bed.gz")
    params:
        in_handle = "{motif}.{depth}.bed",
        d1 = os.path.join(RS_DIR, "all_regions/"),
        d2 = os.path.join(RS_DIR, "outside_peaks/"),
        out_handle = os.path.join(NB_DIR, "{motif}.{depth}")
    shell:
        """
        ionice -c2 -n7 Rscript ../../bin/rawSigNBModel.R \
            -f {params.in_handle} --dir1 {params.d1} --dir2 {params.d2} \
            -o {params.out_handle} -c 7 --writeBed
        """

rule BMO:
    input:
        atac_nb = rules.fit_nbinoms.output,
        motif_counts = os.path.join(
            RS_DIR, "co_occurring", "{motif}.{depth}.bed"
        ),
    output:
        os.path.join(BMO_DIR, "{sample}", "bound", "{motif}.{depth}.bound.bed")
    params:
        bmo_output_dir = prefix_results("bmo", "{sample}")
    shell:
        """
        ionice -c2 -n7 Rscript ../../bin/bmo.R --f1 {input.atac_nb} \
            --f2 {input.motif_counts} -o {params.bmo_output_dir} \
            -n {wildcards.motif}.{wildcards.depth} -p 1
        """
