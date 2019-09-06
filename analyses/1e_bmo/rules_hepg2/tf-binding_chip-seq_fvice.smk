rule remove_non_overlapping:
    input:
        motif = join(TP_DIR, "{chip_handle}.bed"),
    output:
        join(NOOV_TP_DIR, "{chip_handle}.bed"),
    shell:
        """
        ionice -c2 -n7 ../../bin/filter_bed_co-occurring.py -i {input} \
            -d 500 -c 5 -t strict > {output}
        """

rule vsignal:
    input:
        bam = lambda wildcards: bam_from_sample(wildcards.sample),
        motif = rules.remove_non_overlapping.output
    output:
        vsignal = join(VPLOT_dir, "{sample}", "{chip_handle}.vsignal.gz")
    threads: 8
    params:
        "-r 500 -f 3 -F 4 -F 8 -q 30"
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 ../../bin/measure_signal -p {threads} {params} \
            {input.bam} {input.motif} | gzip -c > {output}
        """

rule vplot:
    input:
        rules.vsignal.output
    output:
        join(VPLOT_dir, "{sample}", "{chip_handle}.png"),
    params:
        main_handle = join(VPLOT_dir, "{sample}", "{chip_handle}"),
        name = "{sample}" + "_{chip_handle}"
    shell:
        """
        Rscript ../../bin/makeVplots.R -f {input} \
            -o {params.main_handle} -n {params.name} --ylim 1.0
        """

rule fetch_fvices:
    input:
        join(VPLOT_dir, "{sample}")
    output:
        prefix_results("output", "{sample}.fvice.chip.out")
    shell:
        """
        ionice -c2 -n7 python ../../bin/get_fvice.py {input} {output}
        """