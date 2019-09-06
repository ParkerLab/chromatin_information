def get_max_depth(sample):
    depths = subsample_depths[sample]
    return max(depths)

rule remove_non_overlapping_bmo:
    input:
        motif = lambda wildcards: join(
            BMO_DIR, "{sample}", "bound", 
            "{motif}" + ".{}.bound.bed".format(get_max_depth(wildcards.sample))
        ),
    output:
        join(NOOV_DIR, "{motif}.bed"),
    shell:
        """
        ionice -c2 -n7 ../../bin/filter_bed_co-occurring.py -i {input} \
            -d 500 -c 5 -t strict > {output}
        """

rule vsignal_bmo:
    input:
        bam = lambda wildcards: bam_from_sample(wildcards.sample),
        motif = rules.remove_non_overlapping_bmo.output
    output:
        vsignal = join(VPLOT_BMO_dir, "{sample}", "{motif}.vsignal.gz")
    threads:
        8
    params:
        "-r 500 -f 3 -F 4 -F 8 -q 30"
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 ../../bin/measure_signal -p {threads} {params} \
            {input.bam} {input.motif} | gzip -c > {output}
        """

rule vplot_bmo:
    input:
        rules.vsignal_bmo.output
    output:
        join(VPLOT_BMO_dir, "{sample}", "{motif}.png"),
    params:
        main_handle = join(VPLOT_BMO_dir, "{sample}", "{motif}"),
        name = "{sample}" + "_{motif}"
    shell:
        """
        Rscript ../../bin/makeVplots.R -f {input} \
            -o {params.main_handle} -n {params.name} --ylim 1.0
        """

rule fetch_fvices_bmo:
    input:
        join(VPLOT_BMO_dir, "{sample}")
    output:
        prefix_results("output", "{sample}.fvice.bmo.out")
    shell:
        """
        ionice -c2 -n7 python bin/get_fvice.py {input} {output}
        """
