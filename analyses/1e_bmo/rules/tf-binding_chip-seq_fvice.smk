rule vsignal:
    input:
        bam = "work/mapping/prune/{sample}.pruned.bam",
        motif = "/lab/work/albanus/gm12878/motif_intersects_encode2017/intersections_no_co-occurring/{handle}.bed"
    output:
        vsignal = "work/tf-binding/vplots_chip-seq/{sample}/{handle}.vsignal.gz"
    threads: 8
    params:
        "-r 500 -f 3 -F 4 -F 8 -q 30"
    resources:
        io_limit = 1
    shell:
        "ionice -c2 -n7 measure_signal -p {threads} {params} {input.bam} {input.motif} | gzip -c > {output}"

rule vplot:
    input:
        rules.vsignal.output
    output:
        # "work/tf-binding/vplots_chip-seq/{sample}/{handle}.png"
        "../2019_vplots_new/work/tf-binding/vplots_chip-seq/{sample}/{handle}.png"
    params:
        # main_handle = "work/tf-binding/vplots_chip-seq/{sample}/{handle}",
        main_handle = "../2019_vplots_new/work/tf-binding/vplots_chip-seq/{sample}/{handle}",
        name = "{sample}" + "_{handle}"
    shell:
        "Rscript ~albanus/scripts/makeVplots.R -f {input} -o {params.main_handle} -n {params.name} --ylim 1.0"

rule fetch_fvices:
    input:
        # "work/tf-binding/vplots_chip-seq/{sample}"
        "../2019_vplots_new/work/tf-binding/vplots_chip-seq/{sample}"
    output:
        # "work/tf-binding/output/{sample}.fvice.out"
        "work/tf-binding/output/{sample}.new_fvice.out"
    shell:
        "ionice -c2 -n7 python bin/f-vices/get_fvice.py {input} {output}"