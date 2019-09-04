from os.path import join


def get_samples():
    for sample in config["samples"]:
        yield sample


def get_motifs():
    for motif in config["motifs"]:
        yield motif


RESULTS = config["results"]
VSIGNAL_dir = join(config["vsignal_dir"], "{sample}")
VPLOT_dir = join(RESULTS, "vplots", "{sample}")
VPLOT2_dir = join(RESULTS, "vplots_same_n_motif", "{sample}")


rule all:
    input:
        expand(join(VPLOT_dir, "{motif}.png"), 
               sample=get_samples(), motif=get_motifs()),
        expand(join(VPLOT2_dir, "{motif}.png"), 
               sample=get_samples(), motif=get_motifs())


rule filter_vsignals:
    input:
        join(VSIGNAL_dir, "{motif}.vsignal.gz")
    output:
        clean = join(VPLOT_dir, "{motif}.vsignal.filtered.gz"),
        size = join(VPLOT_dir, "{motif}.size"),
        n = join(VPLOT_dir, "{motif}.nmotifs"),
    shell:
        """
        zcat {input} | awk '$4 > 40 {{print}}' | gzip -c > {output.clean} && \
        zcat {output.clean} | tail -n +2 | wc -l > {output.size} && \
        zcat {output.clean} | tail -n +2 | cut -f 6 | sort -u | \
            wc -l > {output.n}
        """

rule get_smallest:
    input:
        lambda wildcards: expand(
            join(VPLOT_dir, "{motif}.size"), 
            motif=get_motifs(), sample=wildcards.sample
            )
    output:
        join(RESULTS, "min_size.{sample}.txt")
    run:
        sizes = []
        for f in input:
            with open(f, "r") as n:
                n = n.readline().strip()
                sizes.append(int(n))
        min_size = min(sizes)
        shell("echo {min_size} > {output}")
        print(min_size)

rule get_smallest2:
    input:
        lambda wildcards: expand(
            join(VPLOT_dir, "{motif}.nmotifs"), 
            motif=get_motifs(), sample=wildcards.sample
            )
    output:
        join(RESULTS, "min_motifs.{sample}.txt")
    run:
        nmotifs = []
        for f in input:
            print(f)
            with open(f, "r") as n:
                n = n.readline().strip()
                nmotifs.append(int(n))
        min_motif = min(nmotifs)
        shell("echo {min_motif} > {output}")
        print(min_motif)

rule downsample:
    input:
        frags = rules.get_smallest.output,
        motifs = rules.get_smallest2.output,
        vsignal = rules.filter_vsignals.output.clean,
    output:
        join(VPLOT_dir, "{motif}.downsampled.vsignal.gz")
    shell:
        """
        Rscript scripts/downsample.R {input.vsignal} {input.frags} \
            {input.motifs} {output}
        """
    
rule vplots:
    input:
        rules.downsample.output,
    output:
        join(VPLOT_dir, "{motif}.png")
    params:
        out = join(VPLOT_dir, "{motif}"),
        motif = "{motif}"
    shell:
        """
        Rscript ../../bin/makeVplots.R -f {input} \
            -o {params.out} -n {params.motif} --ylim 1.25 --split
    """

rule downsample2:
    input:
        motifs = rules.get_smallest2.output,
        vsignal = rules.filter_vsignals.output.clean,
    output:
        join(VPLOT2_dir, "{motif}.downsampled.vsignal.gz")
    shell:
        """
        Rscript scripts/downsample2.R {input.vsignal} {input.motifs} {output}
        """
    
rule vplots2:
    input:
        rules.downsample2.output,
    output:
        join(VPLOT2_dir, "{motif}.png")
    params:
        out = join(VPLOT2_dir, "{motif}"),
        motif = "{motif}"
    shell:
        """
        Rscript ~albanus/scripts/makeVplots.R -f {input} \
            -o {params.out} -n {params.motif} --ylim 1.25 --split \
            --maxfrags 10000000 --maxpoints 10000000
    """