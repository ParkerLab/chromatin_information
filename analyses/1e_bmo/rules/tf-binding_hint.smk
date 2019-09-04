# HINT
rule hint:
    input:
        bam = os.path.join(BAM_SUB_DIR, "{sample}.{depth}M.bam"),
        peaks = rules.blacklist_filter.output
    output:
        os.path.join(
            HINT_DIR,
            "{sample}",
            "output",
            "{sample}_{depth}_HINT.bed"
        ),
    params:
        output = os.path.join(HINT_DIR, "{sample}", "output"),
        handle = "{sample}_{depth}_HINT"
    shell:
        """
        module load rgt/0.12.1
        ionice -c2 -n7 rgt-hint footprinting --paired-end --atac-seq \
            --output-location={params.output} \
            --output-prefix={params.handle} {input.bam} {input.peaks}
        sortBed -i {output} > {output}.tmp && mv {output}.tmp {output}
        """

rule intersect_with_motifs_hint:
    input:
        footprints = rules.hint.output,
        motif = lambda wildcards: "{}/{}.bed.gz".format(
            config['motif_dir'], wildcards.motif)
    output:
        os.path.join(
            HINT_DIR,
            "{sample}",
            "motif_intersections",
            "{motif}.{depth}.bed"
        )
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 intersectBed -wo -a {input.motif} \
            -b {input.footprints} | sortBed | cut -f 1-6,11 | sort -k 1,1 \
            -k2,2n -k7,7nr | sort -u -k 1,1 -k2,2n | intersectBed -wao \
            -f 1.0 -a {input.motif} -b stdin | cut -f 1-6,13 |\
            sed 's/\t\.$/\t0/g' > {output}
        """
