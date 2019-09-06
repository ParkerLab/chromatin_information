rule make_tp:
    input:
        motif = lambda wildcards: get_motif_from_handle(wildcards.handle),
        chip = lambda wildcards: get_chip_from_handle(wildcards.handle),
    output:
        join(TP_DIR, "{handle}.bed")
    shell:
        """
        ionice -c2 -n7 zcat {input.motif} | \
            intersectBed -f 1.0 -a stdin -b {input.chip} > {output}
        """
