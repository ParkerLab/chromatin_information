## CENTIPEDE
rule make_cut_matrix:
    input:
        motif = os.path.join(MOTIF_DIR, "{motif}.bed.gz"),
        bam = os.path.join(BAM_SUB_DIR, "{sample}.{depth}M.bam"),
    output:
        temp(
            os.path.join(
                CENTI_DIR, "{sample}", "matrices", "{motif}.{depth}.matrix.gz"
            )
        )
    params:
        "-v -d -r 100 -f 3 -F 4 -F 8 -q 30"
    threads: 10
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 make_cut_matrix -p {threads} {params} {input.bam} \
            {input.motif} | gzip -c > {output}
        """

rule centipede:
    input:
        motif = os.path.join(MOTIF_DIR, "{motif}.bed.gz"),
        matrix = rules.make_cut_matrix.output
    output:
        protected(
            os.path.join(
                CENTI_DIR, "{sample}", "posteriors", "{motif}.{depth}.bed.gz"
            )
        )
    params:
        score_column = 5
    shell:
        """
        ionice -c2 -n7 Rscript ~/scripts/run_centipede.R {input.matrix} \
            {input.motif} {output} {params}
        """
