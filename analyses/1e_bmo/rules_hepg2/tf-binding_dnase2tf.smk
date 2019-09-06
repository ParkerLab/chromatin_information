import os

rule make_d2tf_input:
    input:
        lambda wildcards: os.path.abspath(
            subsampled_bam_names_from_depth(wildcards.sample, wildcards.depth)
        )
    output:
        linked = join(D2TF_DIR, "{sample}/input/{sample}.{depth}.bam"),
        ok_file = join(D2TF_DIR, "{sample}/input/{sample}.{depth}/done.ok"),
    params:
        r_script = join(config["dnase2tf_dir"], "dependencies",
                        "bam_compact_split_util",
                        "paired_end_bam2split.r"),
    resources:
        io_limit = 1
    shell:
        """
        ln -sf {input} {output.linked}
        ln -sf {input}.bai {output.linked}.bai
        ionice -c2 -n7 Rscript {params} {output.linked} && \
        touch {output.ok_file}
        """

rule offset_for_atac:
    input:
        done = join(D2TF_DIR, "{sample}/input/{sample}.{depth}/done.ok")
    output:
        join(D2TF_DIR, "{sample}/input/{sample}.{depth}_offset/done.ok")
    params:
        offset = "../../bin/dnase2tf_offset.py",
        indir = join(D2TF_DIR, "{sample}/input/{sample}.{depth}"),
        outdir = join(D2TF_DIR, "{sample}/input/{sample}.{depth}_offset"),
    resources:
        io_limit = 1
    shell:
        """
        for f in `ls {params.indir}/*.txt`; \
        do base=`basename ${{f}}`; \
            out="{params.outdir}/${{base}}";\
            ionice -c2 -n7 python {params.offset} ${{f}} > ${{out}}; \
        done && touch {output}
        """

rule dinucleotide_table:
    input:
        join(BAM_SUB_DIR, "{sample}.{depth}M.bam"),
    output:
        join(D2TF_DIR, "{sample}", "input",
             "dinuc_freq_table_ac_{sample}.{depth}M.txt")
    params:
        calcDFT = "../../bin/calcDFT",
        hg19 = config["hg19_dir"]
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 {params.calcDFT} {params.hg19} {input} &&\
        mv dinuc_freq_table_ac_{wildcards.sample}.{wildcards.depth}M.txt \
            {output}
        """

rule prepare_peaks:
    input:
        rules.blacklist_filter.output
    output:
        os.path.join(D2TF_DIR, "{sample}/input/{sample}.{depth}.peaks.bed")
    resources:
        io_limit = 1
    shell:
        "ionice -c2 -n7 cut -f 1-3 {input} > {output}"

rule DNase2TF:
    input:
        dft = join(D2TF_DIR, "{sample}", "input",
                   "dinuc_freq_table_ac_{sample}.{depth}M.txt"),
        peaks = rules.prepare_peaks.output,
        indir = join(D2TF_DIR, "{sample}", "input",
                     "{sample}.{depth}_offset/done.ok")
    output:
        protected(
            os.path.join(
                D2TF_DIR,
                "{sample}/output/{sample}.{depth}_fdr1.000000.bed"
            )
        )
    params:
        script = "../../bin/run_d2tf.R",
        bam_handle = join(D2TF_DIR, "{sample}",
                          "input/{sample}.{depth}_offset/{sample}.{depth}"),
        out = os.path.join(D2TF_DIR, "{sample}/output/{sample}.{depth}")
    shell:
        """
        ionice -c2 -n7 Rscript {params.script} {params.bam_handle} \
            {params.out} T {input.peaks} {input.dft}
        """

rule fix_d2tf_out:
    input:
        rules.DNase2TF.output
    output:
        os.path.join(D2TF_DIR, "{sample}/output/{sample}.{depth}_fixed.bed")
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 tail -n +2 {input} | sed 's/-\([0-9]*\)/\\1/g; \
            s/Inf/308.0000/g' | sortBed  > {output}
        """

rule intersect_with_motifs_d2tf:
    input:
        footprints = rules.fix_d2tf_out.output,
        motif = lambda wildcards: "{}/{}.bed.gz".format(
            config['motif_dir'], wildcards.motif)
    output:
        os.path.join(
            D2TF_DIR, "{sample}/motif_intersections/{motif}.{depth}.bed"
        )
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 intersectBed -wo -a {input.motif} \
            -b {input.footprints} | sortBed | cut -f 1-6,10 | sort -k 1,1 \
            -k2,2n -k7,7nr | sort -u -k 1,1 -k2,2n | intersectBed -wao \
            -f 1.0 -a {input.motif} -b stdin | cut -f 1-6,13 | \
            sed 's/\t\.$/\t0/g' > {output}
        """
