# PIQ
MATX_DIR = "../../data/pwm"
PIQ_SW = "../../sw/PIQ"
BIN = "../../bin"

rule fetch_scripts:
    input:
        os.path.join(PIQ_SW, "common.r")
    output:
        os.path.join(PIQ_DIR, "common.r")
    shell:
        """
        cp {PIQ_SW}/*.r {PIQ_DIR}
        """

rule convert_pwms:
    input:
        mat = os.path.join(MATX_DIR, "{motif}.mat"),
        script = rules.fetch_scripts.output
    output:
        os.path.join(PIQ_DIR, "motifs", "{motif}.jaspar")
    params:
        script = os.path.join(BIN, "mat2jaspar.R"),
        motif = "{motif}"
    shell:
        """
        Rscript {params.script} {input.mat} {params.motif} > {output}
        """

rule gather_motifs:
    input:
        expand(
            os.path.join(PIQ_DIR, "motifs", "{motif}.jaspar"),
            motif=get_motifs()
        )
    output:
        mat = os.path.join(PIQ_DIR, "motifs", "matrices.txt")
    params:
        outdir = os.path.join(PIQ_DIR, "motifs")
    shell:
        """
        cat {input} > {output.mat}
        """

rule pwm_key:
    input:
        rules.gather_motifs.output,
    output:
        os.path.join(PIQ_DIR, "motifs", "key.txt")
    shell:
        """
        grep "^>" {input} | sed 's/^>//g' > tmp.keys.txt
        sed 's/_//g' tmp.keys.txt > tmp.keys2.txt
        paste tmp.keys.txt tmp.keys2.txt | cat -n > {output} && \
        rm tmp.keys*txt
        """

rule pwm_scan:
    input:
        mat = os.path.join(PIQ_DIR, "motifs", "matrices.txt"),
    output:
        os.path.join(PIQ_DIR, "pwms", "{i}.pwmout.RData")
    params:
        workdir = PIQ_DIR,
        outdir = "pwms",
        i = "{i}",
    shell:
        """
        basedir=`pwd`
        cd {params.workdir} && \
        Rscript pwmmatch.exact.r common.r {input.mat} {params.i} \
            {params.outdir} && \
        cd ${{basedir}}
        """

rule process_bam_file:
    input:
        os.path.join(BAM_SUB_DIR, "{sample}.{depth}M.bam")
    output:
        os.path.join(PIQ_DIR, "bam", "{sample}.{depth}M.RData")
    params:
        out = os.path.join("bam", "{sample}.{depth}M.RData")
    resources:
        io_limit = 1
    shell:
        """
        basedir=`pwd`
        cd {PIQ_DIR} && \
        ionice -c2 -n7 Rscript pairedbam2rdata_FIXED.r common.r {params} \
            {input} && \
        cd ${{basedir}}
        """

rule piq:
    input:
        bam = rules.process_bam_file.output,
        motif = lambda wildcards:
            os.path.join(
                PIQ_DIR, "pwms", 
                "{}.pwmout.RData".format(get_motif_index(wildcards.motif))
            ),
        all_motifs = expand(
            os.path.join(PIQ_DIR, "pwms", "{i}.pwmout.RData"),
            i=[x for x in range(1, nmotifs + 1)]
        ),
    output:
        os.path.join(
            PIQ_DIR, "out", "{sample}", "{depth}M", "out", "{motif}-piq.OK"
        )
    params:
        outdir = os.path.join(PIQ_DIR, "out", "{sample}", "{depth}M"),
        i = lambda wildcards: get_motif_index(wildcards.motif),
        outfile = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}M".format(wildcards.depth), "out", 
                "{}-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
    shell:
        """
        basedir=`pwd`
        cd {PIQ_DIR} && \
        ionice -c2 -n7 Rscript pertf.r common.r pwms/ {params.outdir}/tmp/ \
            {params.outdir}/out/ {input.bam} {params.i} && \
        cd ${{basedir}} && \
        touch {output}
        """

# process_piq1 is taking the call csv files from PIQ (which have the predicted
# bound instances) and intersecting with the original motif file to assign 
# them a 0/1 score to whether there's a footprint there.
rule process_piq1:
    input:
        motif = os.path.join(MOTIF_DIR, "{motif}.bed.gz"),
        ok = rules.piq.output,
    output:
        os.path.join(
            PIQ_DIR, "processed_output_binary", "{sample}", 
            "{depth}M.{motif}.bed"
        ),
    params:
        handle_tmp = "binary_{sample}__{depth}__{motif}",
        f_fw = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}-calls.csv".format(get_piq_handle(wildcards.motif))
            ),
        f_rc = lambda wildcards:
            join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}.RC-calls.csv".format(get_piq_handle(wildcards.motif))
            ),
        f_size_fw = lambda wildcards:
            join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
        f_size_rc = lambda wildcards:
            join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}.RC-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
    shell:
        """
        size_fw=`head -n 25 {params.f_size_fw} | awk '{{print $3 - $2}}' | \
            sort -nr | head -1`
        size_rc=`head -n 25 {params.f_size_rc} | awk '{{print $3 - $2}}' | \
            sort -nr | head -1`
        tmp_fw="{params.handle_tmp}.fw.bed.tmp"
        tmp_rc="{params.handle_tmp}.rc.bed.tmp"
        tail -n+2 {params.f_fw} | sed 's/\"//g; s/,/\\t/g' | \
            awk "{{OFS=\\"\\t\\"; print \\$2,\\$3,\\$3 + ${{size_fw}}}}" \
            > ${{tmp_fw}} && \
        tail -n+2 {params.f_rc} | sed 's/\"//g; s/,/\\t/g' | \
            awk "{{OFS=\\"\\t\\"; print \\$2,\\$3,\\$3 + ${{size_rc}}}}" \
            > ${{tmp_rc}} && \
        ionice -c2 -n7 zcat {input.motif} | \
            intersectBed -c -a stdin -b ${{tmp_fw}} ${{tmp_rc}} | \
            sed 's/\\t[1-9]$/\\t1/g' > {output} && \
        rm ${{tmp_fw}} ${{tmp_rc}}
        """

# process_piq2 is taking the ".all.csv" files from PIQ and extracting purity
# scores (estimated PPV) and then intersecting with the original motif file.
# When there are ties, the highest score will be used.
rule process_piq2:
    input:
        motif = os.path.join(MOTIF_DIR, "{motif}.bed.gz"),
        ok = rules.piq.output,
    output:
        os.path.join(
            PIQ_DIR, "processed_output_scores", "{sample}", 
            "{depth}M.{motif}.bed"
        ),
    params:
        handle_tmp = "scores_{sample}__{depth}__{motif}",
        f_fw = lambda wildcards:
            join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}-calls.all.csv".format(get_piq_handle(wildcards.motif))
            ),
        f_rc = lambda wildcards:
            join(
                PIQ_DIR, "out", "{}".format(wildcards.sample),
                "{}M".format(wildcards.depth), "out", 
                "{}.RC-calls.all.csv".format(get_piq_handle(wildcards.motif))
            ),
        f_size_fw = lambda wildcards:
            join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
        f_size_rc = lambda wildcards:
            os.path.join(
                PIQ_DIR, "out", "{}".format(wildcards.sample), 
                "{}M".format(wildcards.depth), "out", 
                "{}.RC-calls.all.bed".format(get_piq_handle(wildcards.motif))
            ),
    shell:
        """
        size_fw=`head -n 25 {params.f_size_fw} | awk '{{print $3 - $2}}' | \
            sort -nr | head -1`
        size_rc=`head -n 25 {params.f_size_rc} | awk '{{print $3 - $2}}' | \
            sort -nr | head -1`
        tmp_fw="{params.handle_tmp}.fw.bed.tmp"
        tmp_rc="{params.handle_tmp}.rc.bed.tmp"
        tail -n+2 {params.f_fw} | sed 's/\"//g; s/,/\\t/g' | \
            awk "{{OFS=\\"\\t\\"; print \\$2,\\$3,\\$3 + ${{size_fw}}, \
            \\".\\", \\$7, \\"+\\"}}" > ${{tmp_fw}} && \
        tail -n+2 {params.f_rc} | sed 's/\"//g; s/,/\\t/g' | \
            awk "{{OFS=\\"\\t\\"; print \\$2,\\$3,\\$3 + ${{size_rc}}, \
            \\".\\", \\$7, \\"-\\"}}" > ${{tmp_rc}} && \
        ionice -c2 -n7 zcat {input.motif} | \
            intersectBed -s -loj -a stdin -b ${{tmp_fw}} ${{tmp_rc}} | \
            awk '{{OFS="\\t"; print $1,$2,$3,$4,$12,$6}}' | \
            sort -k 1,1 -k2,2n -k5,5nr | sort -u -k 1,1 -k2,2n > {output} && \
        rm ${{tmp_fw}} ${{tmp_rc}}
        """