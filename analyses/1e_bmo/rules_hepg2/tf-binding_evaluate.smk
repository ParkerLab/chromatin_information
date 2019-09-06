import os
import sys
import glob

rule concatenate_predictions:
    input:
        motif = os.path.join(MOTIF_DIR, "{motif}.bed.gz"),
        bmo = os.path.join(
            BMO_DIR, "{sample}", "all", "{motif}.{depth}.all.bed"
        ),
        centipede = rules.centipede.output,
        d2tf = rules.intersect_with_motifs_d2tf.output,
        hint = rules.intersect_with_motifs_hint.output,
        centipede_sb = rules.centipede_sb.output,
        peaks = rules.motifs_in_peaks2.output,
    output:
        os.path.join(CONCAT_DIR, "{sample}", "{motif}.{depth}.bed")
    params:
        motif = "tmp_{sample}.{motif}.{depth}.motif.tmp",
        bmo = "tmp_{sample}.{motif}.{depth}.bmo.tmp",
        centipede = "tmp_{sample}.{motif}.{depth}.centi.tmp",
        d2tf = "tmp_{sample}.{motif}.{depth}.d2tf.tmp",
        hint = "tmp_{sample}.{motif}.{depth}.hint.tmp",
        centipede_sb = "tmp_{sample}.{motif}.{depth}.centi_sb.tmp",
        peaks = "tmp_{sample}.{motif}.{depth}.peaks.tmp",
    resources:
        io_limit = 1
    shell:
        """
        ionice -c2 -n7 zcat {input.motif} > {params.motif} &&
        ionice -c2 -n7 cut -f 13 {input.bmo} > {params.bmo} &&
        ionice -c2 -n7 zcat {input.centipede} | cut -f 7 \
            > {params.centipede} &&
        ionice -c2 -n7 cut -f 7 {input.d2tf} > {params.d2tf} &&
        ionice -c2 -n7 cut -f 7 {input.hint} > {params.hint} &&
        ionice -c2 -n7 zcat {input.centipede_sb} | cut -f 7 \
            > {params.centipede_sb} &&
        ionice -c2 -n7 cut -f 7 {input.peaks} > {params.peaks} &&
        ionice -c2 -n7 paste {params.motif} {params.bmo} {params.centipede} \
            {params.d2tf} {params.hint} {params.centipede_sb} {params.peaks} \
            > {output}  &&
        rm {params.motif} {params.bmo} {params.centipede} {params.d2tf} \
            {params.hint} {params.centipede_sb} {params.peaks}
        """

rule make_bin_names_file:
    output:
        roc = prefix_results("bin_names.txt")
    run:
        print(
            *["BMO", "CENTIPEDE", "DNase2TF", "HINT", "CENTIPEDE_sb", 
              "Peaks_log10"], 
            sep = "\t", file = open(output[0], 'w')
        )

rule get_roc:
    input:
        associations = config["chip_assoc"],
        bins = rules.make_bin_names_file.output
    output:
        ok = prefix_results("roc", "{sample}", "get_roc.ok")
    params:
        prefix_results("roc")
    run:
        with open(input[0], 'r') as associations:
            with open("sh.roc", 'w') as sh_file:
                header = associations.readline()
                for line in associations.readlines():
                    line = line.strip().split("\t")
                    motif  = line[7]
                    tp_name  = line[8]
                    tp_file  = join(TP_DIR, "{}.bed".format(tp_name))
                    for sample, depths in subsample_depths.items():
                        for depth in depths:
                            infile = prefix_results("concatenated_files", "{}/{}.{}.bed".format(sample,motif,depth))
                            outdir = "{}/{}/{}".format(params[0],sample,depth)
                            shell("mkdir -p {outdir}")
                            # print(infile)
                            if os.path.isfile(infile):
                                shell_cmd = "ionice -c2 -n7 Rscript ../..bin/get_ROC_all.R {} 7 {} {} {} {}".format(
                                    infile, tp_name, tp_file, input[1], outdir
                                )
                                # print(shell_cmd)
                                print(shell_cmd, file=sh_file)
                            else:
                                print("get_roc - Could not find {}".format(infile), file=sys.stderr)
        shell("touch {output}")

rule get_f1:
    input:
        associations = config["chip_assoc"],
        bins = rules.make_bin_names_file.output
    output:
        ok = prefix_results("f1/{sample}/get_f1.ok")
    params:
        prefix_results("f1")
    run:
        with open(input[0], 'r') as associations:
            with open("sh.f1", 'w') as sh_file:
                for line in associations.readlines():
                    line = line.strip().split("\t")
                    motif  = line[7]
                    tp_name  = line[8]
                    tp_file  = join(TP_DIR, "{}.bed".format(tp_name))
                    for sample, depths in subsample_depths.items():
                        for depth in depths:
                            infile = prefix_results("concatenated_files/{}/{}.{}.bed".format(sample,motif,depth))
                            peaks  = prefix_results("motifs_in_atac-seq_peaks/{}/{}.{}.bed".format(sample,motif,depth))
                            outdir = "{}/{}/{}".format(params[0],sample,depth)
                            shell("mkdir -p {outdir}")
                            if os.path.isfile(infile):
                                shell_cmd = "ionice -c2 -n7 Rscript ../../bin/calculate_F1.R -f {} -c 7 -n {} -t {} -b {} -p {} -o {}".format(
                                    infile, tp_name, tp_file, input[1], peaks, outdir
                                )
                                print(shell_cmd, file=sh_file)
                            # else:
                            #     print("get_f1 - Could not find {}".format(infile), file=sys.stderr)
        shell("touch {output}")


rule aggregate_f1:
    input:
        rules.get_f1.output
    output:
        prefix_results("f1/{sample}/f1_data.txt")
    run:
        with open(output[0],"w") as out:
            header = ["sample", "depth", "precision", "recall", "f1", 
                      "method", "factor", "total_true_positives"]
            print(*header, sep="\t", file=out)
            for depth in subsample_depths[wildcards.sample]:
                f_pattern = prefix_results("f1/{}/{}/*.out".format(wildcards.sample, depth))
                for f in glob.glob(f_pattern):
                    with open(f, "r") as infile:
                        for line in infile.readlines():
                            extra_info = [wildcards.sample, depth]
                            line = extra_info + line.strip().split("\t")
                            print(*line, sep="\t", file=out)

rule aggregate_roc:
    input:
        rules.get_roc.output
    output:
        prefix_results("roc/{sample}/roc_data.txt")
    run:
        with open(output[0],"w") as out:
            header = ["sample", "depth", "factor", "BMO", "CENTIPEDE", 
                      "DNase2TF", "HINT", "CENTIPEDE_sb", "Peaks_log10",
                      "ChIPseq_size"]
            print(*header, sep="\t", file=out)
            for depth in subsample_depths[wildcards.sample]:
                f_pattern = prefix_results("roc/{}/{}/*.1.0.out".format(wildcards.sample, depth))
                for f in glob.glob(f_pattern):
                    tf_name = os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] # it's late & too lazy for regex tonight
                    with open(f, "r") as infile:
                        for line in infile.readlines():
                            if not line.startswith("BMO"):
                                extra_info = [wildcards.sample, depth, tf_name]
                                line = extra_info + line.strip().split("\t")
                                print(*line, sep="\t", file=out)

rule aggregate_prauc:
    input:
        rules.get_roc.output
    output:
        prefix_results("roc/{sample}/prauc_data.txt")
    run:
        with open(output[0],"w") as out:
            header = ["sample", "depth", "factor", "BMO", "CENTIPEDE", 
                      "DNase2TF", "HINT", "CENTIPEDE_sb", "Peaks_log10",
                      "ChIPseq_size"]
            print(*header, sep="\t", file=out)
            for depth in subsample_depths[wildcards.sample]:
                f_pattern = prefix_results("roc/{}/{}/*.prauc.out".format(wildcards.sample, depth))
                for f in glob.glob(f_pattern):
                    tf_name = os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] # it's late & too lazy for regex tonight
                    with open(f, "r") as infile:
                        for line in infile.readlines():
                            if not line.startswith("BMO"):
                                extra_info = [wildcards.sample, depth, tf_name]
                                line = extra_info + line.strip().split("\t")
                                print(*line, sep="\t", file=out)