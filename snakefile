import os

SAMPLES = ["LS140", "E15", "L86"]
BINARY = ["binary", "binaryError"]
GT16 = ["gt16", "gt16Error"]

wildcard_constraints:
    # Global wildcard constraints, in a regex form; a|b means a OR b.
    # adapted from https://stackoverflow.com/a/59227356/4868692
    sample = '|'.join([x for x in SAMPLES]),
    template = '|'.join([x for x in BINARY+GT16])


rule all:
    input:
        expand("results/{sample}/{template}/done", sample=["E15", "L86"], template=BINARY+GT16),
        expand("results/{sample}/{template}/done", sample=["LS140"], template=BINARY)


rule vcf2fasta:
    input:
        "data/{sample}.vcf",
    output:
        "results/{sample}/{sample}_{template}.fasta"
    params:
        encoding = lambda w: "binary" if w.template in BINARY else "nd16"
    shell:
        "python src/vcf2fasta.py {input} {output} --encoding {params.encoding}"


rule fasta2stats:
    input:
        "results/{sample}/{sample}_{template}.fasta"
    output:
        stats = "results/{sample}/{sample}_{template}.stats.txt",
        fasta = "results/{sample}/{sample}_{template}.filtered.fasta"
    shell:
        "Rscript src/fasta2stats.r {input} {output.stats} --filtered {output.fasta}"


rule fasta2xml:
    input:
        fasta = "results/{sample}/{sample}_{template}.filtered.fasta",
        template = "templates/{template}.xml"
    output:
        "results/{sample}/{template}/{sample}_{template}.xml"
    params:
        datatype = lambda w: "nucleotideDiploid16" if w.template in GT16 else "standard"
    shell:
        "Rscript src/fasta2xml.r {input.template} {input.fasta} {output} --datatype {params.datatype}"


rule beast:
    input:
        "results/{sample}/{template}/{sample}_{template}.xml"
    output:
        "results/{sample}/{template}/{sample}_{template}.trace",
        "results/{sample}/{template}/{sample}_{template}.trees",
        "results/{sample}/{template}/{sample}_{template}.log"
    params:
        beast = os.path.abspath("bin/beast-phylonco.jar"),
        input = "{sample}_{template}.xml",
        log = "{sample}_{template}.log",
        dir = "results/{sample}/{template}"
    shell:
        "cd {params.dir} && java -jar {params.beast} {params.input} > {params.log}"


rule loganalyser:
    input:
        "results/{sample}/{template}/{sample}_{template}.trace"
    output:
        "results/{sample}/{template}/{sample}_{template}.ess.txt"
    shell:
        "loganalyser -b 20 {input} > {output}"


rule treeannotator:
    input:
        "results/{sample}/{template}/{sample}_{template}.trees"
    output:
        "results/{sample}/{template}/{sample}_{template}.tree"
    shell:
        "treeannotator -b 20 -lowMem {input} {output}"


rule done:
    input:
        tree = "results/{sample}/{template}/{sample}_{template}.tree",
        ess = "results/{sample}/{template}/{sample}_{template}.ess.txt",
        log = "results/{sample}/{template}/{sample}_{template}.log"
    output:
        "results/{sample}/{template}/done"
    shell:
        "touch {output}"
