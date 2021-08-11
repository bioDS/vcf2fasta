import os

SAMPLES = ["LS140", "E15", "L86"]
BINARY = ["binary", "binaryError"]
ND16 = ["nd16", "nd16Error"]

rule all:
    input:
        expand("results/{sample}/{template}/done", sample=["E15", "L86"], template=BINARY+ND16),
        expand("results/{sample}/{template}/done", sample=["LS140"], template=BINARY)


rule vcf2fasta:
    input:
        "data/{sample}.vcf",
    output:
        "results/{sample}.fasta"
    shell:
        "python src/vcf2fasta.py {input} {output}"


rule fasta2stats:
    input:
        "results/{sample}.fasta"
    output:
        stats = "results/{sample}_stats.txt",
        fasta = "results/{sample}_filtered.fasta"
    shell:
        "Rscript src/fasta2stats.r {input} {output.stats} --filtered {output.fasta}"


rule fasta2xml:
    input:
        fasta = "results/{sample}_filtered.fasta",
        template = templates/{template}.xml
    output:
        "results/{sample}/{template}/{sample}_{template}.xml"
    params:
        if {template} in ["nd16", "nd16Error"]
            datatype = "--datatype nucleotideDiploid"
        else:
            datatype = ""
    shell:
        "Rscript src/fasta2xml.r {input.template} {input.fasta} {output} {params.datatype}"


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
        "results/{sample}/{template}/{sample}_{template}_ess.txt"
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
        ess = "results/{sample}/{template}/{sample}_{template}_ess.txt",
        log = "results/{sample}/{template}/{sample}_{template}.log"
    output:
        "results/{sample}/{template}/done"
    shell:
        "touch {output}"
