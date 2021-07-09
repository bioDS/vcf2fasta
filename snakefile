SAMPLES = ["LS140", "E15", "L86"]

rule all:
    input: expand("results/{id}/{id}.tree", id=SAMPLES),
           expand("results/{id}/{id}_ess.txt", id=SAMPLES),
           expand("results/{id}_stats.txt", id=SAMPLES),


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
        "results/{sample}_filtered.fasta"
    output:
        "results/{sample}/{sample}.xml"
    shell:
        "Rscript src/fasta2xml.r templates/nucleotideDiploid16.xml {input} {output}"


rule beast:
    input:
        "results/{sample}/{sample}.xml"
    output:
        "results/{sample}/{sample}.trace",
        "results/{sample}/{sample}.trees"
    shell:
        "java -jar beast-phylonco.jar {input}"


rule loganalyser:
    input:
        "results/{sample}/{sample}.trace"
    output:
        "results/{sample}/{samples}_ess.txt"
    shell:
        "loganalyser -b 20 {input} > {output}"


rule treeannotator:
    input:
        "results/{sample}/{sample}.trees"
    output:
        "results/{sample}/{sample}.tree"
    shell:
        "treeannotator -b 20 -lowMem {input} {output}"
