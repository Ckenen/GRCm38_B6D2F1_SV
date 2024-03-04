#!/usr/bin/env runsnakemake
configfile: "config.yaml"
BUILD = config["BUILD"]
SAMPLE = config["SAMPLE"]
METHOD = config["METHOD"]
NAME = "%s_%s_%s" % (BUILD, SAMPLE, METHOD)
FASTA = config["FASTA"]
FASTQ = config["FASTQ"]
CHROMS = config["CHROMS"]
MM2_PRESET = config["MM2_PRESET"]
THREADS = config["THREADS"]
OUTDIR = config["OUTDIR"]


rule all:
    input:
        OUTDIR + "/%s.trimmed.fastq.gz" % NAME,
        OUTDIR + "/%s.primary.fa" % BUILD,
        expand(OUTDIR + "/tandemRepeats/%s.{chrom}.bed" % BUILD, chrom=CHROMS),
        OUTDIR + "/%s.tandemRepeats.bed" % BUILD,
        OUTDIR + "/%s.%s.mmi" % (BUILD, MM2_PRESET),
        OUTDIR + "/%s.mm2.bam" % NAME,
        OUTDIR + "/%s.mm2.flagstat" % NAME,
        OUTDIR + "/%s.sniffles2.vcf.gz" % NAME,


rule cutadapt:
    input:
        fq = FASTQ
    output:
        fq1 = temp(OUTDIR + "/%s.trimmed.tmp.fastq" % NAME),
        fq2 = OUTDIR + "/%s.trimmed.fastq.gz" % NAME
    log:
        OUTDIR + "/%s.trimmed.log" % NAME
    threads:
        THREADS
    shell:
        """(
        cutadapt -j {threads} -e 0.2 -m 400 -g TCGTTCAGTTACGTATTGCT -o {output.fq1} {input.fq}
        cutadapt -j {threads} -e 0.2 -m 400 -a AGCAATACGTAACTGAACGA -o {output.fq2} {output.fq1} ) &> {log}
        """

rule make_fasta:
    input:
        fa = FASTA
    output:
        fa = OUTDIR + "/%s.primary.fa" % BUILD
    shell:
        """
        samtools faidx {input.fa} {CHROMS} > {output.fa}
        samtools faidx {output.fa}
        """

rule find_tandem_repeats:
    input:
        fa = rules.make_fasta.output.fa
    output:
        bed = OUTDIR + "/tandemRepeats/%s.{chrom}.bed" % BUILD
    log:
        OUTDIR + "/tandemRepeats/%s.{chrom}.log" % BUILD
    shell:
        """
        findTandemRepeats --merge --chrom {wildcards.chrom} {input.fa} {output.bed} &> {log}
        """

rule merge_randem_repeats:
    input:
        beds = expand(rules.find_tandem_repeats.output.bed, chrom=CHROMS)
    output:
        bed = OUTDIR + "/%s.tandemRepeats.bed" % BUILD
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n > {output.bed}
        bgzip -c {output.bed} > {output.bed}.gz
        tabix -p bed {output.bed}.gz
        """

rule build_index:
    input:
        fa = rules.make_fasta.output.fa
    output:
        mmi = OUTDIR + "/%s.%s.mmi" % (BUILD, MM2_PRESET)
    log:
        OUTDIR + "/%s.%s.log" % (BUILD, MM2_PRESET)
    threads:
        THREADS
    shell:
        """
        minimap2 -t {threads} -x {MM2_PRESET} -d {output.mmi} {input.fa} &> {log} 
        """

rule minimap2:
    input:
        fq = rules.cutadapt.output.fq2,
        mmi = rules.build_index.output.mmi
    output:
        bam = OUTDIR + "/%s.mm2.bam" % NAME
    log:
        OUTDIR + "/%s.mm2.log" % NAME
    params:
        rg = "@RG\\tID:%s\\tLB:%s\\tSM:%s" % (SAMPLE, SAMPLE, SAMPLE)
    threads:
        THREADS
    shell:
        """(
        minimap2 -ax {MM2_PRESET} --MD -t {threads} -R '{params.rg}' \
            --secondary=no {input.mmi} {input.fq} \
            | samtools view -@ {threads} -u -F 4 - \
            | samtools sort -@ {threads} -T {output.bam}_TMP -o {output.bam} - 
        samtools index -@ {threads} {output.bam} ) &> {log}
        """

rule sniffles2:
    input:
        bam = rules.minimap2.output.bam,
        fa = rules.make_fasta.output.fa,
        bed = rules.merge_randem_repeats.output.bed
    output:
        tmp1 = temp(OUTDIR + "/%s.sniffles2.vcf" % NAME),
        vcf1 = OUTDIR + "/%s.sniffles2.vcf.gz" % NAME,
        tmp2 = temp(OUTDIR + "/%s.sniffles2.filtered.vcf" % NAME),
        vcf2 = OUTDIR + "/%s.sniffles2.filtered.vcf.gz" % NAME
    log:
        OUTDIR + "/%s.sniffles2.log" % NAME
    threads:
        THREADS
    shell:
        """(
        set +u; source activate sniffles2
        sniffles -t {threads} --phase --output-rnames \
            --minsvlen 20 --sample-id {SAMPLE} \
            --tandem-repeats {input.bed} --reference {input.fa} \
            -i {input.bam} -v {output.tmp1}
        bgzip -c {output.tmp1} > {output.vcf1}
        tabix -p vcf -f {output.vcf1}
        cat {output.tmp1} | grep '#' > {output.tmp2}
        cat {output.tmp1} | grep -v '#' | grep -v 'INV' \
            | grep -v 'DUP' | grep -v 'BND' >> {output.tmp2}
        bgzip -c {output.tmp2} > {output.vcf2}
        tabix -p vcf -f {output.vcf2} ) &> {log}
        """

rule flagstat:
    input:
        bam = "{prefix}.bam"
    output:
        txt = "{prefix}.flagstat"
    threads:
        4
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.txt}
        """
