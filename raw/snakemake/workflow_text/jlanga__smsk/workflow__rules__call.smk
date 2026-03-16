rule call_bcftools:
    """
    Generate VCF with genotypes and call SNPs with bcftools
    """
    input:
        fa=RAW + "genome.fa",
        bam=expand(MAP + "{sample}.sorted.bam", sample=SAMPLES),
        bai=expand(MAP + "{sample}.sorted.bam.bai", sample=SAMPLES),
    output:
        protected(CALL + "all.vcf"),
    log:
        CALL + "bcftools.log",
    conda:
        "call.yml"
    shell:
        """
        ( bcftools mpileup \
            -Ou \
            -f {input.fa} \
            {input.bam} \
        | bcftools call \
            -vmO z \
        > {output} \
        ) 2> {log}
        """


rule call:
    """
    Checkpoint rule. Generate all.vcf
    """
    input:
        CALL + "all.vcf",
