"""Dummy version of Snakefile for inferring sample metadata with HTSinfer."""


localrules:
    dummy,
    all,


rule all:
    input:
        config["samples_out"],


rule dummy:
    input:
        config["samples"],
    output:
        config["samples_out"],
    shell:
        'cat {input} > {output}'

