#snakemake --executor slurm --use-conda --use-env --workflow-profile slurm --dry-run
workdir: "/projects/psg/people/pkb156/AW"
configfile: "config.yaml"

#load modules
shell.prefix("module load " + " ".join(config["modules"]) + " || true; ")

# Master Snakefile
include: "rules/refpanel.smk"
include: "rules/downsample.smk"

rule all:
    input:
        expand("refpanel/merged/mergedpanel.nodups.chr{chr}.vcf.gz", chr=config["chromosomes"]),
        rules.all_downsample.input