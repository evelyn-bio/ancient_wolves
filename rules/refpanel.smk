#snakemake --executor slurm --use-conda --use-env --workflow-profile slurm --dry-run
workdir: "/projects/psg/people/pkb156/AW"
configfile: "config.yaml"

#load modules
shell.prefix("module load " + " ".join(config["modules"]) + " || true; ")

rule all_refpanel:
    input:
        expand("refpanel/wolves/phased.all.info08.name.chr{chr}.sites.vcf.gz", chr=config["chromosomes"]),
        expand("refpanel/canid/ref-panel_chr{chr}_sample-snp_filltags_filter.phased.filter.vcf.gz", chr=config["chromosomes"]),
        expand("refpanel/canid/ref-panel_chr{chr}_sample-snp_filltags_filter.phased.filter.vcf.gz.csi", chr=config["chromosomes"]),
        expand("refpanel/wolves/phased.all.info08.name.noancient.chr{chr}.vcf.gz", chr=config["chromosomes"]),
        expand("refpanel/merged/mergedpanel.chr{chr}.vcf.gz", chr=config["chromosomes"]),
        "king/removedups.txt",
        expand("refpanel/merged/mergedpanel.nodups.chr{chr}.vcf.gz", chr=config["chromosomes"])
        
# ---- select sites from wolf panel ----
rule select_sites:
    input: 
      vcf=config["modern_wolves"]
    output: 
      sites="refpanel/wolves/phased.all.info08.name.chr{chr}.sites.vcf.gz"
    resources:
      mem_mb=20000,
      runtime=1*60
    threads: 4
    shell:
        """
        bcftools view --threads {threads} -G -r chr{wildcards.chr} {input.vcf} -Oz -o {output.sites}
        """
        
# ---- subset sites in canid panel ----
rule subset_sites:
    input: 
      canidvcf="canidref/ref-panel_chr{chr}_sample-snp_filltags_filter.phased.vcf.gz",
      wolfvcf="refpanel/wolves/phased.all.info08.name.chr{chr}.sites.vcf.gz"
    output: 
      vcf="refpanel/canid/ref-panel_chr{chr}_sample-snp_filltags_filter.phased.filter.vcf.gz",
      vcfidx="refpanel/canid/ref-panel_chr{chr}_sample-snp_filltags_filter.phased.filter.vcf.gz.csi"
    resources:
      mem_mb=20000,
      runtime=2*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads {threads} -R {input.wolfvcf} {input.canidvcf} | bcftools view --threads {threads} -m2 -M2 -v snps - -o {output.vcf}
        bcftools index -f {output.vcf}
        """
        
# ---- remove ancients from wolf panel ----
rule remove_ancient:
    input:
        vcf=config["modern_wolves"],
        ancient=config["ancientwolflist"]
    output: 
        vcf="refpanel/wolves/phased.all.info08.name.noancient.chr{chr}.vcf.gz",
        vcfidx="refpanel/wolves/phased.all.info08.name.noancient.chr{chr}.vcf.gz.csi",
        vcfstats="refpanel/wolves/phased.all.info08.name.noancient.chr{chr}.stats"
    resources:
        mem_mb=20000,
        runtime=2*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads {threads} -S ^{input.ancient} -r chr{wildcards.chr} --force-samples {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """

# ---- merge two panels ----
rule merge_panel:
    input:
        wolf="refpanel/wolves/phased.all.info08.name.noancient.chr{chr}.vcf.gz",
        canid="refpanel/canid/ref-panel_chr{chr}_sample-snp_filltags_filter.phased.filter.vcf.gz"
    output: 
        vcf="refpanel/merged/mergedpanel.chr{chr}.vcf.gz",
        vcfidx="refpanel/merged/mergedpanel.chr{chr}.vcf.gz.csi",
        vcfstats="refpanel/merged/mergedpanel.chr{chr}.stats"
    resources:
        mem_mb=30000,
        runtime=8*60
    threads: 4
    shell:
        """
        bcftools merge -Oz --force-samples --threads {threads} -r chr{wildcards.chr} {input.wolf} {input.canid} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """
         
# ---- identify duplicate samples ----
rule king_duplicates:
    input: 
        vcf="refpanel/merged/mergedpanel.chr1.vcf.gz"
    output: 
        dups="king/removedups.txt",
        bed="king/mergedpanel.chr1.bed",
        bim="king/mergedpanel.chr1.bim",
        fam="king/mergedpanel.chr1.fam"
    resources:
        mem_mb=80000,
        runtime=1*60
    shell:
        """
        plink --vcf {input.vcf} --double-id --dog --make-bed --out king/mergedpanel.chr1
        king -b king/mergedpanel.chr1.bed --related --sexchr 39 --prefix king/mergedpanel.chr1 --degree 2
        awk '$14=="Dup/MZ"{{print $4}}' king/mergedpanel.chr1.kin0 > {output.dups}
        """
        
# ---- remove duplicates from vcf ----
rule remove_dups:
    input: 
      vcf="refpanel/merged/mergedpanel.chr{chr}.vcf.gz",
      rmfile="king/removedups.txt"
    output: 
      vcf="refpanel/merged/mergedpanel.nodups.chr{chr}.vcf.gz",
      vcfidx="refpanel/merged/mergedpanel.nodups.chr{chr}.vcf.gz.csi",
      vcfstats="refpanel/merged/mergedpanel.nodups.chr{chr}.stats"    
    resources:
      mem_mb=30000,
      runtime=6*60
    threads: 4
    shell:
        """
        bcftools view -Oz --threads 4 -S ^{input.rmfile} --force-samples -r chr{wildcards.chr} {input.vcf} -o {output.vcf}
        bcftools index -f {output.vcf}
        bcftools stats {output.vcf} > {output.vcfstats}
        """