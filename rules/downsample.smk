# downsample bam files and test imputation
import os
shell.prefix("module load " + " ".join(config["modules"]) + " || true;")


# Read TSV file: bam, sample, coverage
bam_paths = {}
coverage_vals = {}
samples = []
with open(config["downsamplebams"]) as f:
    for line in f:
        line = line.strip()
        if not line:  # skip empty lines
            continue
        bam, sample, cov = line.split("\t")
        bam_paths[sample] = bam
        coverage_vals[sample] = cov
        samples.append(sample)
fractions = [str(f) for f in config["downsamplefractions"]]

# ---- final outputs ----
rule all_downsample:
    input:
        expand(["downsample/bams/{sample}_{fraction}.bam","downsample/bams/{sample}_{fraction}.coverage"],sample=samples,fraction=fractions)

# ---- downsample bams ----
rule downsample_bam:
    input:
        bam=lambda wc: bam_paths[wc.sample]
    output:
        bam="downsample/bams/{sample}_{fraction}.bam",
        cov="downsample/bams/{sample}_{fraction}.coverage"
    params:
        prop=lambda wc: float(wc.fraction) / float(coverage_vals[wc.sample])
    resources:
      mem_mb=50000,
      runtime=2*60
    shell:
        """
        [ -f {output.bam} ] && rm {output.bam}
        picard DownsampleSam -I {input.bam} -O {output.bam} -P {params.prop} -CREATE_INDEX True -STRATEGY Chained
        paleomix coverage --overwrite-output {output.bam} {output.cov}
        tail -n +22 {output.cov} | head -1 | awk '{{print "{wildcards.sample}", $NF}}' >> downsample/bams/coverage.txt
        """