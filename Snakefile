SAMPLES = [f"SHi26-{i}" for i in range(1, 13)]

REF_FASTA = "/mnt/speedy/aboylan/ctDNA_2025/ref_genome_mm39/mm39.simple.fa"
GC_WIG = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/gc_mm39_1000kb.wig"
CENTROMERE_TXT = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/mm39_centromere_UCSC-gapTable.txt"

GENMAP_IDX = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_mm39"
GENMAP_RAW = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/100bp_E2.bedgraph"
GENMAP_FILTERED = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/100bp_E2.filtered.bedgraph"
WINDOWS_BED = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/windows_1Mb.bed"


rule all:
    input:
        expand("/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/{sample}.wig", sample=SAMPLES),
        GC_WIG,
        CENTROMERE_TXT,
        GENMAP_IDX,
        GENMAP_RAW,
        GENMAP_FILTERED,
        WINDOWS_BED

# Run readCounter on each sample's deduplicated BAM
rule readcounter:
    input:
        bam="/mnt/speedy/aboylan/ctDNA_2025/experiment_data_2025_01_20/{sample}/mt_coverage_2025_04_21/{sample}_aligned.sorted.rg.markdup.sorted.simple.dedup.bam"
    output:
        wig="/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/{sample}.wig"
    shell:
        """
        mkdir -p $(dirname {output.wig})
        /mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/hmmcopy_utils/bin/readCounter \
          --window 1000000 \
          --quality 20 \
          --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y" \
          {input.bam} > {output.wig}
        """

rule gccounter:
    input:
        fasta=REF_FASTA
    output:
        wig=GC_WIG
    shell:
        """
        mkdir -p $(dirname {output.wig})
        /mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/hmmcopy_utils/bin/gcCounter \
          --window 1000000 \
          -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,X,Y \
          {input.fasta} > {output.wig}
        """

rule get_centromeres:
    output:
        txt=CENTROMERE_TXT
    shell:
        """
        mkdir -p $(dirname {output.txt})
        curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/gap.txt.gz" \
          | gunzip -c \
          | awk -F'\\t' '$8=="centromere" {{print $2"\\t"$3"\\t"$4"\\tcentromere"}}' > {output.txt}
        sed -i 's/^chr//' {output.txt}
        """


rule genmap_index:
    input:
        fasta = REF_FASTA
    output:
        idx_dir = directory(GENMAP_IDX)
    threads: 12
    shell:
        """
        genmap index -F {input.fasta} -I {output.idx_dir}
        """

rule genmap_calculate:
    input:
        idx_dir = GENMAP_IDX          # index folder
    output:
        bed = GENMAP_RAW
    threads: 12
    shell:
        """
        mkdir -p $(dirname {output.bed})

        genmap map \
            -K 100 \
            -E 2 \
            -I {input.idx_dir} \
            -O $(dirname {output.bed}) \
            -T {threads} \
            -bg                 

        mv $(dirname {output.bed})/*.bedgraph {output.bed}
        """

rule genmap_filter_canonical:
    input:
        bed = GENMAP_RAW
    output:
        bed = GENMAP_FILTERED
    shell:
        """
        grep -E '^(1[0-9]|[1-9]|X|Y)[[:space:]]' {input.bed} \
        | sort -k1,1 -k2,2n > {output.bed}
        """

rule make_windows_1Mb:
    input:
        bed = GENMAP_FILTERED
    output:
        bed = WINDOWS_BED
    shell:
        """
        cut -f1,3 {input.bed} \
          | sort -k1,1 -k2,2n \
          | awk '{{len[$1]<$2?len[$1]=$2:1}} END{{for(c in len) print c"\\t"len[c]}}' \
          | sort -k1,1 > chrom.sizes

        bedtools makewindows -g chrom.sizes -w 1000000 > {output.bed}
        rm chrom.sizes
        """



