# Define all sample IDs
SAMPLES = [f"SHi26-{i}" for i in range(1, 13)]

REF_FASTA = "/mnt/speedy/aboylan/ctDNA_2025/ref_genome_mm39/mm39.simple.fa"
GC_WIG = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/gc_mm39_1000kb.wig"
CENTROMERE_TXT = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/mm39_centromere_UCSC-gapTable.txt"

GENMAP_IDX = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_mm39"
GENMAP_RAW = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/100bp_E2.bedgraph"

MAP_WIG = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/map_mm39_1000kb.wig"

# Master rule: build all .wig files
rule all:
    input:
        expand("/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/{sample}.wig", sample=SAMPLES),
        GC_WIG,
        CENTROMERE_TXT,
        GENMAP_IDX,
        GENMAP_RAW,
        MAP_WIG

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
        echo -e "Y\t110000\t3000000\tcentromere" >> {output.txt}
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
            -bg                    # ← use the short flag

        mv $(dirname {output.bed})/*.bedgraph {output.bed}
        """

rule genmap_bin1Mb:
    input:
        bed = GENMAP_RAW
    output:
        wig = MAP_WIG
    shell:
        """
        # make 1 Mb windows per chromosome
        bedtools makewindows -g <(cut -f1,2 {input.bed} | sort -k1,1 -k2,2n \
                                 | awk '{{max[$1]>$2?max[$1]:max[$1]=$2}} END{{for(c in max) print c"\\t"max[c]}}') \
                             -w 1000000 > windows.bed

        # average mappability in each window
        bedtools map -a windows.bed -b {input.bed} -c 4 -o mean \
        | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$4=="."?0:$4}}' \
        | awk 'BEGIN{{chrom=""}} {{if($1!=chrom){{print "fixedStep chrom="$1" start=1 step=1000000"; chrom=$1}} print $4}}' \
        > {output.wig}

        rm windows.bed
        """
