SAMPLES = [f"SHi26-{i}" for i in range(1, 13)]
NORMAL_SAMPLES = [f"SHi26-{i}" for i in range(3, 7)]

REF_FASTA = "/mnt/speedy/aboylan/ctDNA_2025/ref_genome_mm39/mm39.simple.fa"
GC_WIG = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/gc_mm39_1000kb.wig"
CENTROMERE_TXT = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/mm39_centromere_UCSC-gapTable.txt"

GENMAP_IDX = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_mm39"
GENMAP_RAW = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/100bp_E2.bedgraph"
GENMAP_FILTERED = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/100bp_E2.filtered.bedgraph"
WINDOWS_BED = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/windows_1Mb.bed"
BINNED_MAP_BED = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/genmap_raw/mappability_1Mb.bed"
MAP_WIG = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ref/map_mm39_1000kb.wig"

NORMAL_WIG_LIST = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/normal_wig_files.txt"

PON_BASE  = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/pon_output/ichorCNA_PON"
PON_RDS   = PON_BASE + "_median.rds"
PON_TXT   = PON_BASE + "_median.txt"

SAMPLE2GROUP = {
    "SHi26-1":  "group1",
    "SHi26-2":  "group2",
    "SHi26-3":  "group3",
    "SHi26-4":  "group3",
    "SHi26-5":  "group4",
    "SHi26-6":  "group4",
    "SHi26-7":  "group5",
    "SHi26-8":  "group5",
    "SHi26-9":  "group6",
    "SHi26-10": "group6",
    "SHi26-11": "group6",
    "SHi26-12": "group6",
}

NORMAL_GRID = {
    "group1": 'c(0.02,0.05,0.1)',
    "group2": 'c(0.02,0.05,0.1)',
    "group3": 'c(0.97)',
    "group4": 'c(0.97)',
    "group5": 'c(0.85)',
    "group6": 'c(0.85)',
}

rule all:
    input:
        expand("/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/{sample}.wig", sample=SAMPLES),
        GC_WIG,
        CENTROMERE_TXT,
        GENMAP_IDX,
        GENMAP_RAW,
        GENMAP_FILTERED,
        WINDOWS_BED,
        BINNED_MAP_BED,
        MAP_WIG,
        NORMAL_WIG_LIST,
        PON_RDS,
        PON_TXT,
        expand("/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ichor_out/{sample}", sample=SAMPLES)


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
        txt = CENTROMERE_TXT
    shell:
        """
        mkdir -p $(dirname {output.txt})
        (
          # header line ─ four tab‑separated columns
          echo -e "Chr\tStart\tEnd\tGapType"
          curl -s "http://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/gap.txt.gz" \
            | gunzip -c \
            | awk -F'\\t' '$8=="centromere" {{print $2"\\t"$3"\\t"$4"\\tcentromere"}}'
        ) > {output.txt}

        # drop any leading "chr" strings
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
        fasta = REF_FASTA
    output:
        bed = WINDOWS_BED
    shell:
        """
        # build chrom.sizes from the .fai index
        samtools faidx {input.fasta}
        cut -f1,2 {input.fasta}.fai \
          | grep -E '^(1[0-9]|[1-9]|X|Y)\t' \
          | sort -k1,1 -k2,2n > chrom.sizes

        bedtools makewindows -g chrom.sizes -w 1000000 > {output.bed}
        rm chrom.sizes
        """

rule bin_mappability_1Mb:
    input:
        windows = WINDOWS_BED,
        bedgraph = GENMAP_FILTERED
    output:
        bed = BINNED_MAP_BED
    shell:
        """
        bedtools map -a {input.windows} -b {input.bedgraph} -c 4 -o mean \
          | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,($4=="."?0:$4)}}' \
          > {output.bed}
        """

rule bed_to_wig_mappability:
    input:
        bed = BINNED_MAP_BED
    output:
        wig = MAP_WIG
    shell:
        """
        awk 'BEGIN{{chrom=""}}
             {{
               if ($1!=chrom) {{
                 print "fixedStep chrom="$1" start=1 step=1000000 span=1000000";
                 chrom=$1
               }}
               print $4
             }}' {input.bed} > {output.wig}
        """

rule list_normal_wig_files:
    input:
        expand("/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/{sample}.wig", sample=NORMAL_SAMPLES)
    output:
        list = NORMAL_WIG_LIST
    shell:
        """
        ls {input} > {output.list}
        """

rule create_pon:
    input:
        filelist   = NORMAL_WIG_LIST,
        gcwig      = GC_WIG,
        centromere = CENTROMERE_TXT
    output:
        rds = PON_RDS,
        txt = PON_TXT
    shell:
        """
        mkdir -p $(dirname {output.rds})

        Rscript /mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ichorCNA/scripts/createPanelOfNormals.R \
          --filelist    {input.filelist} \
          --gcWig       {input.gcwig} \
          --centromere  {input.centromere} \
          --outfile     {PON_BASE} \
          --method      median \
          --genomeStyle NCBI \
          --chrs        'c(1:19,"X","Y")' \
          --chrNormalize 'c(1:19)'
        """

rule run_ichorCNA:
    input:
        wig = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/wig_output/{sample}.wig",
        pon = PON_RDS,
        gc  = GC_WIG,
        cen = CENTROMERE_TXT
    output:
        dir = directory("/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ichor_out/{sample}")
    params:
    	normal = lambda wc: NORMAL_GRID[SAMPLE2GROUP[wc.sample]],
    	outdir = "/mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ichor_out",
    	id     = lambda wc: wc.sample,
    	group  = lambda wc: SAMPLE2GROUP[wc.sample],
    	ploidy = lambda wc: '"c(2)"' if SAMPLE2GROUP[wc.sample] in ["group3", "group4"] else '"c(3)"',
    	est_normal = lambda wc: "True" if SAMPLE2GROUP[wc.sample] in ["group3", "group4"] else "True",
    	est_ploidy = lambda wc: "True" if SAMPLE2GROUP[wc.sample] in ["group3", "group4"] else "True",
    	est_sc = lambda wc: "False" if SAMPLE2GROUP[wc.sample] in ["group3", "group4"] else "True",
    	sc_states = lambda wc: "NULL" if SAMPLE2GROUP[wc.sample] in ["group3", "group4"] else '"c(1,3)"'

    shell:
    	r"""
    	mkdir -p {params.outdir}

    	Rscript /mnt/speedy/aboylan/ctDNA_2025/ichorCNA_2025_05_07/ichorCNA/scripts/runIchorCNA.R \
      	--id          {params.id} \
      	--WIG         {input.wig} \
      	--ploidy      {params.ploidy} \
      	--normal      "{params.normal}" \
      	--maxCN       7 \
      	--gcWig       {input.gc} \
      	--centromere  {input.cen} \
      	--normalPanel {input.pon} \
      	--includeHOMD False \
      	--chrs        'c(1:19)' \
      	--chrTrain    'c(1:19)' \
      	--estimateNormal       {params.est_normal} \
      	--estimatePloidy       {params.est_ploidy} \
      	--estimateScPrevalence {params.est_sc} \
      	--scStates    {params.sc_states} \
      	--outDir      {params.outdir}
    	"""

