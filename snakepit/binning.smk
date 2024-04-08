
rule all:
    input:
        expand('triobin/{sample}.tagged_reads.csv.gz',sample=config['samples']),
        expand('kmers/{sample}_genomescope2',sample=config['samples'])

rule get_tagged_alignments:
    input:
        cram = '/cluster/work/pausch/alex/NEW_ASSEMBLY/alignments/{sample}.mm2.cram' 
    output:
        'triobin/{sample}.tagged_reads.csv.gz'
    resources:
        mem_mb = 5000
    shell:
        '''
        samtools view -F 2304 --reference {config[reference]} {input.cram} |\
        awk 'BEGIN {{print "chromosome position flag HP"}} {{for(i=1;i<=NF;i++) {{if($i ~/HP:i:/) {{$i=substr($i,6);print $3,$4,$2,$i}} }} }}' |\
        pigz -p 2 -c > {output}
        '''

rule meryl_count:
    input:
        '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long_reads/{sample}/alignment/{sample}.mm2.cram'
    output:
        directory('kmers/{sample}.meryl')
    threads: 24
    resources:
        mem_mb = 3000
    shell:
        '''
        samtools fastq -@ {threads} --reference {config[reference]} -0 $TMPDIR/reads.fq.gz {input}
        meryl count k=21 threads={threads} $TMPDIR/reads.fq.gz output {output}
        '''

rule meryl_histogram:
    input:
        rules.meryl_count.output
    output:
        'kmers/{sample}.hist'
    localrule: True
    shell:
        '''
        meryl histogram {input} > {output}
        '''

rule genomescope2:
    input:
        rules.meryl_histogram.output
    output:
        directory('kmers/{sample}_genomescope2')
    conda: 'R'
    shell:
        '''
        /cluster/work/pausch/alex/software/genomescope2.0/genomescope.R -i {input} -o {output} -k 21
        '''
