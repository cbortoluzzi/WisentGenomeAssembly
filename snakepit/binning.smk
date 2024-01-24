
rule all:
    input:
        expand('triobin/{sample}.tagged_reads.csv.gz',sample=config['samples'])

rule get_tagged_alignments:
    input:
        cram = '/cluster/work/pausch/alex/NEW_ASSEMBLY/alignments/{sample}.mm2.cram' #'/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long_reads/{sample}/alignment/{sample}.mm2.cram'
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

