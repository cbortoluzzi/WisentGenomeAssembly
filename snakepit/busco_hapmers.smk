#bedtools genomecov -i hap1.hap1_hapmers.bed -g ../../REF_DATA/ARS-UCD1.2_Btau5.0.1Y.fa.fai -bg > hap1.hap1_hapmers.wig
from pathlib import PurePath

rule all:
    input:
        expand('busco/{sample}/summary.txt',sample=config['assemblies'])

rule minibusco_download:
    output:
        directory(config.get('minibusco_download','/cluster/work/pausch/alex/REF_DATA/BUSCO'))
    localrule: True
    shell:
        '''
        minibusco download cetartiodactyla -L {output}
        '''

rule minibusco_run:
    input:
        fasta = lambda wildcards: config['assemblies'][wildcards.sample],
        db = rules.minibusco_download.output
    output:
        summary = 'busco/{sample}/summary.txt',
        table = 'busco/{sample}/cetartiodactyla_odb10/full_table_busco_format.tsv',
        _dir = directory('busco/{sample}')
    threads: 4
    resources:
        mem_mb = 12000
    shell:
        '''
        minibusco run -a {input.fasta} -o {output._dir} -l cetartiodactyla -L {input.db} -t {threads}
        '''

rule scrape_complete_buscos:
    input:
        rules.minibusco_run.output['table']
    output:
        'busco/{sample}.complete.bed'
    localrule: True
    shell:
        '''
        awk -v OFS='\t' '$2=="Complete" {{print $3,$4,$5,$1}}' {input} | bedtools sort -i - > {output}
        '''


rule sort_missing_busco:
    input:
        expand(rules.scrape_complete_buscos.output,sample=config['assemblies'])
    output:
        'busco/complete.upset'
    localrule: True
    shell:
        '''
        bedtools multiinter -cluster -i {input} > {output}
        '''

