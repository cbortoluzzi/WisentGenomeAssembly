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
        awk -v OFS='\t' '$2=="Complete" {{print $3,$4,$5,$1}}' {input} | sort -k1,1V -k2,2n > {output}
        '''


rule sort_missing_busco:
    input:
        expand(rules.scrape_complete_buscos.output,sample=config['assemblies'])
    output:
        'busco/complete.upset'
    localrule: True
    shell:
        '''

        LC_ALL=C
        #grep -Ff <(cut -f 4 {input[0]}) {input[1]} | sort -k1,1V -k2,2n > hap1_primary.match.bed
        #grep -vFf <(cut -f 4 {input[0]}) {input[1]} | sort -k1,1V -k2,2n > hap1_primary.missing.bed
        #grep -Ff <(cut -f 4 {input[2]}) {input[1]} | sort -k1,1V -k2,2n > hap2_primary.match.bed
        #grep -vFf <(cut -f 4 {input[2]}) {input[1]} | sort -k1,1V -k2,2n > hap2_primary.missing.bed
    
        grep -Ff <(cut -f 4 {input[0]}) {input[1]} | bedtools sort -i - > hap1_primary.match.bed
        grep -vFf <(cut -f 4 {input[0]}) {input[1]} | bedtools sort -i - > hap1_primary.missing.bed
        grep -Ff <(cut -f 4 {input[2]}) {input[1]} | bedtools sort -i - > hap2_primary.match.bed
        grep -vFf <(cut -f 4 {input[2]}) {input[1]} | bedtools sort -i - > hap2_primary.missing.bed

        comm -3 <(sort hap1_primary.match.bed) <(sort hap2_primary.match.bed) > {output}
        '''

rule odgi_pav:
    input:
        'graphs/{chromosome}.pggb.og'
    output:
        'graphs/{chromosome}.pggb.pav'
    shell:
        '''
        bedtools makewindows -g WIS_primary#0#17.fa.fai -w 100000 > windows.bed
        odgi pav -i {input} -t 1 -b windows.bed -M > {output}
        '''
