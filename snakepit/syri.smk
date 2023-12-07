from pathlib import PurePath
from itertools import product

localrules: generate_genomes

if 'assemblies' not in config:
    config['assemblies'] = {'Ref':config['reference']} | {f'{sample}_{haplotype}':f'{haplotype}.scaffolds.fasta' for (sample,haplotype) in product(config['samples'],('hap1','hap2'))}

#wildcard_constraints:
#    asm = r'Ref|' + r'|'.join([f'{sample}_{haplotype}' for (sample,haplotype) in product(config['samples'],('hap1','hap2'))]),
#    A = r'Ref|' + r'|'.join([f'{sample}_{haplotype}' for (sample,haplotype) in product(config['samples'],('hap1','hap2'))]),
#    B = r'Ref|' + r'|'.join([f'{sample}_{haplotype}' for (sample,haplotype) in product(config['samples'],('hap1','hap2'))])

rule all:
    input:
        'syri.png',

rule samtools_faidx:
    input:
        lambda wildcards: config['assemblies'][wildcards.asm]
    output:
        'syri/{asm}.fa'
    params:
        regions = ' '.join(list(map(str,range(1,30))))
    shell:
        'samtools faidx --continue {input} {params.regions} | seqtk rename - chr > {output}'

rule minimap2_align:
    input:
        A = 'syri/{A}.fa',
        B = 'syri/{B}.fa'
    output:
        'syri/{A}_{B}.paf'
    threads: 4
    resources:
        mem_mb = 10000
    shell:
        'minimap2 -t {threads} -cxasm5 --cs --eqx {input.A} {input.B} > {output}'

rule syri_annotate:
    input:
        A = 'syri/{A}.fa',
        B = 'syri/{B}.fa',
        paf = 'syri/{A}_{B}.paf'
    output:
        'syri/{A}_{B}.syri.out'
    params:
        prefix = lambda wildcards, output: PurePath(output[0]).with_suffix('').with_suffix('').name + '.',
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    threads: 1
    resources:
        mem_mb = 25000
    conda:
        'syri'
    shell:
        '''
        syri -c {input.paf} --hdrseq --all -F P -r {input.A} -q {input.B} --dir {params._dir} --prefix {params.prefix}
        '''

rule generate_genomes:
    output:
        'syri/genomes.txt'
    run:
        with open(output[0],'w') as fout:
            fout.write('\t'.join(('#file','name','tags'))+'\n')
            for name,asm in config['assemblies'].items():
                fout.write('\t'.join((f'syri/{name}.fa',name,'lw:1.5'))+'\n')

def chain_input(wildcards):
    assemblies = list(config['assemblies'])
    for A,B in zip(assemblies,assemblies[1:]):
        yield f'syri/{A}_{B}.syri.out'

rule plotsr:
    input:
        outs = chain_input,
        genomes = rules.generate_genomes.output[0]
    output:
        'syri.png'
    params:
        SRs = lambda wildcards, input: ' '.join((f'--sr {out}' for out in input.outs))
    conda:
        'syri'
    threads: 1
    resources:
        mem_mb = 7500
    shell:
        '''
        plotsr \
        {params.SRs} \
        --genomes {input.genomes} \
        -o {output}
        '''
