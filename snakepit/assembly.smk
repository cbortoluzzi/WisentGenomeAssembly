from pathlib import PurePath

rule all:
    input:
        #expand('mdbg_assembly/hap{haplotype}.fa',haplotype=(1,2))
        expand('hifiasm/{sample}/asm.bp.hap1.p_ctg.gfa',sample=(1,2))

rule samtools_fastq:
    input:
        '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/long_reads/WISENTM_Urano_F1/alignment/WISENTM_Urano_F1.mm2.cram'
    output:
        temp('reads/hap{sample}.fa.gz')
    params:
        reference = '/cluster/work/pausch/inputs/ref/BTA/UCD1.2/ARS-UCD1.2_Btau5.0.1Y.fa'
    shell:
        '''
        samtools fasta --reference {params.reference} -d HP:0 -d HP:{wildcards.sample} --threads {threads} -0 {output} {input}
        '''

rule hifiasm:
    input:
        rules.samtools_fastq.output
    output:
        gfa = expand('hifiasm/{{sample}}/asm.bp.{haplotype}.p_ctg.gfa',haplotype=('hap1','hap2')),
        intermediates = temp(multiext('hifiasm/{sample}/asm','.ec.bin','.ovlp.source.bin','.ovlp.reverse.bin','.bp.r_utg.gfa','.bp.p_utg.gfa'))
    params:
        out = lambda wildcards, output: PurePath(output['gfa'][0]).with_suffix('').with_suffix('').with_suffix('').with_suffix('')
    threads: 32
    resources:
        mem_mb = lambda wildcards, input, threads: max(int(input.size_mb*2.25/threads),8000),
        walltime = '4h'
    shell:
        'hifiasm -o {params.out} -t {threads} -r 3 -a 5 -n 5 {input}'
        
rule mdbg:
    input:
        rules.samtools_fastq.output
    output:
        'mdbg_assembly/hap{haplotype}.wa'
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 16
    resources:
        mem_mb = 1500,
        constraint = "EPYC_7H12"
    conda:
        'rusty'
    shell:
        '''
        /cluster/work/pausch/alex/software/rust-mdbg/target/release/rust-mdbg --threads {threads} -k 21 -d 0.003 -l 14 --prefix {params} {input}
        /cluster/work/pausch/alex/software/rust-mdbg/utils/magic_simplify {params}
        '''

