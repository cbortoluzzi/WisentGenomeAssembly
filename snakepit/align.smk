from pathlib import PurePath

local_samples, SRA_samples = [], {}
with open(config['samples'],'r') as fin:
    for line in fin:
        name, ID, *_ = line.split()
        if ID == 'local':
            local_samples.append(name)
        elif name[0] == '#':
            continue
        else:
            SRA_samples[name] = ID

samples = local_samples + list(SRA_samples.keys())

rule all:
    input:
        expand('publicSamples/{sample}.{reference}.cram',reference=('ARS',),sample=samples)


rule fastq_dl:
    output:
        temp(expand('publicSamples/fastq/{accession}_R{N}.fastq.gz', N=(1,2), allow_missing=True)),
    params:
        _dir = lambda wildcards, output: PurePath(output[0]).parent
    threads: 1
    resources:
        mem_mb = 5000,
        proxy_load = 1
    envmodules:
        'eth_proxy'
    conda: 'fastq-dl'
    shell:
        '''
        fastq-dl --accession {wildcards.accession} --cpus {threads} --silent --outdir {params._dir} --prefix {wildcards.accession} --group-by-sample && [[ -s {output[0]} ]] && [[ -s {output[1]} ]]
        '''

rule fastq_dl_fixed:
    input:
        fastq = 'publicSamples/fastq/{accession}_R{N}.fastq.gz'
    output:
        temp('publicSamples/fastq/{accession}_R{N}.fixed.fastq.gz')
    threads: 4
    resources:
        mem_mb = 1500
    shell:
        '''
        pigz -dc -p {threads} {input.fastq} | awk '{{print (NR%4==1 && NF>1) ? "@"$2  : $0}}' | pigz -c -p {threads} > {output}
        '''

rule fastp_filter:
    input:
        lambda wildcards: expand(rules.fastq_dl_fixed.output,accession=SRA_samples[wildcards.sample],N=(1,2),allow_missing=True) if wildcards.sample in SRA_samples else expand(config['local_bams']+'{{sample}}_R{N}.fastq.gz',N=(1,2))
    output:
        fastq = temp(expand('publicSamples/fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True))
    params:
        min_quality = config.get('fastp',{}).get('min_quality',15),
        unqualified = config.get('fastp',{}).get('unqualified',40),
        min_length  = config.get('fastp',{}).get('min_length',15)
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        '''
        fastp -q {params.min_quality} -u {params.unqualified} -g --length_required {params.min_length} --thread {threads} -i {input[0]} -o {output.fastq[0]} -I {input[1]} -O {output.fastq[1]} --json /dev/null --html /dev/null
        '''

def generate_SR_aligner_command(wildcards,input):
    match wildcards.mapper:
        case 'strobe':
            return f'strobealign {input.reference} {input.fastq} --rg-id {wildcards.sample}'
        case 'bwa':
            return f'bwa-mem2 mem -R "@RG\\tID:{wildcards.sample}\\tCN:UNK\\tLB:{wildcards.sample}\\tPL:illumina\\tSM:{wildcards.sample}" -Y {input.reference} {input.fastq}'
        case _:
            raise('Unknown aligner')

rule short_read_align:
    input:
        fastq = expand(rules.fastp_filter.output,allow_missing=True),
        reference = lambda wildcards: config['reference'][wildcards.reference]
    output:
        bam = multiext('publicSamples/{sample}.{reference}.{mapper,strobe|bwa}.cram','','.crai'),
        dedup_stats = 'publicSamples/{sample}.{reference}.{mapper}.dedup.stats'
    params:
        aligner_command = lambda wildcards, input: generate_SR_aligner_command(wildcards,input)
    threads: lambda wildcards: 24 if wildcards.mapper == 'bwa' else 12
    resources:
        mem_mb = 3000,
        scratch = '50g',
        walltime = '24h'
    shell:
        '''
        {params.aligner_command} -t {threads} |\
        samtools collate -u -O -@ {threads} - |\
        samtools fixmate -m -u -@ {threads} - - |\
        samtools sort -T $TMPDIR -u -@ {threads} |\
        samtools markdup -T $TMPDIR -S -@ {threads} --write-index -f {output.dedup_stats} --reference {input.reference} --output-fmt-option version=3.0 - {output.bam[0]}
        '''

rule strobealign:
    input:
        fastq = expand('publicSamples/fastq/{sample}.R{N}.fastq.gz',N=(1,2),allow_missing=True),
        reference = lambda wildcards: config['reference'][wildcards.reference]
    output:
        bam = temp(multiext('publicSamples/{sample}.{reference}.cram','','.crai')),
        dedup_stats = 'publicSamples/{sample}.{reference}.dedup.stats'
    params:
        rg = '"@RG\\tID:{sample}\\tCN:UNK\\tLB:{sample}\\tPL:illumina\\tSM:{sample}"',
    threads: 24
    resources:
        mem_mb = 2000,
        scratch = '50g',
        walltime = '4h'
    shell:
        '''
        strobealign {input.reference} {input.fastq} -t {threads} --rg-id {wildcards.sample} |\
        samtools collate -u -O -@ {threads} - |\
        samtools fixmate -m -u -@ {threads} - - |\
        samtools sort -T $TMPDIR -u -@ {threads} |\
        samtools markdup -T $TMPDIR -S -@ {threads} --write-index -f {output.dedup_stats} --reference {input.reference} --output-fmt-option version=3.1 - {output.bam[0]}
        '''

rule samtools_depth:
    input:
        cram = rules.strobealign.output['bam'],
        reference = lambda wildcards: config['reference'][wildcards.reference]
    output:
        'publicSamples/{sample}.{reference}.cov'
    resources:
        walltime = '15m'
    shell:
        '''
        #19664000-19667000
        samtools depth -aa -r 29:17990000-18000000 --reference {input.reference} {input.cram[0]} | awk -v S={wildcards.sample} '{{print S,$2,$3}}' > {output}
        '''

rule collate_depths:
    input:
        expand(rules.samtools_depth.output,sample=samples,reference=('ARS',))
    output:
        'publicSamples/coverage.csv.gz'
    localrule: True
    shell:
        '''
        {{ echo "sample position depth" ; cat {input} ; }} | pigz -c -p 2 > {output}
        '''
