samples = ['Bison_bison','Bison_bonasus','Bos_gaurus','Bos_grunniens','Bos_indicus','Bos_mutus','Bos_taurus']

rule all:
    input:
         expand('graphs/minigraph/{chromosome}.gfa',chromosome=range(1,30))

rule split_fasta:
    input:
        'graph_assemblies/{sample}.fa.gz'
    output:
        fa = 'assemblies/{chromosome}/{sample}.fa',
        fai = 'assemblies/{chromosome}/{sample}.fa.fai',
    shell:
        '''
        echo "{wildcards.chromosome}" | seqtk subseq {input} - |\
        sed 's/>.*/>{wildcards.sample}/' > {output.fa}
        samtools faidx {output.fa}
        '''

rule mash_triangle:
    input:
        lambda wildcards: expand('graph_assemblies/{sample}.fa.gz',sample=samples)
    output:
        'tree/lower_triangle.txt'
    threads: 6
    resources:
        mem_mb = 2500,
        walltime = '4h'
    shell:
        '''
        mash triangle -s 10000 -k 25 -p {threads} {input} | awk 'NR>1' > {output}
        '''

import numpy as np
from scipy.spatial.distance import squareform

def read_mash_triangle(mash_triangle,strip_leading=False):
    names, vals = [], []
    with open(mash_triangle,'r') as fin:
        for i,line in enumerate(fin):
            parts = line.rstrip().split()
            if strip_leading:
                names.append(parts[0].split('/')[1].split('.')[0])
            else:
                names.append(parts[0])
            vals.append(parts[1:]+[0])
    Q = np.asarray([np.pad(a, (0, len(vals) - len(a)), 'constant', constant_values=0) for a in vals],dtype=float)
    return names, (Q+Q.T)

def make_minigraph_order(mash_triangle,sequences=None):
    ref_ID = 'Bos_taurus'
    names, dists = read_mash_triangle(mash_triangle)
    sequences = sequences or names
    ref_ID_index = [i for i,s in enumerate(sequences) if ref_ID in s][0]

    return ' '.join(map(lambda k:k[0], sorted(zip(sequences,dists[ref_ID_index]),key=lambda k: k[1])))

rule minigraph_construct:
    input:
        mash_distances = rules.mash_triangle.output,
        assemblies = expand('assemblies/{{chromosome}}/{sample}.fa', sample=samples)
    output:
        temp('graphs/minigraph/{chromosome,\d+}.basic.gfa')
    threads: 1
    resources:
        mem_mb = 20000,
        walltime = '4h'
    params:
        sample_order = lambda wildcards, input:make_minigraph_order(input.mash_distances[0],input.assemblies),
        L = 50,
        j = 0.2
    shell:
        '''
        minigraph -t {threads} -cxggs -j {params.j} -L {params.L} {params.sample_order} > {output}
        '''

rule minigraph_call:
    input:
        gfa = rules.minigraph_construct.output,
        sample = 'assemblies/{chromosome}/{sample}.fa'
    output:
        'graphs/minigraph/{chromosome}.{sample}.bed'
    params:
        L = 50,
        j = 0.2
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        '''
        minigraph -t {threads} -cxasm --call -j {params.j} -L {params.L} {input.gfa} {input.sample} > {output}
        '''

rule minigraph_path:
    input:
        paths = expand('graphs/minigraph/{{chromosome}}.{sample}.bed',sample=filter(lambda x: x != 'Bos_taurus',samples)),
        gfa = 'graphs/minigraph/{chromosome}.basic.gfa'
    output:
        'graphs/minigraph/{chromosome}.gfa'
    params:
        samples = '\\n'.join(filter(lambda x: x != 'Bos_taurus',samples))
    localrule: True
    shell:
        '''
        #needs special branch of mgutils.js
        {{ vg convert -r 0 -g {input.gfa} -f ; paste {input.paths} | /cluster/work/pausch/alex/software/minigraph/misc/mg_old.js path <(echo -e "{params.samples}") - ; }} > {output}
        '''

rule vg_deconstruct:
    input:
        rules.minigraph_path.output
    output:
        vcf = '{chromosome}.vcf',
        SVs = 'SVs.{chromosome}.count'
    resources:
        mem_mb = 7500
    shell:
        '''
        vg deconstruct -p Bos_taurus -a -d 1 -e {input} |\
        bcftools norm -m -any |\
        bcftools view -i 'abs(ILEN)>=50&&abs(ILEN)<=100000' |\
        bcftools +setGT - -- -t . -n 0 |\
        tee {output.vcf} |\
        bcftools query -f '%REF %ALT [%GT]\\n' |\
        awk '{{length($1)>length($2)?++A["0"$3]:++B["0"$3]}} END {{ for (k in A) {{print "DEL",k,A[k]}}; for (k in B) {{print "INS",k,B[k]}}; }}' > {output.SVs}
        '''

rule count_SVs:
    input:
        expand(rules.vg_deconstruct.output['SVs'],chromosome=range(1,30))
    output:
        'SVs.count'
    localrule: True
    shell:
        '''
        awk '{{A[$1][$2]+=$3}} END {{ for (k in A) {{ for (V in A[k]) {{ print k,V,A[k][V] }} }} }}' {input} > {output}
        '''

rule vep:
    input:
        vcf = rules.vg_deconstruct.output['vcf'],
        gff = config['gff'],
        reference = config['reference']
    output:
        stats = 'vep/{chromosome}.vep.stats',
        vcf = 'vep/{chromosome}.vep.vcf.gz'
    threads: 1
    resources:
        mem_mb = 5000
    shell:
        '''
        perl -e 'use DBI;'
        vep --plugin SpliceRegion \
        -i {input.vcf} \
        -gff {input.gff} \
        -fasta {input.reference} \
        --hgvs --symbol -o stdout \
        --stats_file {output.stats}
        --vcf | bgzip -p {threads} -c > {output.vcf}
        '''

rule vep_filter:
    input:
        rules.vep.output['vcf']
    output:
        'vep.HIGH.vcf'
    localrule: True
    shell:
        '''
        zcat {input} | awk '$11==1&&$10==1&&($11 + $10 + $12 +$13 +$14 +$15+$16)==2&&$0~/HIGH/ {{print $1, $2, $3, $4, $5, $6, $7, $8, $9}}' > {output}
        #bcftools +split-vep -c IMPACT {input.vcf} | bcftools view -s Bison_bison,Bison_bonasus -x -i 'INFO/IMPACT=="HIGH"' | cut -f -9 > {output}
        '''
