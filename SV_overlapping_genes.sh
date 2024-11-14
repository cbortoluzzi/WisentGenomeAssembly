for T in muscle liver
do
  for C in extreme normal
  do
    ## get coordinates from ID
    zgrep -wf <(cut -f 1 ${T}_genes_${C}.csv) GCF_002263795.3_ARS-UCD2.0_genomic.gtf.gz |\
    awk -v OFS='\t' '$3=="exon" {print $1,$4,$5,$10}' |\
    bedtools slop -i /dev/stdin -g GCA_002263795.4_ARS-UCD2.0_genomic.hard_masked_PAR_X.fa.gz.fai -b 0 |\
    sort -k1,1V -k2,2n -u > ${T}_genes_${C}.bed
    
    ##check for overlaps
    echo -en "Checking for overlaps in ${T}_${C}\n$(bedtools intersect -a ${T}_genes_${C}.bed -b TPM_SV/non_cattle_SVs.1Kb.vcf.gz | cut -f 4 | sort -u | wc -l) out of $(cat ${T}_genes_${C}.bed | cut -f 4 | sort -u | wc -l)\n"

    echo -en "Checking for overlaps in ${T}_${C}\n$(bedtools intersect -a ${T}_genes_${C}.bed -b TPM_SV/bison.1Kb.vcf.gz | cut -f 4 | sort -u | wc -l) out of $(cat ${T}_genes_${C}.bed | cut -f 4 | sort -u | wc -l)\n"

    echo -en "Checking for overlaps in ${T}_${C}\n$(bedtools intersect -a ${T}_genes_${C}.bed -b SVs/SV_bisons.vcf.gz | wc -l) out of $(cat ${T}_genes_${C}.bed | wc -l)\n"

    echo -en "Checking for overlaps in ${T}_${C}\n$(bedtools intersect -a ${T}_genes_${C}.bed -b SVs/SVs.vcf.gz | wc -l) out of $(cat ${T}_genes_${C}.bed | wc -l)\n"

  done
done
