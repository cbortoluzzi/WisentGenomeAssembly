## genome-wide analysis
echo "chromosome total Bison_bison Bison_bonasus Bos_gaurus Bos_grunniens Bos_indicus Bos_mutus"
for C in {1..29}
do
  echo -n "${C} $(bcftools view ${C}.vcf.gz | bcftools stats | awk '$1=="SN"&&$5=="records:" {printf $6" "}')"
  for S in Bison_bison Bison_bonasus Bos_gaurus Bos_grunniens Bos_indicus Bos_mutus
  do
    bcftools view -s ${S} -x ${C}.vcf.gz | bcftools stats | awk '$1=="SN"&&$5=="records:" {printf $6" "}'
  done
  echo 
done

## specific region analysis for wisent
bedtools makewindows -g <(head -n 29 ARS-UCD1.2_Btau5.0.1Y.fa.fai) -w 100000 > ARS.100kb.bed

bcftools query -f '%CHROM\t%POS\t%POS' autosomes.vcf.gz > autosomes.all.bed
bedtools coverage -a ARS.100kb.bed -b autosomes.all.bed > intersect.all.bed

for S in Bison_bison Bison_bonasus Bos_gaurus Bos_grunniens Bos_indicus Bos_mutus
do
  bcftools view -s $S -x autosomes.vcf.gz | bcftools query -f '%CHROM\t%POS\t%POS' > autosomes.${S}.bed
  bedtools coverage -a ARS.100kb.bed -b autosomes.${S}.bed > intersect.${S}.bed

  echo "${S}: "
  paste <(cut -f -4 intersect.all.bed) <(cut -f 5 intersect.${S}.bed) |  awk -v OFS='\t' '($4-$5)<=-10 {print $1,$2,$3,$4-$5}' | sort -k1,1n -k2,2n |  bedtools merge -d 100000 -i - | awk '($3-$2)>=300000'
  echo ""
done

