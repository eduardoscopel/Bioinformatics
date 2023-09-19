from cyvcf2 import VCF

chr = [0]*17
i=0

for variant in VCF("ZP862.vcf"):
    if i == 0:
        chr[i] = variant.CHROM
        i+=1
    if chr[i-1] != variant.CHROM:
        chr[i] = variant.CHROM
        i+=1





    variant.CHROM, variant.POS, variant.ALT, variant.INFO.get('AF1')
