#!/bin/python3

import argparse

def get_AF(input_AF_vcf, input_nonAF_vcf, output_vcf):
    ''' this function takes the AF values from the VCF pre-phased fies and adds them to the post phased files. also adds the AF tag in FORMAT'''


    with open(input_AF_vcf, 'r') as afvcf:
        AFdict={}
        for line in  afvcf:
            if line.startswith("#"):
                continue

            variant = line.strip().split("\t")
            chrom_pos = (variant[0], variant[1])
            AF = variant[9].split(":")[3]
            is_somatic = ":SOMATIC" in variant[9]
            if is_somatic:
                AFdict[chrom_pos] = AF


    with open(input_nonAF_vcf, 'r') as nonafvcf, open(output_vcf, 'w') as outvcf:
        for line in nonafvcf:
            if line.startswith("#"):
                outvcf.write(line)
                continue

            variant = line.strip().split("\t")
            chrom_pos = (variant[0], variant[1])
            is_somatic = ":SOMATIC:" in line.strip().split("\t")[9]
            is_phased = "|" in variant[9].split(":")[0]

            if is_somatic and is_phased:
                variant[8] = ":".join(variant[8].split(":")[:2] + ["AF"] + variant[8].split(":")[2:])
                variant[9] = ":".join(variant[9].split(":")[:2] + [AFdict[chrom_pos]] + variant[9].split(":")[2:])
                line = "\t".join(variant) + "\n"
                outvcf.write(line)

            else:
                outvcf.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add AF to a VCF based on somatic phased variants.")
    parser.add_argument("-a", "--input_AF_vcf", required=True, help="Path to the input AF VCF file")
    parser.add_argument("-n", "--input_nonAF_vcf", required=True, help="Path to the input non-AF VCF file")
    parser.add_argument("-o", "--output_vcf", required=True, help="Path to the output VCF file")

    args = parser.parse_args()

    get_AF(args.input_AF_vcf, args.input_nonAF_vcf, args.output_vcf)



