#!/bin/python3

import argparse

def filter_germline_variants(input_vcf, output_vcf, distance):
    somatic_variants=[]
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            chrom, pos, gt, block = columns[0], int(columns[1]), columns[9].split(":")[0], columns[9].split(":")[-1]

            if columns[9].split(":")[-2] == "SOMATIC":
                somatic_variants.append([chrom, pos, gt, block])

        infile.seek(0)

        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            columns = line.strip().split("\t")
            chrom, pos, gt, block = columns[0], int(columns[1]), columns[9].split(":")[0], columns[9].split(":")[-1]

            if columns[9].split(":")[-2] == "SOMATIC":
                outfile.write(line)

            if columns[9].split(":")[-2] == "GERMLINE":
                for som_var in somatic_variants:
                    if som_var[0] == chrom and abs(som_var[1] - pos) <= distance:
                        if gt == "1/1" or gt == "1|1":
                            outfile.write(line)
                            break
                        elif "|" in gt and gt == som_var[2] and block == som_var[3]:
                            outfile.write(line)
                            break

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter germline variants from a VCF file based on proximity to somatic variants and phasing.")
    parser.add_argument('-i', '--input_vcf', required=True, help="Input VCF file")
    parser.add_argument('-o', '--output_vcf', required=True, help="Output VCF file")
    parser.add_argument('-d', '--distance', type=int, required=True, help="Maximum allowed distance between somatic and germline variants (in base pairs)")

    args = parser.parse_args()

    filter_germline_variants(args.input_vcf, args.output_vcf, args.distance)
