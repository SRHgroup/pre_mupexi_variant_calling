#!/bin/python3

import argparse

def move_AF_to_FORMAT_column(input_vcf, output_vcf):
    with open(input_vcf, 'r') as invcf, open(output_vcf, 'w') as outvcf:
        for line in invcf:
            if line.startswith("#"):
                outvcf.write(line)
                continue
            if line.strip().split("\t")[9].split(":")[-2] == "SOMATIC":
                outvcf.write(line)
                continue
            if line.strip().split("\t")[9].split(":")[-2] == "GERMLINE":
                line = line.split("\t")
                AF = float(line[7].split(";")[1].split("=")[1])
                line[8] = ":".join(line[8].split(":")[:2] + ["AF"] + line[8].split(":")[2:])
                line[9] = ":".join(line[9].split(":")[:2] + [str(AF)] + line[9].split(":")[2:])
                outvcf.write("\t".join(line))
                continue


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Copy AF value and tag to samples and format columns")
    parser.add_argument("-i", "--input_vcf", required=True, help="Path to the input vcf file.")
    parser.add_argument("-o", "--output_vcf", required=True, help="Path to the output vcf file.")

    args = parser.parse_args()

    move_AF_to_FORMAT_column(args.input_vcf, args.output_vcf)
