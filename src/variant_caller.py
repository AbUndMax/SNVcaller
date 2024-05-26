import argparse
from VCFwriter import VCFwriter
from SNV import SNV


def parse_pileup_file(args):
    with open(args.pileup_file, 'r') as file:
        found_snv = []

        for line in file:
            line.strip()
            chrom, pos, ref, depth, bases, qualities = line.split('\t')
            maybe_snv = SNV(chrom, pos, ref, depth, bases, qualities)
            is_snv = maybe_snv.detect_snv(args.min_depth, args.min_alt_count, args.min_alt_freq, args.min_base_qual)
            if is_snv:
                maybe_snv.filter = 'PASS'
                found_snv.append(maybe_snv)

    return found_snv


def main():
    parser = argparse.ArgumentParser(description="Variant Caller. Reads Pileup file "
                                                 "and calls SNV's into VCF file.")
    parser.add_argument("-pf", "--pileup_file", required=True, help="Pileup file to read")
    parser.add_argument("-o", "--output_folder", help="Output folder to write VCF file")
    parser.add_argument("-md", "--min_depth", type=int, default=10, help="Minimum depth to call a variant")
    parser.add_argument("-mbq", "--min_base_qual", type=int, default=20, help="Minimum quality to call a variant")
    parser.add_argument("-mac", "--min_alt_count", type=float, default=4, help="Minimum count for alternative base allele")
    parser.add_argument("-maf", "--min_alt_freq", type=float, default=0.2, help="Minimum alternative frequency to call a variant")
    args = parser.parse_args()

    print("\n>>> Variant Caller running...\n")
    print(">>> set parameters: ")
    print("> minimum depth set:", args.min_depth)
    print("> minimum base quality set:", args.min_base_qual)
    print("> minimum alternative allele count set:", args.min_alt_count)
    print("> minimum alternative frequency set:", args.min_alt_freq, "\n")

    found_snv = parse_pileup_file(args)
    number_of_snv = 0
    number_of_insertions = 0
    number_of_deletions = 0

    for snv in found_snv:
        if len(snv.ref) > 1:
            number_of_deletions += 1
        elif len(snv.alt) > 1:
            number_of_insertions += 1
        else:
            number_of_snv += 1

    print(">>> Found SNV's: \t", len(found_snv))
    if number_of_insertions > 0:
        print(">>> Found Insertions: \t", number_of_insertions)
    if number_of_deletions > 0:
        print(">>> Found Deletions: \t", number_of_deletions)
    print()

    vcf_writer = VCFwriter(args)
    vcf_writer.write_all_entries(found_snv)
    print(">>> VCF file written to:", vcf_writer.get_vcf_file_name(), "\n")

    print(">>> Variant Caller finished.")


if __name__ == "__main__":
    main()
