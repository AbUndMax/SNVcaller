import sys
import traceback
import SNV
import datetime
from typing import List


class VCFwriter:

    chrom_space = 10
    pos_space = 15
    id_space = 5
    ref_space = 5
    alt_space = 5
    qual_space = 6
    filter_space = 10
    info_space = 70
    format_space = 10
    sample_space = 0

    def __init__(self, args):
        self.filename = f"{args.output_folder}/output.vcf"
        self.__write_header(args)

    def __write_header(self, args):
        try:
            with open (self.filename, 'w') as file:
                file.write("##fileformat=VCFv4.2\n")
                file.write("##source=variant_caller.py by Niklas G.\n")
                file.write(f"##fileDate={datetime.datetime.now().strftime('%Y%m%d')}\n")
                file.write(f"##thresholds: min_depth={args.min_depth}, min_base_qual={args.min_base_qual}, min_alt_count={args.min_alt_count}, min_alt_freq={args.min_alt_freq}\n")
                file.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
                file.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
                file.write("##INFO=<ID=SB,Number=1,Type=Float,Description=\"Fisher Strand Bias\">\n")
                file.write("##INFO=<ID=FC, Number=1, Type=String, Description=\"Annotation of the Functional Consequence (via VEP)\">\n")
                file.write("##INFO=<ID=IMP, Number=1, Type=String, Description=\"Annotation of the Impact Level\">\n")
                file.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                header = ("#CHROM".ljust(self.chrom_space) + "POS".ljust(self.pos_space) + "ID".ljust(self.id_space) +
                "REF".ljust(self.ref_space) + "ALT".ljust(self.alt_space) + "QUAL".ljust(self.qual_space) +
                "FILTER".ljust(self.filter_space) + "INFO".ljust(self.info_space) + "\tFORMAT".ljust(self.format_space) +
                "SAMPLE".ljust(self.sample_space) + "\n")
                file.write(header)
        except:
            print("Error writing VCF header")
            sys.exit(1)

    def __generate_vcf_entry(self, snv: SNV):
        return (snv.chrom.ljust(self.chrom_space) + str(snv.pos).ljust(self.pos_space) + ".".ljust(self.id_space) +
                snv.ref.ljust(self.ref_space) + snv.alt.ljust(self.alt_space) + str(snv.qual).ljust(self.qual_space) +
                snv.filter.ljust(self.filter_space) + snv.info().ljust(self.info_space) +
                "\tGT".ljust(self.format_space) + snv.genotype() + "\n")

    def write_all_entries(self, snv_list: List[SNV]):
        try:
            with open(self.filename, 'a') as file:
                for snv in snv_list:
                    vcf_entry = self.__generate_vcf_entry(snv)
                    file.write(vcf_entry)
        except:
            print("Error writing VCF entry")
            traceback.print_exc()
            sys.exit(1)

    def get_vcf_file_name(self):
        return self.filename
