import sys
import traceback

import SNV
import datetime
from typing import List

class VCFwriter:

    def __init__(self, args):
        self.filename = f"{args.output_folder}/output.vcf"
        self.__write_header(args)

    def __write_header(self, args):
        try:
            with open (self.filename, 'w') as file:
                file.write("##fileformat=VCFv4.2\n")
                file.write("##source=VariantCaller\n")
                file.write(f"##fileDate={datetime.datetime.now().strftime('%Y%m%d')}\n")
                file.write(f"##thresholds: min_depth={args.min_depth}, min_base_qual={args.min_base_qual}, min_alt_count={args.min_alt_count}, min_alt_freq={args.min_alt_freq}\n")
                file.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
                file.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
                file.write("##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand Bias\">\n")
                header = "{:<6}\t{:<10}\t{:<3}\t{:<3}\t{:<3}\t{:<4}\t{:<6}\t{}\t{}\t{}\n".format(
                    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"
                )
                file.write(header)
        except:
            print("Error writing VCF header")
            sys.exit(1)

    def __generate_vcf_entry(self, snv: SNV):
        return f"{snv.chrom:<{6}}\t{snv.pos:<{10}}\t{'.':<{3}}\t{snv.ref:<{3}}\t{snv.alt:<{3}}\t{snv.qual:<{4}}\t{snv.filter:<{6}}\t{snv.info()}\n"

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
