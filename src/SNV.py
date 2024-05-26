import scipy.stats as stats
import VariantAnnotator


class SNV:

    def __init__(self, chrom, pos, ref, depth, bases, qualities):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.depth = int(depth)
        self.bases = bases.strip()
        self.qualities = qualities.strip()

        # these fields will be filled by detect_snv method
        self.alt = None
        self.alt_count = 0
        self.alt_allele_frequency = 0.0
        self.qual = 0
        self.filter = None
        self.infoDic = {}

    def info(self):
        if not self.infoDic:
            return "."
        return ";".join(f"{key}={value}" for key, value in self.infoDic.items())

    def genotype(self):
        if self.alt_count == 0:
            return "0/0"  # no variant detected
        elif 0.3 <= self.alt_allele_frequency <= 0.7:
            return "0/1"  # heterozygous
        elif self.alt_allele_frequency > 0.9:
            return "1/1"  # homozygous
        else:
            return "unsure"


    def detect_snv(self, min_depth, min_alt_count, min_alt_frequency, min_base_quality):

        # initial depth check if coverage is sufficient
        if self.depth < min_depth:
            return False

        alternatives = {
            'A': Alt('A'),
            'C': Alt('C'),
            'G': Alt('G'),
            'T': Alt('T')
        }

        # count the number of alternative and reference alleles by using only high quality bases
        bases_above_threshold = []
        qualities_above_threshold = []
        depth = 0

        q = 0
        b = 0
        while q != len(self.qualities) and b != len(self.bases):
            quality = self.qualities[q]
            base = self.bases[b]
            base_in_upper = base.upper()

            # skip base if it is marking new read
            if base == '^':
                # synch the b iterator with the q iterator by skipping the new Read quality after ^
                b += 2
                q += 1
                continue

            phred = self.__calculate_phred(quality)
            if phred >= min_base_quality:

                if base in '+-':
                    size_of_indel = self.bases[b + 1]
                    base = self.bases[b: b + 2 + int(size_of_indel)]
                    base_in_upper = base.upper()
                    # set the next index for the base (since we have to skip the indel "information" in the bases string)
                    b = b + 1 + int(size_of_indel)

                    if base_in_upper in alternatives:
                        alternatives[base_in_upper].count += 1
                        alternatives[base_in_upper].quality_sum += phred
                    else:
                        alternatives[base_in_upper] = Alt(base_in_upper, 1, phred)

                elif base_in_upper in 'ACTGN' and base_in_upper != self.ref:
                    alternatives[base_in_upper].count += 1
                    alternatives[base_in_upper].quality_sum += phred

                if base.isupper():
                    alternatives[base_in_upper].forward_counter += 1
                elif base.islower():
                    alternatives[base_in_upper].reverse_counter += 1


                bases_above_threshold.append(base)
                qualities_above_threshold.append(quality)
                depth += 1

            b += 1
            q += 1

        self.bases = ''.join(bases_above_threshold)
        self.qualities = ''.join(qualities_above_threshold)
        self.depth = depth

        # calculate alternative base and alternative allele count
        highest_count = 0
        most_abundant_alt = None
        for alt in alternatives.values():
            if alt.count > highest_count:
                most_abundant_alt = alt
                highest_count = self.alt_count = alt.count
                self.qual = alt.quality_sum

                if alt.name[0] == '+':
                    self.alt = self.ref + alt.name[2:]
                elif alt.name[0] == '-':
                    self.alt = self.ref
                    self.ref = self.ref + alt.name[2:]
                else:
                    self.alt = alt.name

        # if there is no alternative allele, return False
        if most_abundant_alt is None:
            return False

        # calculate alternative allele frequency
        self.alt_allele_frequency = self.alt_count / self.depth

        # calculate strand bias
        ref_rev = self.bases.count(',')
        ref_fwd = self.bases.count('.')
        alt_rev = most_abundant_alt.reverse_counter
        alt_fwd = most_abundant_alt.forward_counter

        _, pvalue = stats.fisher_exact([[ref_fwd, ref_rev], [alt_fwd, alt_rev]])

        self.infoDic['DP'] = self.depth
        self.infoDic['AF'] = self.alt_allele_frequency
        self.infoDic['SB'] = pvalue

        # since base quality, depth and altternative allele count are already checked, we only need to check for allele frequency
        return (self.alt_allele_frequency >= min_alt_frequency
                and self.alt_count > min_alt_count
                and self.depth > min_depth)

    def __calculate_phred(self, char):
        return ord(char) - 33

    def annotate(self):
        annotator = VariantAnnotator.VariantAnnotator()
        annotation = annotator.annotate_variant(self.chrom, self.pos, self.ref, self.alt)

        if annotation is None:
            self.infoDic["Annotation"] = "No annotation found"
        else:
            self.infoDic["Annotation"] = annotation

class Alt:

    def __init__(self, name="None", count=0, quality_sum=0):
        self.name = name
        self.count = count
        self.quality_sum = quality_sum
        self.reverse_counter = 0
        self.forward_counter = 0
