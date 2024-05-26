class SNV:

    def __init__(self, chrom, pos, ref, depth, bases, qualities):
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.depth = int(depth)
        self.bases = bases
        self.qualities = qualities

        # these fields will be filled by detect_snv method
        self.alt = None
        self.alt_count = 0
        self.alt_allel_frequency = 0.0
        self.qual = 0
        self.filter = None
        self.infoDic = None

    def info(self):
        if self.infoDic is None:
            return "NO SNV DETECTED"

        return ";".join(f"{key}={value}" for key, value in self.infoDic.items())

    def detect_snv(self, min_depth, min_alt_count, min_alt_frequency, min_base_quality):

        # initial depth check if coverage is sufficient
        if self.depth < min_depth:
            return False

        base_count = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }
        alt_phred_sum = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0
        }

        # count the number of alternative and reference alleles by using only high quality bases
        bases_above_threshold = []
        qualities_above_threshold = []
        depth = 0

        for base, quality in zip(self.bases, self.qualities):

            # check if base quality is above threshold
            phred = self.__calculate_phred(quality)

            if phred >= min_base_quality:
                possible_alts = ['A', 'C', 'G', 'T']
                possible_alts.remove(self.ref.upper())
                if base.upper() in possible_alts:
                    base_count[base.upper()] += 1
                    alt_phred_sum[base.upper()] += phred

                bases_above_threshold.append(base)
                qualities_above_threshold.append(quality)
                depth += 1

        self.bases = ''.join(bases_above_threshold)
        self.qualities = ''.join(qualities_above_threshold)
        self.depth = depth

        # calculate alternative base and alternative allele count
        high_score = 0
        for key, value in base_count.items():
            if value > high_score:
                high_score = value
                self.alt = key
                self.alt_count = value
                self.qual = alt_phred_sum[key]

        # calculate alternative allele frequency
        self.alt_allel_frequency = self.alt_count / self.depth

        self.infoDic = {
            'DP': self.depth,
            'AF': self.alt_allel_frequency
        }

        # since base quality, depth and altternative allele count are already checked, we only need to check for allele frequency
        return (self.alt_allel_frequency >= min_alt_frequency
                and self.alt_count > min_alt_count
                and self.depth > min_depth)

    def __calculate_phred(self, char):
        return ord(char) - 33
