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
        self.alt_allel_frequency = 0.0
        self.qual = 0
        self.filter = None
        self.infoDic = {
            'DP': 0,
            'AF': 0.0
        }

    def info(self):
        if self.infoDic is None:
            return "NO SNV DETECTED"

        return ";".join(f"{key}={value}" for key, value in self.infoDic.items())

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
        while q != len(self.qualities):
            quality = self.qualities[q]
            base = self.bases[b]

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
                    base = alt = self.bases.upper()[b: b + 2 + int(size_of_indel)]
                    # set the next index for the base (since we have to skip the indel "information" in the bases string)
                    b = b + 1 + int(size_of_indel)

                    if alt in alternatives:
                        alternatives[alt].count += 1
                        alternatives[alt].quality_sum += phred
                    else:
                        alternatives[alt] = Alt(alt, 1)
                        alternatives[alt].quality_sum = phred

                elif base.upper() in 'ACTGNactgn' and base.upper() != self.ref:
                    alternatives[base.upper()].count += 1
                    alternatives[base.upper()].quality_sum += phred

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
        for alt in alternatives.values():
            if alt.count > highest_count:
                highest_count = self.alt_count = alt.count
                self.qual = alt.quality_sum

                if alt.name[0] == '+':
                    self.alt = self.ref + alt.name[2:]
                elif alt.name[0] == '-':
                    self.alt = self.ref
                    self.ref = self.ref + alt.name[2:]
                else:
                    self.alt = alt.name


        # calculate alternative allele frequency
        self.alt_allel_frequency = self.alt_count / self.depth

        self.infoDic['DP'] = self.depth
        self.infoDic['AF'] = self.alt_allel_frequency

        # since base quality, depth and altternative allele count are already checked, we only need to check for allele frequency
        return (self.alt_allel_frequency >= min_alt_frequency
                and self.alt_count > min_alt_count
                and self.depth > min_depth)

    def __calculate_phred(self, char):
        return ord(char) - 33

class Alt:

    def __init__(self, name="None", count=0, quality_sum=0):
        self.name = name
        self.count = count
        self.quality_sum = quality_sum
