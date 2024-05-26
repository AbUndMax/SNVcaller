import requests


class VariantAnnotator:
    def __init__(self, server="https://rest.ensembl.org"):
        self.server = server

    def get_variant_annotation(self, chrom, pos, ref, alt):
        chrom = chrom.replace("chr", "")
        ext = f"/vep/human/region/{chrom}:{pos}-{pos}/{alt}"
        headers = {"Content-Type": "application/json"}
        try:
            response = requests.get(self.server + ext, headers=headers)
            response.raise_for_status()
        except requests.exceptions.HTTPError:
            return None
        if not response.ok:
            response.raise_for_status()
            return None
        return response.json()

    def parse_annotation(self, annotation):
        if "transcript_consequences" in annotation:
            for consequence in annotation["transcript_consequences"]:
                if "amino_acids" in consequence and "codons" in consequence:
                    return True
        return False

    def annotate_variant(self, chrom, pos, ref, alt):
        annotation = self.get_variant_annotation(chrom, pos, ref, alt)
        if annotation:
            for result in annotation:
                if self.parse_annotation(result):
                    return "Protein change"
        return "No protein change"