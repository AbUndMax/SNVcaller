import requests


class VariantAnnotator:
    def __init__(self, server="https://grch37.rest.ensembl.org"):
        self.server = server

    def send_request(self, chrom, pos, ref, alt):

        if len(ref) > 1:
            # deletion
            alt = "-"
            end_pos = pos + len(ref) - 1

        elif len(alt) > 1:
            # insertion
            alt = alt[1:]
            pos = pos - 1

        else:
            end_pos = pos

        chrom = chrom.replace("chr", "")
        ext = f"/vep/human/region/{chrom}:{pos}-{end_pos}:1/{alt}"
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

    def _parse_consequence_annotation(self, response):
        annotation = response[0]

        consequence = "unknown"
        impact_level = "unknown"

        if "transcript_consequences" in annotation:

            for dic in annotation["transcript_consequences"]:
                if "impact" in dic:
                    impact_level = dic["impact"]
                break

        if "most_severe_consequence" in annotation:
            consequence = annotation["most_severe_consequence"]

        return consequence, impact_level

    def annotate_consequence(self, chrom, pos, ref, alt):
        response = self.send_request(chrom, pos, ref, alt)

        if response:
            consequence, impact = self._parse_consequence_annotation(response)

            return {
                "MSC": consequence,
                "IMP": impact
            }

        return {}
