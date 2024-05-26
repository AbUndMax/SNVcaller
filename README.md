# SNVcaller
A batch program for calling filter SNVs from pileup files to Variant Call Format (VCF) files.

## Running
start the program by running the following command:
```
python3 snv_caller.py -pf <path to the pileup file> -o <output directory> 
```

#### Optional parameters
```
-md  --min_depth        <Minimum depth to call a variant | defult: 10>
-mbq --min_base_qual    <Minimum quality to call a variant | default: 20>
-mac --min_alt_count    <Minimum count for alternative base allele | default: 4>
-maf --min_alt_freq     <Minimum alternative frequency to call a variant | default: 0.2>
-a   --annotate         <Annotate the the SNV's with effect prediction information | default: False>
```

The annotation is done via [VEP REST API](https://rest.ensembl.org/documentation/info/vep_region_get) from the Ensemble 
genome browser. 
(Peter W Harrison et al. Ensembl 2024 Nucleic Acids Res. 2024, 52(D1):D891–D899 PMID: 37953337 [10.1093/nar/gkad1049](https://academic.oup.com/nar/article/52/D1/D891/7416379?login=false))
