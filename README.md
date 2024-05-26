# SNVcaller
A pipeline for calling SNVs from RNA-seq data

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
```
