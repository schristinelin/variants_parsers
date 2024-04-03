```
Usage: variant_parsers.py [OPTIONS] GENE_NAME

Options:
  -of, --regen_output_files       Set to true to regenerate output files
  -am, --regen_alphamissense_data
                                  Set to true to regenerate alphamissense data
  -h, --help                      Show this message and exit.

```

Example:
``` python3 variant_parsers.py BRCA1 test_merged_output.csv ``` 


### References:

#### BRCA1
1. Functional dataset: Supplementary Table 1 - Findlay, G.M., Daza, R.M., Martin, B. et al. Accurate classification of BRCA1 variants with saturation genome editing. Nature 562, 217–222 (2018). https://doi.org/10.1038/s41586-018-0461-z
2. Predictor dataset:
   a. Clinvar - https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1%5Bgene%5D&redir=gene

#### MSH2
1. Functioal dataset: Data S1. Tables S5 - Xiaoyan Jia, Bala Bharathi Burugula, Victor Chen, Rosemary M. Lemons, Sajini Jayakody, Mariam Maksutova, Jacob O. Kitzman, Massively parallel functional testing of MSH2 missense variants conferring Lynch syndrome risk, The American Journal of Human Genetics, Volume 108, Issue 1, 2021, Pages 163-175, ISSN 0002-9297, https://doi.org/10.1016/j.ajhg.2020.12.003.
2. Predictor dataset:
   a. Clinvar - https://www.ncbi.nlm.nih.gov/clinvar/?term=MSH2%5Bgene%5D
