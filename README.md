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
   b. popEVE - https://pop.evemodel.org/protein/NP_009225-1

#### MSH2
1. Functional dataset: Data S1. Tables S5 - Xiaoyan Jia, Bala Bharathi Burugula, Victor Chen, Rosemary M. Lemons, Sajini Jayakody, Mariam Maksutova, Jacob O. Kitzman, Massively parallel functional testing of MSH2 missense variants conferring Lynch syndrome risk, The American Journal of Human Genetics, Volume 108, Issue 1, 2021, Pages 163-175, ISSN 0002-9297, https://doi.org/10.1016/j.ajhg.2020.12.003.
2. Predictor dataset:
   a. Clinvar - https://www.ncbi.nlm.nih.gov/clinvar/?term=MSH2%5Bgene%5D
   b. popEVE - https://pop.evemodel.org/protein/NP_000242-1

#### TP53
1. Functional datasets:

   a. Supplementary Table 1 - Giacomelli, A.O., Yang, X., Lintner, R.E. et al. Mutational processes shape the landscape of TP53 mutations in human cancer. Nat Genet 50, 1381–1387 (2018). https://doi.org/10.1038/s41588-018-0204-y
2. Predictor dataset:
   a. Clinvar - https://www.ncbi.nlm.nih.gov/clinvar/?term=TP53%5Bgene%5D&redir=gene
   b. popEVE - https://pop.evemodel.org/protein/NP_000537-3

#### PTEN
1. Functional datasets:

   a. Table S2 - Taylor L. Mighell, Sara Evans-Dutson, Brian J. O’Roak, A Saturation Mutagenesis Approach to Understanding PTEN Lipid Phosphatase Activity and Genotype-Phenotype Relationships, The American Journal of Human Genetics, Volume 102, Issue 5, 2018, Pages 943-955, ISSN 0002-9297, https://doi.org/10.1016/j.ajhg.2018.03.018 (https://www.sciencedirect.com/science/article/pii/S0002929718301071)
   b. Supplementary Dataset 1 - Matreyek, K.A., Starita, L.M., Stephany, J.J. et al. Multiplex assessment of protein variant abundance by massively parallel sequencing. Nat Genet 50, 874–882 (2018). https://doi.org/10.1038/s41588-018-0122-z
2. Predictor dataset:
   a. Clinvar - https://www.ncbi.nlm.nih.gov/clinvar/?term=PTEN%5Bgene%5D&redir=gene

#### AlphaMissense data for all genes:
Jun Cheng et al., Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science381, eadg7492(2023). DOI:10.1126/science.adg7492

Retrieved from https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?prefix=&forceOnObjectsSortingFiltering=false

### popEVE data for all genes: 
Deep generative modeling of the human proteome reveals over a hundred novel genes involved in rare genetic disorders. Rose Orenbuch, Aaron W. Kollasch, Hansen D. Spinner, Courtney A. Shearer, Thomas A. Hopf, Dinko Franceschi, Mafalda Dias, Jonathan Frazer, Debora S. Marks, medRxiv 2023.11.27.23299062; doi: https://doi.org/10.1101/2023.11.27.23299062

Retrieved from https://pop.evemodel.org/
