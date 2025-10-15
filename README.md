# Acoustic Detections Increase in Citizen Science Data with the use of Automated Detection Apps
This repository is for the submission of the manuscript titled above. Note that this repository contains **code only** aside from the calculated audio visual index. Information about each related dataset can be found in the manuscript but a short blurb is given below about data access.

### Merlin Sound Identification Data (MSID)
MSID data are not freely available and requests to use such data should be directed to the Cornell Lab of Ornithology whom own the dataset.

### eBird Basic Dataset (EBD)
We analyse the March 2025 release of the EBD (Cornell Lab of Ornithology, 2025). The related sampling file from this dataset is also used to process the number of checklists each observer has recorded prior to the study period start data of March 2023. To access the full eBird Basic Dataset required to run the initital data cleaning script please see here https://ebird.org/data/download.

### eBird Status and Trends
To process species occurrence across the study region we use occurrence estimates from eBird Status and Trends. These are not available in a data repository as this can be downloaded from inside R with a relevant eBird access token which can be found here https://ebird.org/st/request.

### Macaulay Library Records
The audio visual index used in this dataset is processed from records from the online Macaulay Library. We publish the final audio visual index used in the analysis but not the raw data associated with Macaulay Library records. This file can be found in folder `AVI` below.

## Information about code
The following is a breakdown of scripts for the analysis. All scripts contain a description of their purpose in the header too.

| Script                   | Purpose                                                                                                   |
|--------------------------|-----------------------------------------------------------------------------------------------------------|
| 00 eBird pre-processing  | Process the full March 2025 eBird Basic Dataset                                                           |
| 01 eBird processing      | Filter processed script to be only relevant observers, years and months                                   |
| 02 OH processing         | Calculates observer histories from the processed eBird sampling dataset                                   |
| 03 MSID processing       | Calculates observers proportional use of MSID                                                             |
| 04 ML processing         | Calculates the audio visual index from the raw Macualay Library data                                      |
| 05 SA processing         | Calculates species availability from eBird Status product                                                 |
| 06 RR processing         | Calculates observer level species reporting rates and aggregates all previous covariates                  |
| 07 data anonymisation    | Anonymises the observers for all relevant analysis scripts                                                |
| 08 individual models     | Runs models for each individual species                                                                   |
| 09a community models     | Runs one model for all species in individual models                                                       |
| 09b community models max | Runs one model for all species, including those not in individual models                                  |
| 10 model comparison      | Compares species results for individual models and community model                                        |
| 11 phylogenetic trees    | Plots phylogenetic tree plots for individual, community and community max species results                 |
| 12 phylogenetic analysis | Regression with respect to phylogenetic distance on species traits for individual models' species results |

## References
Cornell Lab of Ornithology. (2025). eBird Basic Dataset (Version EBD_relMar-2025) [Data set]
