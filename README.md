
For applications, we provide .Rmd files to conduct analysis on synthetic data. For simulation, we provide scripts and wrappers to reproduce results shown in the manuscript. We have developed an R package `multipleDL` for the method proposed. 

# Applications 

The data will not be made publicly available for ethical, privacy reasons. The dataset contains sensitive information. Furthermore, we do not own the data, and we made an agreement with those who collected the data that we will not share it. We are happy to provide others with the relevant information for how they might be able to obtain the data under reasonable requests.

We provide synthetic data here for reproducibility. 

Access to the original data in the multi-detection limit study requires submitting a concept sheet to the Caribbean, Central, and South America network for HIV Epidemiology (CCASAnet). The concept sheet is then assigned to a relevant working group that approves your concept sheet and then sends it to the CCASAnet Executive Committee who give the final approval. Each study site also has an opt-out option. Information on this process and the concept sheet are available at ccasanet.org.

Access to the original data used in the single detection limit study requires approval from the principal investigator of that study, Dr. John Koethe. 

Dr. Shepherd (bryan.shepherd@vanderbilt.edu) is happy with data requests for either of these studies.


## Multiple Detection Limits 

The data is from a multi-center HIV study in Latin America. The data include adults living with HIV starting antiretroviral therapy at one of 5 study centers in Latin America. Subjects in the study have HIV viral load observations fell below varying detection limits.

- Code: applications_multiple.Rmd

- Data: synthetic_data_multiple.Rda

- Data dictionary: data_dictionary.txt

## Single Detection Limits

- Code: applications_single.Rmd

- Data: synthetic_data_single.Rda

- Data dictionary: data_dictionary.txt

# Simulation

## Multiple Detection Limits

- Overall: multiple_dl_wrapper.R

- Table 2: multiple_scenario1.R, multiple_scenario2.R, multiple_scenario3.R, multiple_scenario4.R, multiple_scenario5.R

- Table 3: multiple_correct_incorrect_transformation.R

- Table 4: comparison_cai_shen.R (Please note the code is very slow. We set the number of replication to 10 for this one file.)

- Figure S5: multiple_p_values.R

- Table S7: multiple_mis_link.R

- Table S8: multiple_left_out.R

## Single Detetion Limits

- Overall: single_dl_wrapper.R

- Table S10: single_dl_scenario1.R, single_dl_scenario2.R, single_dl_scenario3.R, single_dl_scenario4.R, single_dl_scenario5.R, single_dl_scenario6.R