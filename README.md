Do file 1_Preparing longitudinal datafile:   
      -> Prepares the longitudinal dataset for analysis, using 'poliovaccinestudy_Tanzania_aim1_baselinedata_092425.dta'; 'Polio vaccine study follow up dates.excl'; 'stool_serumneutralization_09242025.dta'  
      -> Applies WHO antropometic analysis package to generate wasting, stunting, overweight and underweight variables  
      -> Merges dates of follow-up, neutralisation data, and binding data
      -> Saves in long format with one row per visit
      -> Cleans data for tables (e.g. groups co-infections, calculates number of observations, tags missing data or age outliers etc)  
      -> Description of sample for table 1 and supplementary descriptive tables, heatmap  
      -> Output = 'merged_dataset_long.dta'  
  
Do file 2_Preparing cross sectional datafile:  
      -> Prepares the cross-sectional dataset for analysis, using 'poliovaccinestudy_cleandataset_aim2_10_20_2025.dta'  
      -> Applies WHO antropometic analysis package and cleans data for tables in the same way as above 
      -> Merges binding data
      -> Description of sample for table 1  
      -> Output = 'cross-sectional_clean.dta'  
  
Do file 3_Merging cross-sectional and longitudinal  
      -> Uses 'cross-sectional_clean.dta' and 'merged_dataset_long.dta'  
      -> Appends cross-sectional and longitudinal files in long format  
      -> Generates a unified visit variable (1=baseline, 2=1 month, 3=3 months, 4=6 months, 5=1 year, 6=2 years, 7=3 years, 8= 8-9 years)  
      -> Generates summary outcome variables for the regression analyses (any positive response between months 1-6, maximum titre between 1-6)  
      -> Estimates half-lives for neutralisation titres  
      -> Regressions with co-infections  
      -> Prepares correlations dataset by keeping only cross-sectional or month 1 for the longitudinal 
      -> Prepares more managable dataset for boxplots in long format (one row per assay performed)
      -> Prepares infections dataset by keeping only visits at month 1, without age outliers or missing infection status data and reshaping as long according to assay (PV1/3; stool/serum)
      -> Output = 'Both_datsets_long.dta'; 'correlations_dataset.dta'; 'infections_dataset.dta', regression outputs in excel saved as 'reg_long_results.xlsx'  
