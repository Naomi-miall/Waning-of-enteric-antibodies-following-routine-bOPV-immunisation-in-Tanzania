Do file 1_Preparing longitudinal datafile:   
      -> Prepares the longitudinal dataset for analysis, using 'poliovaccinestudy_Tanzania_aim1_baselinedata_092425.dta'; 'Polio vaccine study follow up dates.excl'; 'stool_serumneutralization_09242025.dta'  
      -> Applies WHO antropometic analysis package to generate wasting, stunting, overweight and underweight variables  
      -> Merges dates of follow-up, neutralisation data, and binding data
      -> Saves in long format with one row per visit
      -> Cleans data for tables (e.g. groups co-infections, calculates number of observations, tags missing data or age outliers etc)  
      -> Output = 'merged_dataset_long.dta'  
  
Do file 2_Preparing cross sectional datafile:  
      -> Prepares the cross-sectional dataset for analysis, using 'poliovaccinestudy_cleandataset_aim2_10_20_2025.dta'  
      -> Applies WHO antropometic analysis package and cleans data for tables in the same way as above 
      -> Merges binding data
      -> Output = 'cross-sectional_clean.dta'  
  
Do file 3_Merging cross-sectional and longitudinal  
      -> Uses 'cross-sectional_clean.dta' and 'merged_dataset_long.dta'  
      -> Appends cross-sectional and longitudinal files in long format  
      -> Generates a unified visit variable (1=baseline, 2=1 month, 3=3 months, 4=6 months, 5=1 year, 6=2 years, 7=3 years, 8= 8-9 years)  

Do file 4_Main tables and figures
      -> Table 1 descriptives
      -> Prepares more managable dataset for boxplots in long format (one row per assay performed)
      -> Forest plot:
            -> Generates summary outcome variables for the regression analyses (any positive response between months 1-6, maximum titre between 1-6) 
            -> Runs regressions for each covariate to include in the forest plot
      
Do file 5_Supplementary materials  
      -> Flow diagram
      -> Median and IQR for outcomes plus p-values for supplementary tables 2 and 4
      -> Estimates half-lives for neutralisation titres  
      -> Comparing the full longitudinal dataset to the complete data sample
      -> Prepares correlations dataset by keeping only cross-sectional or month 1 for the longitudinal
      -> Type 2 analysis

Do file 6_Total IgA analysis
      -> Estimates total IgA stratified by timepoint and study group


R file_Fig 2 Box plots script.R
      -> Uses Box_plots.dta to create a box plot for neut, IgA and IgG titers for PV1 and PV3 (main text - saved as Figure 1.png) and PV2 (supplement - saved as Figure 1.2_PV2.png) in stool and serum, stratified by timepoint and study group.


