*Cleaning cross sectional dataset for merging
*to get pvx_stool; pvxneuttitersserum; no_obs; visit; date; femalesex; ageweek; age_mother; complete_edu_adult; minutes_towatersource; hiv_status_mother; preceding_birth_interval_months; infant_feeding_type; currently_breastfed; birthweight; anthro; 


**# Anthro measures
*Anthro measures:
/* Indicate to the Stata compiler where the igrowup_restricted.ado file is stored*/
adopath + "Stata code\WHO anthro macro/"

/* Load the data file */
use "poliovaccinestudy_cleandataset_aim2_10_20_2025.dta", clear
describe

list pid dob_child bopv14w_date interview_date if pid=="POL 292"
replace interview_date = 23450 if pid=="POL 292"
list pid dob_child bopv14w_date interview_date if pid=="POL 292"

/* Generate the first three parameters reflib, datalib & datalab	*/

	//reflib - specify the file where the five STATA data sets containing the WHO Child Growth Standards are stored
	gen reflib="Stata code\WHO anthro macro/"
	lab var reflib "Directory of reference tables"

	//datalib - specify the file where the input STATA data set is stored.
	gen datalib="Stata code\Alex's code/"
	lab var datalib "Directory for datafiles"

	//datalab - specify the name that will prefix the output files (datalab_z_rc and datalab_prev_rc).
	gen datalab="Anthrooutputs_crosssectional"
	lab var datalab "Working file"

//For all variables ensure any missing are coded to ==.
/*	check the variable for "sex"	1 = male, 2=female */
gen sex = .
replace sex = 1 if femalesex == 0 // Can I double check coding of childgender? femalesex variable makes it look like 0=male?
replace sex = 2 if femalesex == 1  
label define sexlbl 1 "Male" 2 "Female"
label values sex sexlbl
tab sex femalesex, m

/*	check the variable for "age" */
**need to convert age in days to compare norms 
	gen agedays = interview_date - dob_child
	label var agedays "Age in days at time of interview"
	summ agedays, d
	*checking:
	gen agemonth2 = agedays/30
    count if agemonth!=agemonth2
	list agemonth agemonth2 dob_child bopv14w_date interview_date in 1/10

	/*	define your ageunit	*/
	gen ageunit="days"
	lab var ageunit "=days or =months"

/*	check the variable for body "weight" which must be in kilograms*/
/* 	NOTE: if not available, please create as [gen weight=.]*/
	summ weight_kgs, d

/* 	check the variable for "height" which must be in centimeters	*/ 
/* 	NOTE: if not available, please create as [gen height=.]*/
	summ height_cm, d

/*	check the variable for "measure"*/
//specify the name of a character variable indicating whether recumbent length or standing height was measured. The values of this variable must be "l" or "L" for recumbent length, and "h" or "H" for standing height. Users must code any missing values as " "
/* 	NOTE: if not available, please create as [gen str1 measure=" "]*/
	gen measure = ""
	replace measure ="l" if height_measurement_method==1
	replace measure = "h" if height_measurement_method == 2
	tab measure height_measurement_method, m
	label var measure "Length (l) or height (h) measured?"

//The macro will use height, measure and age to create _clenhei - the adjusted length/height with 0.7 where relevant.

/* 	check the variable for "oedema"*/
//values of this variable must be "n" or "N" for non-oedema, and "y" or "Y" for oedema. Macro assumes those with missing values are non-oedema. For oedema cases, weight-related z-scores (_zwei, _zwfl and _zbmi) are not calculated.
/* 	NOTE: if not available, please create as [gen str1 oedema="n"]*/
	tab legedema, m
	tab legedema, nolab //slight oedema == 1; none == 0
	gen oedema = ""
	replace oedema = "n" if legedema==0
	replace oedema = "y" if legedema==1
	tab legedema oedema, m

/*	check the variable for "sw" for the sampling weight*/
/* 	NOTE: if not available, please create as [gen sw=1] */
	gen sw=1
	
/* 	Fill in the macro parameters to run the command */
igrowup_restricted reflib datalib datalab sex agedays ageunit weight_kgs height_cm measure oedema sw

use "Stata code\Alex's code\Anthrooutputs_crosssectional_z_rc.dta"

*check height adjustment worked
*recommended that lenght measured for children < 731 days, if height measured for children < 731 days then add 0.7 cm. 
	bysort measure: summ agedays
	sort pid
	tab pid if measure=="h" & agedays<731
	list pid agedays measure height_cm _clenhei if height_cm!=_clenhei // looks good

*checking z scores worked
list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl in 1/15 //looks good

/* exploring extremes
_fwei=1 if (_zwei < -6 | _zwei >5) - extreme weight for age
_flen=1 if (_zlen < -6 | _zlen >6) - extreme height for age
_fbmi=1 if (_zbmi < -5 | _zbmi >5) - extreme bmi (for age)
_fwfl=1 if (_zwfl < -5 | _zwfl >5) - extreme weight for height
*/

list pid agedays weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if legedema==1

	summ _fwei _flen _fbmi _fwfl
	
	summ height_cm, d //median amongst participants is 76.3cm
	summ weight_kgs, d //median amongst participants is 10kgs
	
list pid agemonth weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if _flen==1
	
list pid agemonth weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if _fbmi==1

list pid agemonth weight_kgs _clenhei legedema _cbmi _zwei _zlen _zbmi _zwfl if _fwfl==1

/*Generating binary variables
WHO Indicator Definitions: 
	Stunting Height-for-age < -2SD
	Severe wasting Weight for Height < -3 SD
	Wasting Weight for Height < -2 SD
	Overweight Weight for Height > +2 SD
	Underweight Weight-for-age < -2SD
 SOURCE:https://cdn.who.int/media/docs/default-source/child-growth/child-growth-standards/software/anthro-survey-analyser-quickguide.pdf?sfvrsn=dc7ddc6f_6
*/

gen stunting = (_zlen<-2) if !missing(_zlen)
	label var stunting "Stunting defined as height-for-age below -2SD"
	tab stunting
	list pid sex agemonth _clenhei _zlen if stunting==1

gen wasting = (_zwfl<-2) if !missing(_zwfl)
	label var wasting "Wasting defined as weight-for-height below -2SD"
	tab wasting
	list pid sex agedays height_cm weight_kgs _zwfl if wasting==1

gen severe_wasting = (_zwfl<-3) if !missing(_zwfl)
	label var severe_wasting "Severe wasting defined as weight-for-height below -3SD"
	tab severe_wasting
	list pid sex agedays height_cm weight_kgs _zwfl wasting if severe_wasting==1
	
gen underweight = (_zwei<-2) if !missing(_zwei)
	label var underweight "Underweight defined as weight-for-age below -2SD"
	tab underweight
	list pid sex agedays weight_kgs _zwei if underweight==1	

gen overweight = (_zwei>2) if !missing(_zwei)
	label var overweight "Overweight defined as weight-for-age above 2SD"
	tab overweight
	list pid sex agedays weight_kgs _zwei if overweight==1	
	
	drop sw measure oedema ageunit sex datalab datalib reflib _agedays _clenhei _cbmi _fwfl _flen _fwei _fbmi
	drop _zwei _zlen _zbmi
	
save "Stata code\Alex's code\Anthrooutputs_crosssectional_z_rc.dta", replace
//////////////////////////////////////////////////////////////////

**# Data cleaning for coinfections
	*Data cleaning: 
	*Creating co-infection exposure groupings and dropping repeat variables
use "Stata code\Alex's code\Anthrooutputs_crosssectional_z_rc.dta", clear

*HIV
tab hiv_pregnancy hiv_status_mother,m //230 missing on hiv_pregnancy are HIV negative. Of 16 monthers HIV positive, 12 were positive during pregnancy
gen hiv_exposure = 0
replace hiv_exposure = 1 if hiv_pregnancy==1
drop hiv_pregnancy hiv_status_mother
tab hiv_exposure

*Amoeba
tab katokatz_results parastie_identified_from_kato, m //repeats so keep only 1 and rename to coinfect_parasite for easy searching
rename katokatz_results coinfect_parasite
label variable coinfect_parasite "Collective parasite infection"
tab coinfect_parasite, m nol
// parasite_identified_from_kato gives exact parasites - most common is amoeba

*Schisto
tab serum_caa_result,m  // 3 missing, 33 positive
summ caa_concentration, d
bysort serum_caa_result: summ caa_concentration
list pid agemonth caa_concentration if serum_caa_result==.
//rename binary variable to coinfect_caa for easy seraching
rename serum_caa_result coinfect_caa
label variable coinfect_caa "Schisto infection (caa concentration > 30 in serum)"

gen any_infection = 0
replace any_infection = . if coinfect_caa == . | coinfect_parasite==.
replace any_infection = 1 if coinfect_caa == 1 | coinfect_parasite==1
tab any_infection, m
list coinfect* if any_infection ==.

tab coinfect_parasite coinfect_caa, m

*Age
gen ageweek = int(agedays/7)

*Age mother - generating age of mother at birth
gen age_mother_days = (age_mother*365.25) + 183
gen mother_age_at_birth = .
replace mother_age_at_birth = int((age_mother_days - agedays)/365.25)
list age_mother agemonth mother_age_at_birth in 1/10

*Cat birthweight
summarize birthweight, detail
gen bw_bucket = .
replace bw_bucket = 1 if birthweight < 1.5        
replace bw_bucket = 2 if birthweight >= 1.5 & birthweight < 2.5 
replace bw_bucket = 3 if birthweight >= 2.5
label define bw_cat 1 "Very low (<1.5)" 2 "<2.5" 3 "Normal (≥2.5)"
label values bw_bucket bw_cat
tab bw_bucket

*Breastfeeding at 4 months (because we have baseline data for longitudinal sample feeding status which is on average at 16 weeks so 4 months)
gen bf_4months = .
tab child_feeding_status, m
replace bf_4months = 1 if child_feeding_status==1
tab month_child_stopped_bf, m
replace bf_4months = 1 if month_child_stopped_bf>4 & month_child_stopped_bf<40
replace bf_4months=0 if month_child_stopped_bf<=4
replace bf_4months=0 if child_feeding_status>1 & ageweek<16
replace bf_4months=0 if child_feeding_status>1 & agemonth<4
tab child_feeding_status bf_4months, m //lots of missing data
tab agemonth if bf_4months==. 

tabulate vaccines_not_received 

*bOPV
rename dob_child infant_dob
gen age_bopv0_days = opv_date - infant_dob
gen delay_bopv0 = .
replace delay_bopv0 = 1 if age_bopv0_days > 7
replace delay_bopv0 = 0 if age_bopv0_days <= 7 & !missing(age_bopv0_days)
label define delay_lbl 0 "On time (≤7d)" 1 "Delayed (>7d)"
label values delay_bopv0 delay_lbl
tabulate delay_bopv0, m

**calculating age in days at the first bOPV vaccination 
gen age_opv0_days = opv_date - infant_dob
tabulate age_opv0_days
summarize age_opv0_days, detail

**weeks from birth till bOPV 4 at 14 weeks
**calculating age in weeks at the time of week 14 vaccination
gen age_opv14_days = bopv14w_date - infant_dob
gen age_opv14_weeks = round(age_opv14_days / 7, 0.1)
tabulate age_opv14_weeks
summarize age_opv14_weeks, detail

**calculating age in weeks at the time of IPV vaccine
gen ipv_days = ipv_date - infant_dob
gen ipv_weeks = round(ipv_days / 7, 0.1)
tabulate ipv_weeks
summarize ipv_weeks, detail 
list ipv_weeks if vaccines_not_received ==2
list ipv_weeks vaccines_not_received all_vaccines_received if pid=="POL 584"
gen ipv_bin = 1
replace ipv_bin = 0 if ipv_weeks==.
tab ipv_bin all_vaccines_received
tab ipv_bin vaccines_not_received
tab ipv_bin

drop opv_date bopv6w_date bopv10w_date ipv_date ipv_days age_opv0_days age_opv6_days age_opv10_days age_opv14_days

**ROTAVIRUS data
**Looking at the days/months since rotavirus vaccine. 
rename interview_date date
gen days_since_rota6 = date - rota_vaccine_w6
gen days_since_rota10 = date - rota_vaccine_w10
gen weeks_since_rota10 = int(days_since_rota10/7)
summarize days_since_rota6 days_since_rota10, detail

list pid rota_vaccine_w6 rota_vaccine_w10 if vaccines_not_received == 4
tab vaccines_not_received if pid=="POL 184"
replace vaccines_not_received = . if pid=="POL 184"
replace all_vaccines_received = 1 if pid=="POL 184"

gen rota_bin = 1
replace rota_bin = 0 if vaccines_not==4
tab rota_bin

*Make variable for follow-up length in days
gen follow_up_days = .
replace follow_up_days = date - bopv14w_date
list bopv14w_date date follow_up_days age_opv14_weeks in 1/10
histogram follow_up_days, bin(200)

*Flagging age outliers - will remove from future analysis
summ age_opv14_weeks, d
hist age_opv14_weeks, bin(100)
gen flag_ageoutlier = 0
replace flag_ageoutlier = 1 if age_opv14_weeks > 40 
tab flag_ageoutlier // 5 >40 weeks

*Data cleaning of mucosal immunity data:
**convert titers to log(2) scale. Higher titer indicates higher amount of antibodies. The oral poliovirus vaccine generates antibodies to PV1 and PV3 (types 1 and 3 poliovirus)
**Using the equation log2(x) = ln(x)/ln(2)

generate logPV1_stool = ln(pv1_stool)/ln(2)
generate logPV2_stool = ln(pv2_stool)/ln(2)
generate logPV3_stool = ln(pv3_stool)/ln(2)

generate logPV1_serum = ln(pv1_serum)/ln(2)
generate logPV2_serum = ln(pv2_serum)/ln(2)
generate logPV3_serum = ln(pv3_serum)/ln(2)

*Generating binary neuts data:
foreach x in 1 2 3 {
	gen pv`x'_stool_positive = 0
	replace pv`x'_stool_positive = 1 if pv`x'_stool>=8
	replace pv`x'_stool_positive = . if pv`x'_stool==.
	bysort pv`x'_stool_positive: summ pv`x'_stool 
	tab pv`x'_stool_positive if pv`x'_stool==., m
	gen pv`x'_stool_detectable = 0
	replace pv`x'_stool_detectable = 1 if pv`x'_stool>=4
	replace pv`x'_stool_detectable = . if pv`x'_stool==.
	bysort pv`x'_stool_detectable: summ pv`x'_stool
}
tab1 pv1_stool_positive pv2_stool_positive pv3_stool_positive pv1_stool_detectable pv2_stool_detectable pv3_stool_detectable
tab pv1_stool_positive pv3_stool_positive, col
//None missing stool data

//Tag if missing any serum neut data
 misstable summarize logPV1_serum logPV2_serum logPV3_serum
tab logPV1_serum if logPV3_serum == ., m
tab logPV3_serum if logPV1_serum == ., m
tab logPV2_serum if logPV3_serum == ., m
 //again same 5 particpants are missing on all
tab agemonth if logPV3_serum == ., m 
tab logPV3_stool if logPV3_serum==., m //All have low neuts - 1 detectable at 6.7 and rest are under 4.
 
bysort pid: egen flag_missing_serum_neuts = max(logPV1_serum), missing
recode flag_missing_serum_neuts 0/50=0 .=1
sort pid agemonth
list agemonth flag_missing_serum_neuts logPV1_serum logPV3_serum logPV3_serum if flag_missing_serum_neuts==1, sepby(pid)

//Binary variables to simplify regression analyses
replace minutes_towatersource = "15" if minutes_towatersource=="15 mintues"
destring minutes_towatersource, replace
summ minutes_, d

gen sampletype = 1
label define sampletype 1 "cross-sectional" 2 "longitudinal"
label values sampletype sampletype

gen follow_up_months = follow_up_days/30 

gen notbreastfeeding_bin = currently_breastfed
recode notbreastfeeding_bin 1=0 
label define notbreastfeeding_binl 0 "Exclusive breastfeeding" 1 "Not exclusive"
label values notbreastfeeding_bin notbreastfeeding_binl
tab notbreastfeeding_bin

gen watersource_bin = 0
replace watersource_bin = 1 if minutes_towatersource>=15 //upper quartile
label define watersource_binl 0 "<15 minutes" 1 ">=15 minutes, Upper quartile"
label values watersource_bin watersource_binl
tab watersource_bin

save "cross-sectional_clean", replace
use "cross-sectional_clean", clear

//////////////////////////////////////////////
*Merging in Ig data
use "stool_IgA_withlabids_studyvisitmonths.dta", clear
replace pid = "POL 482" if pid=="POL. 482"
replace pid = "POL 282" if pid=="POL 282 "
keep if study_visit_month=="aim2"
drop study_visit_month
duplicates tag pid, gen(duplicated_flag)
tab duplicated_flag

foreach x in 1 2 3{
	destring IPV`x', replace
	rename IPV`x' PV`x'_IgA_stool
	gen logPV`x'_IgA_stool = ln(PV`x'_IgA_stool)/ln(10)
}
list in 1/10, sepby(pid)
save "crosssectional_stool_IgA", replace

use "serum_IgG_1000f_with_labids_studyvisits.dta", clear
replace pid = "POL 482" if pid=="POL. 482"
keep if study_visit_month=="Aim2"
drop study_visit_month
duplicates tag pid, gen(duplicated_flag)
list if duplicated>0, sepby(pid)

foreach x in 1 2 3{
	destring IPV`x', replace
	bysort pid: egen PV`x'_IgG_serum  = mean(IPV`x')
	list in 1/10, sepby(pid)
	list if duplicated>0, sepby(pid)
	drop IPV`x'
	gen logPV`x'_IgG_serum = ln(PV`x'_IgG_serum )/ln(10)
}
list if duplicated>0, sepby(pid)
duplicates drop pid, force
list in 1/10, sepby(pid)
save "crosssectional_serum_IgG", replace

use "serum_IgA_with_labids_studyvisits.dta", clear
replace pid = "POL 482" if pid=="POL. 482"
keep if study_visit_month=="Aim2"
drop study_visit_month
duplicates tag pid, gen(duplicated_flag)
list if duplicated>0, sepby(pid)

foreach x in 1 2 3{
	replace IPV`x'="." if IPV`x'=="NaN"
	destring IPV`x', replace
	bysort pid: egen PV`x'_IgA_serum  = mean(IPV`x')
	list in 1/10, sepby(pid)
	list if duplicated>0, sepby(pid)
	drop IPV`x'
	gen logPV`x'_IgA_serum = ln(PV`x'_IgA_serum)/ln(10)
}
list if duplicated>0, sepby(pid)
duplicates drop pid, force
list in 1/10, sepby(pid)
save "crosssectional_serum_IgA", replace

use "crosssectional_serum_IgA", clear

merge 1:1 pid using "crosssectional_serum_IgG.dta"
drop _merge

merge 1:1 pid using "crosssectional_stool_IgA.dta"
list pid duplicate logPV1* if _merge!=3, sepby(pid)
drop _merge

save "crosssectional_Ig_data", replace

use "cross-sectional_clean", clear
merge 1:1 pid using "crosssectional_Ig_data"
list pid age_opv14_weeks duplicate logPV1* if _merge!=3, sepby(pid)
drop if _merge==2 //No demographic data available
drop _merge
list pid duplicate logPV1* in 1/10

save "cross-sectional_clean", replace

