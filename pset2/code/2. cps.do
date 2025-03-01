	clear all
	set more off


		global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS2"
		global data "$path\Data\IPUMS"
		global results "$path\Results"

		
		* Execute the IPUMS codebook
		/*
		do "$data\cps_00004.do"		
		drop if year==2020
		save "$data\cps_2021_2023.dta", replace
		*/
		use "$data\cps_2021_2023.dta", clear


		* Generate three-month forward matching variables
		gen month_3 = month + 3

		* Adjust for months that exceed 12
		replace month_3 = month_3 - 12 if month_3 > 12

		* Generate year_3 variable (shift forward only when crossing into a new year)
		gen year_3 = year
		replace year_3 = year_3 + 1 if month_3 < month
		
		sort cpsid cpsidp year_3 month_3

		save "$data\cps_data_temp.dta", replace
	

		use "$data\cps_data_temp.dta", clear
		drop month_3 year_3 
		rename (empstat ind) (empstat_3 ind_3)
		rename (year month) (year_3 month_3)
		
		sort cpsid cpsidp year_3 month_3
		
		save "$data\cps_data_temp_clean.dta", replace				
		
		
		use "$data\cps_data_temp.dta", clear
		merge cpsid cpsidp year_3 month_3 using "$data\cps_data_temp_clean.dta"
		keep if _merge==3
		

				keep if age>=25 & age<=65

		
		* Employment transitions 
		gen emp_transition = .

		* Stayed employed
		replace emp_transition = 1 if inlist(empstat, 10, 12) & inlist(empstat_3, 10, 12) 

		* Employed → Unemployed
		replace emp_transition = 2 if inlist(empstat, 10, 12) & inlist(empstat_3, 21, 22)

		* Employed → NILF
		replace emp_transition = 3 if inlist(empstat, 10, 12) & inlist(empstat_3, 32, 34, 36)

		* Unemployed → Employed
		replace emp_transition = 4 if inlist(empstat, 21, 22) & inlist(empstat_3, 10, 12)

		* NILF → Employed
		replace emp_transition = 5 if inlist(empstat, 32, 34, 36) & inlist(empstat_3, 10, 12)

		* Stayed Unemployed
		replace emp_transition = 6 if inlist(empstat, 21, 22) & inlist(empstat_3, 21, 22)

		* Unemployed → NILF
		replace emp_transition = 7 if inlist(empstat, 21, 22) & inlist(empstat_3, 32, 34, 36)

		* NILF → Unemployed
		replace emp_transition = 8 if inlist(empstat, 32, 34, 36) & inlist(empstat_3, 21, 22)

		* Stayed NILF
		replace emp_transition = 9 if inlist(empstat, 32, 34, 36) & inlist(empstat_3, 32, 34, 36)

		* Exclude individuals who are in the Armed Forces or NIU (Not in Universe)
		drop if empstat == 0 | empstat == 1 | empstat_3 == 0 | empstat_3 == 1
		
		
		
		* Check https://cps.ipums.org/cps/codes/ind_2020_codes.shtml 


		
		gen naics_sector = ""

				
		* 1. Food, Beverage, and Tobacco Products
		replace naics_sector = "01 - Food, Beverage, and Tobacco Products" if inlist(ind, 1070, 1080, 1090, 1170, 1180, 1190, 1270, 1280, 1290, 1370, 1390)

		* 2. Textiles and Apparel Products
		replace naics_sector = "02 - Textiles and Apparel Products" if inlist(ind, 1470, 1480, 1490, 1570, 1590, 1670, 1691, 1770, 1790)

		* 3. Wood, Paper, Printing, and Related Products
		replace naics_sector = "03 - Wood, Paper, Printing and Related Products" if inlist(ind, 1870, 1880, 1890, 1990, 3770, 3780, 3790, 3875)

		* 4. Petroleum and Coal Products
		replace naics_sector = "04 - Petroleum and Coal Products" if inlist(ind, 2070, 2090)

		* 5. Chemicals
		replace naics_sector = "05 - Chemicals" if inlist(ind, 2170, 2180, 2190, 2270, 2280, 2290)

		* 6. Plastics and Rubber Products
		replace naics_sector = "06 - Plastics and Rubber Products" if inlist(ind, 2370, 2380, 2390)

		* 7. Nonmetallic Mineral Products
		replace naics_sector = "07 - Nonmetallic Mineral Products" if inlist(ind, 2470, 2480, 2490, 2570, 2590)

		* 8. Primary and Fabricated Metal Products
		replace naics_sector = "08 - Primary and Fabricated Metal Products" if inlist(ind, 2670, 2680, 2690, 2770, 2780, 2790, 2870, 2880, 2890, 2970, 2980, 2990)

		* 9. Machinery
		replace naics_sector = "09 - Machinery" if inlist(ind, 3070, 3080, 3095, 3170, 3180, 3291)

		* 10. Computer, Electrical, and Appliance
		replace naics_sector = "10 - Computer, Electrical, and Appliance" if inlist(ind, 3365, 3370, 3380, 3390, 3470, 3490)

		* 11. Transportation Equipment
		replace naics_sector = "11 - Transportation Equipment" if inlist(ind, 3570, 3580, 3590, 3670, 3680, 3690)

		* 12. Furniture and Miscellaneous Products
		replace naics_sector = "12 - Furniture and Miscellaneous Products" if inlist(ind, 3895, 3960, 3970, 3980, 3990)

		* 13. Transport Services
		replace naics_sector = "13 - Transport Services" if inlist(ind, 6070, 6080, 6090, 6170, 6180, 6190, 6270, 6280, 6290, 6370, 6380, 6390)

		* 14. Information and Communication
		replace naics_sector = "14 - Information and Communication" if inlist(ind, 6470, 6480, 6490, 6570, 6590, 6670, 6672, 6680, 6690, 6695, 6770, 6780)

		* 15. Finance and Insurance
		replace naics_sector = "15 - Finance and Insurance" if inrange(ind, 6870, 6992)

		* 16. Real Estate
		replace naics_sector = "16 - Real Estate" if inrange(ind, 7071, 7080)

		* 17. Education
		replace naics_sector = "17 - Education" if inrange(ind, 7860, 7890)

		* 18. Health Care
		replace naics_sector = "18 - Health Care" if inrange(ind, 7970, 8290)

		* 19. Accommodation and Food Services
		replace naics_sector = "19 - Accommodation and Food Services" if inrange(ind, 8660, 8690)

		* 20. **Other Services (Updated Classification)**
		replace naics_sector = "20 - Other Services" if inlist(ind, 0170, 0180, 0190, 0270, 0280, 0290, 0370, 0380, 0390, 0470, 0490, 7790, 8070, 8080, 8090, 8170, 8180, 8380, 8390, 8470, 8570, 8580, 8590, 8970, 8980, 8990, 9070, 9080, 9090, 9160, 9170, 9180, 9190)

		* 21. Wholesale and Retail Trade
		replace naics_sector = "21 - Wholesale and Retail Trade" if inrange(ind, 4070, 5790)

		* 22. Construction
		replace naics_sector = "22 - Construction" if ind == 0770

		* First, define "Unemployed" as a separate category
		replace naics_sector = "00 - Unemployed" if inlist(empstat, 21, 22)

				
		gen naics_sector_3 = ""
	


		replace naics_sector_3 = "01 - Food, Beverage, and Tobacco Products" if inlist(ind_3, 1070, 1080, 1090, 1170, 1180, 1190, 1270, 1280, 1290, 1370, 1390)

		replace naics_sector_3 = "02 - Textiles and Apparel Products" if inlist(ind_3, 1470, 1480, 1490, 1570, 1590, 1670, 1691, 1770, 1790)

		replace naics_sector_3 = "03 - Wood, Paper, Printing and Related Products" if inlist(ind_3, 1870, 1880, 1890, 1990, 3770, 3780, 3790, 3875)

		replace naics_sector_3 = "04 - Petroleum and Coal Products" if inlist(ind_3, 2070, 2090)

		replace naics_sector_3 = "05 - Chemicals" if inlist(ind_3, 2170, 2180, 2190, 2270, 2280, 2290)

		replace naics_sector_3 = "06 - Plastics and Rubber Products" if inlist(ind_3, 2370, 2380, 2390)

		replace naics_sector_3 = "07 - Nonmetallic Mineral Products" if inlist(ind_3, 2470, 2480, 2490, 2570, 2590)

		replace naics_sector_3 = "08 - Primary and Fabricated Metal Products" if inlist(ind_3, 2670, 2680, 2690, 2770, 2780, 2790, 2870, 2880, 2890, 2970, 2980, 2990)

		replace naics_sector_3 = "09 - Machinery" if inlist(ind_3, 3070, 3080, 3095, 3170, 3180, 3291)

		replace naics_sector_3 = "10 - Computer, Electrical, and Appliance" if inlist(ind_3, 3365, 3370, 3380, 3390, 3470, 3490)

		replace naics_sector_3 = "11 - Transportation Equipment" if inlist(ind_3, 3570, 3580, 3590, 3670, 3680, 3690)

		replace naics_sector_3 = "12 - Furniture and Miscellaneous Products" if inlist(ind_3, 3895, 3960, 3970, 3980, 3990)
		
		replace naics_sector_3 = "13 - Transport Services" if inlist(ind_3, 6070, 6080, 6090, 6170, 6180, 6190, 6270, 6280, 6290, 6370, 6380, 6390)
		
		replace naics_sector_3 = "14 - Information and Communication" if inlist(ind_3, 6470, 6480, 6490, 6570, 6590, 6670, 6672, 6680, 6690, 6695, 6770, 6780)

		replace naics_sector_3 = "15 - Finance and Insurance" if inrange(ind_3, 6870, 6992)
		
		replace naics_sector_3 = "16 - Real Estate" if inrange(ind_3, 7071, 7080)

		replace naics_sector_3 = "17 - Education" if inrange(ind_3, 7860, 7890)

		replace naics_sector_3 = "18 - Health Care" if inrange(ind_3, 7970, 8290)

		replace naics_sector_3 = "19 - Accommodation and Food Services" if inrange(ind_3, 8660, 8690)

		replace naics_sector_3 = "20 - Other Services" if inlist(ind_3, 0170, 0180, 0190, 0270, 0280, 0290, 0370, 0380, 0390, 0470, 0490, 7790, 8070, 8080, 8090, 8170, 8180, 8380, 8390, 8470, 8570, 8580, 8590, 8970, 8980, 8990, 9070, 9080, 9090, 9160, 9170, 9180, 9190)

		replace naics_sector_3 = "21 - Wholesale and Retail Trade" if inrange(ind_3, 4070, 5790)

		replace naics_sector_3 = "22 - Construction" if ind_3 == 0770
	
		replace naics_sector_3 = "00 - Unemployed" if inlist(empstat_3, 21, 22)
		
		
				gen sector_cdp_1 = real(substr(naics_sector, 1, 2))

				gen sector_cdp_3 = real(substr(naics_sector_3, 1, 2))

		
		gen industry_transition = "Same Industry"
		replace industry_transition = "Industry Change" if naics_sector != naics_sector_3
		
		tab naics_sector
		tab naics_sector_3
		tab industry_transition  // excludes unemployment
		
		gen sector_transition = naics_sector + " → " + naics_sector_3
		
		gen sector = real(substr(naics_sector, 1, 2))
		gen sector_3 = real(substr(naics_sector_3, 1, 2))

	* Export transition matrices
	
	preserve

local categories "Unemployed Food_Beverage_Tobacco Textiles_Apparel Wood_Paper_Printing Petroleum_Coal Chemicals Plastics_Rubber Nonmetallic_Mineral Primary_Fabricated_Metal Machinery Computer_Electrical_Appliance Transportation_Equipment Furniture_Miscellaneous Transport_Services Information_Communication Finance_Insurance Real_Estate Education Health_Care Accommodation_Food_Services Other_Services Wholesale_Retail_Trade Construction"

* Get the unique state codes that actually exist in the dataset
levelsof statefip, local(state_list)

foreach s in `state_list' {  // Loop over existing state codes only
    di "Processing State: `s'..."
    
    keep if statefip == `s'

    * Ensure all 23 categories appear even if zero transitions exist
    tempname transition_matrix
    matrix `transition_matrix' = J(23, 23, 0)  // 23×23 matrix filled with zeros
    
    * Populate the transition matrix	
    tabulate naics_sector naics_sector_3, matcell(temp_matrix)

    * Fill the pre-created 23×23 matrix with real values
    forvalues i = 1/23 {
        forvalues j = 1/23 {
            matrix `transition_matrix'[`i', `j'] = temp_matrix[`i', `j']
        }
    }

    * Convert the matrix to a dataset
    clear
    svmat `transition_matrix', names(col)

    * Add sector names as row labels
    gen row_sector = ""
    forvalues i = 1/23 {
        replace row_sector = word("`categories'", `i') in `i'
    }

    * Save the state-specific transition matrix
    gen statefip = `s'
    save "$results/transition_matrix_`s'.dta", replace

    restore, preserve
}

	



foreach s in `state_list' {
    di "Processing State Probabilities: `s'..."
    
    * Load the transition matrix
    use "$results/transition_matrix_`s'.dta", clear

    drop row_sector statefip  // Drop unnecessary variables

    * Compute column sums for normalization
    foreach var of varlist c1-c23 {
        qui sum `var'  // Get column total
        local col_total = r(sum)  // Store the total
        
        * Normalize each cell by dividing by the column total
        if `col_total' > 0 {
            replace `var' = `var' / `col_total'
        }
    }

    * Save the probability transition matrix
    save "$results/transition_matrix_probs_`s'.dta", replace
}



		// Those above are the within state matrix transitions but we are missing past industry of interstate movers.
		
		// Full transition matrix: 
		
		// 1. Assume workers who move across states and in second period are in state i and sector j have past industry = workers who did not switch states and are in the second period in state i sector j
		
		