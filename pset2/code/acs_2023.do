clear all
set more off


		/// global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS2"
		/// global data "$path\Data"
		/// global results "$path\Results"

	global path "/Users/henriquemota/Dropbox/BigDataFiles/pset2"
	global data "$path/ACS"
	

		use "$data/acs_personal_raw.dta"
		
		// Clean data: restrict the sample to people between 25 and 65yearsof agewhoareeither non-employed or employed in one of the sectors included in the analysis.
		
		keep if agep>=25 & agep<=65
		
		
		*Create a clean numeric version of NAICSP
		gen naicsp_clean = regexs(0) if regexm(naicsp, "^[0-9]+")  // Extract numeric prefix - Generates the exact same numbers of missing naics


		replace naicsp_clean = substr(naicsp, 1, 3) if regexm(naicsp, "^[0-9]{3}[A-Z]")  // Keep 3-digit codes
		replace naicsp_clean = substr(naicsp, 1, 4) if regexm(naicsp, "^[0-9]{4}[A-Z]")  // Keep 4-digit codes
		replace naicsp_clean = substr(naicsp, 1, 5) if regexm(naicsp, "^[0-9]{5}[A-Z]")  // Keep 5-digit codes

		* destring naicsp_clean, generate(naicsp_num) force  // Convert to numeric

		* Extract the first two/three digits of NAICS for broad sector classification
		gen naicsp_2digit = real(substr(naicsp_clean, 1, 2)) 
		gen naicsp_3digit = real(substr(naicsp_clean, 1, 3)) 

		* Create sector classification variable
		gen sector_classification = ""

		// **Manufacturing (12 Sectors)**
		replace sector_classification = "01 - Food, Beverage, and Tobacco" if inrange(naicsp_3digit, 311, 312) // ok
		replace sector_classification = "02 - Textile, Apparel, Leather" if inrange(naicsp_3digit, 313, 316) // ok
		replace sector_classification = "03 - Wood, Paper, Printing" if inrange(naicsp_3digit, 321, 323) // ok 
		replace sector_classification = "04 - Petroleum and Coal" if naicsp_3digit == 324 // ok
		replace sector_classification = "05 - Chemical" if naicsp_3digit == 325 // ok
		replace sector_classification = "06 - Plastics and Rubber" if naicsp_3digit == 326 // ok
		replace sector_classification = "07 - Nonmetallic Mineral" if naicsp_3digit == 327 // ok
		replace sector_classification = "08 - Primary and Fabricated Metal" if inrange(naicsp_3digit, 331, 332) // ok 
		replace sector_classification = "09 -  Machinery" if naicsp_3digit == 333 // ok
		replace sector_classification = "10 - Computer, Electronic, Electrical Equipment" if inrange(naicsp_3digit, 334, 335) // ok 
		replace sector_classification = "11 - Transportation Equipment" if naicsp_3digit == 336 //ok 
		replace sector_classification = "12 - Furniture, Miscellaneous Manufacturing" if inrange(naicsp_3digit, 337, 339) // ok 
		
		// **Wholesale & Retail Trade (Fixing the Issue)**
		replace sector_classification = "13 - Wholesale and Retail Trade" if inrange(naicsp_2digit, 42, 45) // ok 
		
		// **Construction**
		replace sector_classification = "14 - Construction" if naicsp_2digit == 23  // ok 

		// **Services (8 Sectors)**
		replace sector_classification = "15 - Transport Services" if inrange(naicsp_3digit, 481, 488) // ok 
		replace sector_classification = "16 - Information Services" if inrange(naicsp_3digit, 511, 518) // ok
		replace sector_classification = "17 - Finance and Insurance" if inrange(naicsp_3digit, 521, 525) // ok 
		replace sector_classification = "18 - Real Estate" if inrange(naicsp_3digit, 531, 533) // ok 
		replace sector_classification = "19 - Education" if naicsp_2digit == 61 // ok 
		replace sector_classification = "20 - Health Care" if inrange(naicsp_3digit, 621, 624) //ok 
		replace sector_classification = "21 - Accommodation and Food Services" if inrange(naicsp_3digit, 721, 722) //ok 
		replace sector_classification = "22 - Other Services" if inlist(naicsp_3digit, 493, 541, 561, 562, 711, 712, 713, 811, 812, 813, 814) // ok, actually mistake at 55...
		replace sector_classification = "22 - Other Services" if naicsp_2digit == 55 // ok
		

//		keep if sector_classification!="" // This is probably a mistake 

		replace sector_classification = "00 - Unemployed" if naicsp == "999920" 

		
		* Dropping the employed in sectors that are not classified according to CDP 
		
		gen not_CDP = 0
		replace not_CDP = 1 if sector_classification == "" && naicsp != ""
		
		tabout naicsp if not_CDP == 1 using "$data/not_CDP.csv", style(csv) replace
		
		drop if not_CDP == 1
		drop not_CDP 
		
		tab esr // All the Armed Forces members
		
		* Generate unemployed (Oddly enought the original data does not have it)
		replace sector_classification = "00 - Unemployed" if sector_classification == "" 
		replace sector_classification = "00 - Unemployed" if esr == 6 | esr == 3
		
		
		tabstat pwgtp, stat(sum) by(sector_classification)
		
		
		gen sector_cdp = real(substr(sector_classification, 1, 2))
		
		
		*** Clean State 
		
		merge m:1 state using "$data/state_code.dta"
		drop _merge
		
		
		save "$data/acs_personal_clean.dta", replace
		
		clear all 
		
		* Load your original dataset
		use "$data/acs_personal_clean.dta", clear

		* Step 1: Extract unique states and save as a temporary file
		preserve
			keep state_new 
			duplicates drop
			tempname temp_states
			tempfile temp_states_file
			save `temp_states_file'
		restore

		* Step 2: Extract unique sectors and save as a temporary file
		preserve
			keep sector_cdp
			duplicates drop
			tempname temp_sectors
			tempfile temp_sectors_file
			save `temp_sectors_file'
		restore

		* Step 3: Generate partially right collapse 

		collapse  (first) state_name (first) sector_classification (sum) employment = pwgtp ,by(state_new sector_cdp) 
				
				egen temp = sum(employment)
				gen L0 = employment/temp
				drop temp
				
		tempname temp_collapsed
		tempfile temp_collapsed_file		
		save `temp_collapsed_file'

		use `temp_states_file', clear
		cross using `temp_sectors_file'   // Cross-join the temporary sector file

		* Step 4: Merge with the original dataset
		merge m:1 state_new sector_cdp using `temp_collapsed_file'

		replace employment = 0 if employment == .
		replace L0 = 0 if L0 == . 
		
		bysort state_new : gen NAME = state_name[_N]
		replace state_name = NAME if state_name =="" 
		
		bysort sector_cdp : gen SECTOR = sector_classification[_N]
		replace sector_classification = SECTOR if sector_classification =="" 
		
		sort state_new sector_cdp 
		
		drop NAME SECTOR 
		
		
		gen state_sector = state_name + "-" + sector_classification
		
		* Diagnosis, zero sector.
		
		tabout state_sector if _merge == 1 using "$data/zero_sectors.csv", style(csv) replace
		
		drop _merge state_sector 
		
		save "$path/L_0/L0_complete.dta", replace 
		
		preserve
		
		keep L0
		
		export delimited using "$path/L_0/L0.csv", replace
		
		restore 
		
		
		