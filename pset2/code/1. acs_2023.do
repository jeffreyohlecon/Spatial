	clear all
	set more off


		global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS2"
		global data "$path\Data"
		global results "$path\Results"

		/*
		import delimited "$data\csv_pus\psam_pusa.csv", clear
		save "$data\acs_personal.dta", replace

		import delimited "$data\csv_pus\psam_pusb.csv", clear
		append using "$data\acs_personal.dta"
		save "$data\acs_personal.dta", replace
		*/
		
		// Clean data: restrict the sample to people between 25 and 65yearsof agewhoareeither non-employed or employed in one of the sectors included in the analysis.
		use "$data\acs_personal.dta", clear
		
		keep if agep>=25 & agep<=65
		
		
		*Create a clean numeric version of NAICSP
		gen naicsp_clean = regexs(0) if regexm(naicsp, "^[0-9]+")  // Extract numeric prefix

		
		destring naicsp_clean, generate(naicsp_num) force  // Convert to numeric
		
		* Extract the first two/three digits of NAICS for broad sector classification
		gen naicsp_2digit = real(substr(naicsp_clean, 1, 2)) 
		gen naicsp_3digit = real(substr(naicsp_clean, 1, 3)) 

* Create sector classification variable
		gen sector_classification = ""

		// **Manufacturing (12 Sectors)**
		replace sector_classification = "01 - Food, Beverage, and Tobacco" if inrange(naicsp_3digit, 311, 312) // ok
		replace sector_classification = "02 - Textile, Apparel, Leather" if inrange(naicsp_3digit, 313, 316) // ok
		replace sector_classification = "02 - Textile, Apparel, Leather" if naicsp == "31M"
		
		
		replace sector_classification = "03 - Wood, Paper, Printing" if inrange(naicsp_3digit, 321, 323) // ok 
		replace sector_classification = "04 - Petroleum and Coal" if naicsp_3digit == 324 // ok
		replace sector_classification = "05 - Chemical" if naicsp_3digit == 325 // ok
		replace sector_classification = "06 - Plastics and Rubber" if naicsp_3digit == 326 // ok
		replace sector_classification = "07 - Nonmetallic Mineral" if naicsp_3digit == 327 // ok
		replace sector_classification = "08 - Primary and Fabricated Metal" if inrange(naicsp_3digit, 331, 332) // ok 
		replace sector_classification = "08 - Primary and Fabricated Metal" if naicsp == "33MS"
		
		replace sector_classification = "09 -  Machinery" if naicsp_3digit == 333 // ok
		replace sector_classification = "10 - Computer, Electronic, Electrical Equipment" if inrange(naicsp_3digit, 334, 335) // ok 
		replace sector_classification = "11 - Transportation Equipment" if naicsp_3digit == 336 //ok 
		replace sector_classification = "12 - Furniture, Miscellaneous Manufacturing" if inrange(naicsp_3digit, 337, 339) // ok 
		
		// **Wholesale & Retail Trade (Fixing the Issue)**
		replace sector_classification = "13 - Wholesale and Retail Trade" if inrange(naicsp_2digit, 42, 45) // ok 43 group doesn't exist anymore - it used to be a special kind of wholesale
		replace sector_classification = "13 - Wholesale and Retail Trade" if naicsp == "4MS"
		
		// **Construction**
		replace sector_classification = "14 - Construction" if naicsp_2digit == 23  // ok 

		// **Services (8 Sectors)**
		replace sector_classification = "15 - Transport Services" if inrange(naicsp_3digit, 481, 488) // ok 
		
		* replace sector_classification = "16 - Information Services" if inrange(naicsp_3digit, 511, 518) // ok Do we want to exclude 519? I don't think so! 
		
		replace sector_classification = "16 - Information Services" if naicsp_2digit == 51 // ok
		
		replace sector_classification = "17 - Finance and Insurance" if naicsp_2digit == 52 // ok 
		
		replace sector_classification = "18 - Real Estate" if naicsp_2digit == 53 // ok 
		
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
		
		tabout sector_classification using "$data/observations_sector.csv", style(csv) replace
		
		gen sector_cdp = real(substr(sector_classification, 1, 2))
				
		

		*keep if sector_classification!=""

		


				**********************************
		
		* Employment status (esr = employment status for 16+)
		*gen unemployed=.
		*replace unemployed=1 	if esr==3  
		*replace unemployed=0 	if esr!=3 & esr!=6 & esr!=.  // 6 = not in the labor force

		
		* Current state (state)
		
		* State of residence one year ago (mig = lived here one year ago, migsp = state where they lived 1 year ago)
		
		gen mover_house=.
		replace mover_house=1 	if mig==3
		replace mover_house=0 	if mig!=3 & mig!=2   // the model abstracts from international migration
		
		gen state_0=migsp
		replace state_0=. 		if state_0>=72 & state_0!=. // moved from other countries
		
		* I verified that this matches - mover_house & state_1 excluding 0
		
		* CDP:  "We find that around 2% of the U.S. population moves across states in a year in this time period."
		
		gen state_mover=.
		replace state_mover=0 	if state_0==state
		replace state_mover=1 	if state_0!=state & state_0!=0 & state_0!=.
		

		* Tabulate % of the U.S. population that moved across states
		* Set survey weights correctly
		svyset [pweight = pwgtp]

		* Generate weighted table
		recode state_mover(.=0) // only for tab purposes

		svy: tab state_mover, percent
		 
		// NOTE: ACS does not have information on workers’ past employment status or the industries in which people worked during the previous period => CPS data
		
		
	preserve
	
		* Step 1: Compute Industry Totals for Non-Movers (Diagonal y's)
		gen non_mover = 1 if state_0 == state  

		* Generate an integer weight variable from pwgtp (since collapse does not support pweight)
		gen weight = round(pwgtp)

		collapse (sum) non_mover [fw=weight], by(state sector_classification)

		gen sector_cdp = real(substr(sector_classification, 1, 2))


		* Create the column variable based on sector_cdp mapping
		gen sector_col = ""

		* Map sector_cdp values to "c1", "c2", ..., "c23"
		forvalues i = 0/22 {   // Assuming sector_cdp values are from 0 to 22
			replace sector_col = "c" + string(`i' + 1) if sector_cdp == `i'
		}

		
		save "$results/acs_nonmover_totals.dta", replace

		
	restore
	
	
		* Keep only those who moved across states
		keep if state_mover == 1 & state != state_0

		* Convert pweights to integer fweights
		gen weight_int = round(pwgtp)

		* Collapse the data by origin state_0, destination state, and current sector_classification
		collapse (sum) mover_count=state_mover [fw=weight_int], by(state_0 state sector_classification)
		gen sector_cdp = real(substr(sector_classification, 1, 2))


		* Create the column variable based on sector_cdp mapping
		gen sector_col = ""

		* Map sector_cdp values to "c1", "c2", ..., "c23"
		forvalues i = 0/22 {   // Assuming sector_cdp values are from 0 to 22
			replace sector_col = "c" + string(`i' + 1) if sector_cdp == `i'
		}

		* Save the collapsed dataset
		save "$results/acs_off_diagonal_totals.dta", replace
			
			
	