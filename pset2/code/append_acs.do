	clear all
	set more off


		/// global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS2"
		/// global data "$path\Data"
		/// global results "$path\Results"

	global path "/Users/henriquemota/Dropbox/BigDataFiles/pset2"
	global data "$path/ACS"
	
	
    import delimited "/Users/henriquemota/Dropbox/BigDataFiles/pset2/ACS/psam_pusa.csv", clear
    save "$data/acs_personal_raw.dta", replace

    import delimited "$data/psam_pusb.csv", clear
    append using "$data/acs_personal_raw.dta"
    save "$data/acs_personal_raw.dta", replace

	
	
