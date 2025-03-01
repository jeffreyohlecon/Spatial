clear all
	set more off


		global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS2"
		global data "$path\Data"
		global results "$path\Results"
		
	

	use "$results/acs_nonmover_totals.dta", clear
	
	levelsof state, local(state_list)


	

foreach s in `state_list' {  
   
            * Load ACS nonmovers totals for this specific transition
            use "$results/acs_nonmover_totals.dta", clear
			keep if state == `s'
            keep sector_col non_mover  

            * Convert to row vector
            capture drop temp
            gen temp = 1  
            reshape wide non_mover, i(temp) j(sector_col) string	

            * Rename columns
            foreach var of varlist non_moverc* {
                local newvar = substr("`var'", 10, .)  
                rename `var' `newvar'
            }
            
            drop temp
            tempfile acs_data
            save `acs_data', replace

            * Load transition probabilities for state `s`
            use "$results/transition_matrix_probs_`s'.dta", clear
            tempfile transition_matrix
            save `transition_matrix', replace

            * Append ACS totals as last row
            use `transition_matrix', clear
            append using `acs_data', force  

            local row_count = _N  // Get row number of the last row (mover totals)

            * Multiply transition probabilities by the corresponding ACS mover count
            foreach var of varlist c1-c23 {
                replace `var' = `var' * `var'[`row_count'] if _n != `row_count'
            }

            * Replace all missing values with 0
            foreach var of varlist c1-c23 {
                replace `var' = 0 if missing(`var')
            }

            * Drop the last row (mover_count row is not needed anymore)
            drop if _n == `row_count'

            * Save the final adjusted transition matrix
            save "$results/interstate_matrix_`s'_to_`s'.dta", replace
        }

    
	
	