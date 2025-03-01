	clear all
	set more off


		global path "C:\Users\yulia\OneDrive - The University of Chicago\Winter 2025\Spatial Economics\Problem Sets\PS2"
		global data "$path\Data"
		global results "$path\Results"
		
		

	
	clear all

* Step 1: Identify all valid (state_0, state) pairs with transitions
use "$results/acs_off_diagonal_totals.dta", clear
gen pair = string(state_0) + "_" + string(state)
levelsof pair, local(valid_pairs)  // Store all valid pairs

levelsof state, local(state_list)


* Step 2: Process each state pair
foreach s in `state_list' {  
    foreach origin of numlist 1/56 {
        if `origin' == `s' continue  // Skip intra-state transitions (diagonal)

        local current_pair "`origin'_`s'"  // Create state pair identifier

        * Check if this pair is in the list of valid transitions
        local is_valid = 0
        foreach p in `valid_pairs' {
            if "`p'" == "`current_pair'" local is_valid = 1
        }

        di "Processing Transition Matrix: `origin' → `s'..."

        if `is_valid' {  
            * Step 3: Load ACS movers totals for this specific transition
            use "$results/acs_off_diagonal_totals.dta", clear
            keep if state == `s' & state_0 == `origin'
            keep sector_col mover_count  

            * Convert to row vector
            capture drop temp
            gen temp = 1  
            reshape wide mover_count, i(temp) j(sector_col) string	

            * Rename columns
            foreach var of varlist mover_countc* {
                local newvar = substr("`var'", 12, .)  
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
            save "$results/interstate_matrix_`origin'_to_`s'.dta", replace
        }
        else {
            * Step 4: Create a 23×23 zero matrix if no transitions exist
            di "No recorded movers from `origin' to `s'. Creating zero matrix..."
            clear
            set obs 23  // 23 rows (each industry)
            foreach var of numlist 1/23 {
                gen c`var' = 0  // Generate all-zero columns
            }

            * Save zero transition matrix
            save "$results/interstate_matrix_`origin'_to_`s'.dta", replace
        }
    }
}

	