clear all
set more off

* Define list of available states
local state_list "1 2 4 5 6 8 9 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56"

* Get the number of states
local N_states : word count `state_list'
local N_sectors = 23  // Number of sectors

* Initialize an empty Mata matrix of size (N_states * 23) × (N_states * 23)
mata: M = J(`N_states' * `N_sectors', `N_states' * `N_sectors', 0)

* Convert state_list into a Stata-local list for indexing
local i = 1
foreach origin of local state_list {
    local j = 1
    foreach dest of local state_list {
        
        * Load the corresponding transition matrix
        capture confirm file "$results/interstate_matrix_`origin'_to_`dest'.dta"
        if _rc == 0 {
            use "$results/interstate_matrix_`origin'_to_`dest'.dta", clear
        }
        else {
            * If no matrix exists, create a 23x23 zero matrix
            clear
            set obs `N_sectors'
            forvalues k = 1/`N_sectors' {
                gen c`k' = 0
            }
        }

        * Convert the loaded transition matrix into Mata format
        mata: X = st_data(., "c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15 c16 c17 c18 c19 c20 c21 c22 c23")

        * Assign it to the corresponding block in the large matrix
        mata: M[(`i'-1) * `N_sectors' + 1 .. `i' * `N_sectors', (`j'-1) * `N_sectors' + 1 .. `j' * `N_sectors'] = X

        local j = `j' + 1
    }
    local i = `i' + 1
}

* Store the final stacked matrix into Stata
mata: st_matrix("final_matrix", M)

* Convert the final matrix into a dataset
clear
svmat final_matrix, names(col)

* Save the full block matrix
save "$results/full_transition_matrix.dta", replace








clear all
set more off

* Load the full transition matrix
use "$results/full_transition_matrix.dta", clear

* Compute the row totals
egen row_sum = rowtotal(c1-c1150)

* Normalize each column by dividing by the row sum
foreach var of varlist c1-c1150 {
    replace `var' = `var' / row_sum if row_sum > 0
}

* Drop the row_sum variable (optional)
drop row_sum

* Save the normalized transition probability matrix
save "$results/mu.dta", replace


