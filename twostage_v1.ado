/*
twostage.ado file

Estimation of two-stage models for individual participant data.

version 1.1 

with prefix command

Author: R Thiesmeier
*/ 

cap prog drop twostage
program define twostage, properties(prefix)
	version 17
	local version : di "version " string(_caller()) ":"

	set prefix twostage
	
	capture _on_colon_parse `0'
	if c(rc) | `"`s(after)'"' == "" {
		if !c(rc) {
			local 0 `"`s(before)'"'
		}
		if replay() {
			if "`e(prefix)'" != "twostage" {
				error 301
			}
			`version' Display `0'
			exit
		}
	}

	`version' TwoStage `0'
end 

cap prog drop TwoStage
prog def TwoStage , properties(mi) eclass
	version 17
	local version : di "version " string(_caller()) ":"
	local cmdline : copy local 0

	_on_colon_parse `0'
	local command `"`s(after)'"'
	local 0 `"`s(before)'"'
	
	/*
	if replay() {
		if ("`e(cmd)'" != "twostage") {
			di `0'
		}
	}
	*/
	
syntax [if] [in]  ///
	[, ID(string) ///
	FRame(string) /// 
	fixed /// 
	random(string) *]
	
	************ checks *************
	
	// check observations
	
	marksample touse, strok novarlist
	
	qui count if `touse'
	local noobs = r(N)
	local N = _N
	
	if `noobs' == 0 {
		di as err "no observations"
		exit 2000
	}
	
	if `noobs' == 1 {
		di as err "insufficient number of observations"
		exit 2001
	}
	
	// check ID and number of studies; ID cannot be optional
	
	if "`id'" == "" {
		di as err "{bf:id} is missing"
		exit 198
	}
	
	qui tab `id' if `touse'
		local studies = r(r)
			if (`studies' == 1) {
				di in red "More than one study is required"
				exit 198
			}
	
	// run the model after ":" once

	local model `command'
	
	qui `model'
		if _rc != 0 {
			di as err `"Model specification is wrong. Please check the model"'
			exit 198
		}
	
		local title = e(title)
		local cmdname = e(cmd)
		local cmdline = e(cmdline)
		local twostage_predict = e(predict)
		local twostage_estat = e(estat_cmd) 
		local df_m = e(df_m)		
		local cmd = e(cmd) 
		mac shift 
		local cmdbody = "`*'" 
		if regexm("`cmdline'", ",") != 1  local cmdline "`cmdline' , "
		
		local row_names : rownames e(b)
		local col_names : colnames e(b)
			
		local nr_p : colsof e(b) 
		mat est_beta = e(b)
		
		local list_name_b ""
		forv i = 1/`nr_p' {
						local list_name_b "`list_name_b' y`i'"
		}
		
		local list_name_V ""
		local tot_nr_V = `nr_p'*(`nr_p'+1)/2
		forv i = 1/`tot_nr_V' {
						local list_name_V "`list_name_V' V`i'"
		}
		
							// names for frame 
							local col_names_fr "`col_names'"
							local col_names_fr : subinstr local col_names_fr "_cons" "" , all
							local c_ : word count `col_names_fr'
							local _names_fr_se ""
							foreach w of local col_names_fr {
								local _names_fr_se "`_names_fr_se' `w'_se"
							}

	************ twostage *************	
	
	// Step 1 
	
	local c = 1
	quietly levelsof `id' if `touse', local(l_studies) // check
	foreach s of local l_studies {
		
		if "`if'" != "" {
						local ourif `"`if' & `id'==`s'"'
						local modelc : subinstr local model "`if'" "`ourif'"

		}
		else {
				local ourif `" if `id'==`s'"'
				local modelc : subinstr local model "," " `ourif' , "
		}
		
		quietly `modelc'
		estimates store _est_study_`s'
		mat e_b = e(b)	
		mat e_V = e(V)
		local df_m = e(df_m)
		
				// post study specifc estimates into frame			
				local post_beta ""
				local post_beta_se ""
				forv k = 1 / `c_'{
					local post_beta = "`post_beta' (`=r(table)[1,`k']')"
					local post_beta_se = "`post_beta_se' (`=r(table)[2,`k']')"
				}
				
				if `c' == 1 {
							if "`frm'" == ""{
								cap frame drop _est_studies_2stage
								frame create _est_studies_2stage  `id' `col_names_fr' `_names_fr_se'
								local frame_name "_est_studies_2stage"
							}
							else if "`frm'" != "" {
								cap frame drop `frm'
								frame create `frm'  `id' `col_names_fr' `_names_fr_se'
								local frame_name "`frm'"
							}
				}
				
				frame post _est_studies_2stage (`s') `post_beta' `post_beta_se'
			
		if `c'==1 {
					cap frame drop _est_fr_two_stage 
					frame create _est_fr_two_stage  study `list_name_b' `list_name_V'
		}

			local post_est_b ""
			forv i = 1/`nr_p' {
				local post_est_b "`post_est_b' (`=e_b[1,`i']')"
			}
			local post_est_V ""
			forv i = 1/`nr_p' {
					forv j = `i'/`nr_p' {
						local post_est_V "`post_est_V' (`=e_V[`j',`i']')"
					}
			}
		
		    frame post _est_fr_two_stage (`s') `post_est_b' `post_est_V'
			local c = `c' + 1
			
	}
	
	// Step 2 
	
    frame _est_fr_two_stage  {
		qui ds 
		local list_var "`r(varlist)'"
		
		foreach x of local list_var {
			if `x'[1] == 0 drop `x'
		}
		
		// meta mvregress with option for fixed or random effect
			
		if "`fixed'" == "fixed" {
		
				qui meta mvregress y* ,  wcovvariables(V*) fixed `options'
				
				mat est_b_mv = e(b)
				mat est_V_mv = e(V)
				
				local effect `fixed'
		}
		
		else {
		
				qui meta mvregress y* ,  wcovvariables(V*) random(`random') `options'
				
				mat est_b_mv = e(b)[1, 1..`nr_p']
				mat est_V_mv = e(V)[1..`nr_p', 1..`nr_p']
				
				local effect "`random'"
		}
		
			local studies = _N
	} 
	
	
	// Step 3
	
	mat est_b_c = e_b 
	mat est_V_c = e_V
	
	local nr_c_mv : rowsof est_V_mv
	
	local f = 1
	forv i = 1/ `nr_p' {
		if est_b_c[1,`i'] != 0 {
						mat est_b_c[1, `i'] = est_b_mv[1, `f']
						local f = `f' + 1
		}
	}
 
	local f = 1
	local p = 1
 	forv i = 1/`nr_p' {
					forv j = `i'/`nr_p' {
					
						if est_V_c[`j',`i'] != 0 {
												 
												mat est_V_c[`j',`i'] = est_V_mv[`f',`p']
												
												if `f' != `nr_c_mv' {
															local f = `f' + 1
												}
												else {
														local f = `p' + 1
														local p = `p' + 1
												}
						}					
					}
	}
	
 	mata: clean_V = makesymmetric(lowertriangle(st_matrix("est_V_c")))
	mata: st_matrix("est_V_c2", clean_V)
	
	// Step 4

	tempname b V
	
	mat `b' = est_b_c
	mat `V' =  est_V_c2 
	
	local row_eq_names: coleq `b'

	mat colnames `b' = `col_names'
	mat rownames `V' = `col_names'
	mat colnames `V' = `col_names'
	mat coleq `V' = `row_eq_names'
	mat roweq `V' = `row_eq_names'
	
	************ post results *************
	eret post `b' `V', esample(`touse') obs(`N')
	di _newline as text "Multivariate two-stage individual participant data meta-analysis" 
	di _newline as text "Number of studies: `studies'"
	di as text "Number of participants: `N'"
	di as text "Estimation method: `effect'"
	eret display 
	
	************ return results ************* 
		ereturn local prefix twostage
		eret local predict `twostage_predict'
		eret local estat `twostage_estat'
		eret local properties "b V"
		eret local mtype `cmd'
		eret local effect `effect'
		eret local frame_name `frame_name'
		eret local title "Multivariate two-stage individual participant data meta-analysis"
		ereturn local cmdline `"twostage `0': `cmdline'"'
		ereturn local cmdname `cmdname'
		ereturn local cmd twostage
		ereturn local command `"`command'"'
		
		eret scalar studies = `studies'
		eret scalar df_m = `df_m'
end 
