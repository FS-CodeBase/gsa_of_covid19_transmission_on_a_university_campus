Global Sensitivity Analysis (GSA) of Covid-19 Transmission on a University Campus
** All code is written in MATLAB R2021b

-----------
Directories
-----------
1)	contact_matrices/
> Folder contains the contact matrices for both campus contacts and dorm living contacts under different classroom and dorm capacity scenarios. The following key indicates the specifics of each contact scenario based on a case number 1 through 24.
	                               Class Caps
Number of Students|  None	100 students	50 students	25 students
      in Dorms    |------------------------------------------------
              All | 24	    23	            22	        21
	          2.5K| 20	    19	            18	        17
              2.0K|	16	    15	            14	        13
              1.5K|	12	    11	            10	        9
              1.0K|	8	    7	            6	        5
              0.5K|	4	    3	            2	        1

2)	model_sols/
> Model solutions created by the script compute_all_model_solutions.m are saved in this folder.

3)	sobol_indices_dt/
> Sobol indices computed using the doubling time metric, with the script compute_sobol_indices_doubling_time.m, are saved to this folder.

4)	sobol_indices_in_time/
> Sobol indices computed in time, using the cumulative number of cases computed with the script compute_sobol_indices_in_time.m are saved to this folder.

-------
Scripts
-------
1)	compute_all_model_solutions.m
> Computes all of the model solutions for particular contact scenarios (see contact_matrices) and model parameter ranges (see fun_model_parameter_ranges.m). All solutions are saved in the model_sols/ folder.

2)	compute_sobol_indices_doubling_time.m
> Computes the first and total-effect Sobol indices using the disease doubling time metric (see fun_doubling_time.m)
and model solutions found in the model_sols/ folder. The sobol indices are saved in the sobol_indices_dt/ folder.

3)	compute_sobol_indices_in_time.m
> Computes the first and total-effect Sobol indices using the cumulative number of cases (see fun_sobol_indices_by_solution_in_time.m) and the model solutions found in the model_sols/ folder. The resulting sobol indices are saved in the sobol_indices_in_time/ folder.

---------
Functions
---------
1)	covid_model_ode.m
4)	> SEIR model of disease transmission dynamics on a university campus. This function is used by the compute_all_model_solutions.m script.

2)	fun_doubling_time.m
> Computes the disease doubling time over a desired number of consecutive days since the start of a semester.

3)	fun_initialize_contact_matrices.m
> Combines contact matrices found in the contact_matrices/ folder into arrays that can be accessed by case number (see Directories section). 

4)	fun_lhs_sampling_of_parameters.m
> Performs Latin hypercube sampling (lhs) of desired parameters in the ranges specified in the fun_model_parameter_ranges.m function. Model solutions are then computed using the sets of sampled parameter values.

5)	fun_model_parameter_ranges.m
> Stores the user indicated model parameters ranges and has a binary option for whether a parameter is to be varied (true) or not (false).

6)	fun_sobol_indices_by_measure.m
> Computes the first and total-effect Sobol indices using a measure function that is applied to model solutions stored in the model_solutions/ folder.

7)	fun_sobol_indices_by_solution_in_time.m
> Computes the first and total-effect Sobol indices using the raw model solutions saved in the model_solutions/ folder.  
