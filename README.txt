The code is for designing switching LPV controller under inexact measurement of scheduling parameters. The methodology is explained in the following paper:

P. Zhao and R. Nagamune, “Switching LPV controller design under uncertain scheduling parameters,"Automatica, 2017, 76: 243-250.

The code can also be used for designing (non-switching) LPV controllers when there is no measurement uncertainty. 

File list and function
   1. "SLPV_uncert_main" - main file 
   2. "SLPV_uncert" - creat and solve the LMI problems involved in controller design
   3. "AugPltEval" - obtain plant matrices by plugging value of scheduling parameters  
   4. "AdmRegGrid" - grid the admissible region for solving parameter-dependent LMIs
   5. "LMI_SwSurf" - generate LMIs for switching surface conditions
  
  
CONTACT:
  Please send comments/bug reports to boranzhao9@gmail.com or panzhao@mech.ubc.ca

  Pan Zhao, 
  Control Engineering Lab,
  Department of Mechanical Engineering
  University of British Columbia
  Vancouver, BC, Canada