# VAR_Project

Written in Matlab. 

This code constructs a bivariate VAR using UK data on unemployment 
and GDP per capita. (See GDP PC.xls and UR.xls for raw data or 
project_data.xlsx for the combined speadsheet) 

The long run identification scheme found in Blanchard and Quah (1989)
is followed to determine the effects of supply and demand shocks on the 
economy. 

The baseline model uses the window 1971-2019 with 4 lags, a VAR(4). 

-------------------------------------------------------------------
Download project_code.m for the code.

The three related functions: 
inpaint_nans.m
boundedline.m
companion.m
(ShadePlotForEmphasis.m is not currently used).

The data are in excel:
Main sheet - project_data.xlsx

(Raw data from ONS: GDP PC.xls + UR.xls)

