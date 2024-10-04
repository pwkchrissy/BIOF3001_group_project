# BIOF3001_group_project
Codes for illustrating the estimation of reproductive numbers in R and MATLAB

To estimate the R_t  values in our example in Hong Kong, we would first need to deconvolute the time series of cases by dates of symptom onset to obtain the time series of cases by dates of infection, for all the reported cases in Hong Kong. 

The file “2023_09_23_mpox_CHP_press_release.xlsx” stored the symptom onset dates per reported patient in the 2023 mpox outbreak in Hong Kong. The Matlab file “main_deconvol.m” reads the former file and deconvolute the data of onset time. It could be done through the self-defined function stored in the file “DeconvolutionIncidence1.m”, which provides the time series of the infection dates. 

Following the above steps, we could then estimate the values of R_t  over time from the times series of infection dates using the R file “main_mpox_hk.R”. The output of the R codes is the mean and SD of the serial interval and that of the R_t  estimates. Please read the codes in the files mentioned above to have a more detailed understanding towards the steps for estimations of R_t  . Please also add annotations to your codes as shown in the example. 

To facilitate your understanding, we further provided an R version for the MATLAB codes for deconvolution (DeconvolutionIncidence1.R for DeconvolutionIncidence1.m, main_deconvol.R for main_deconvol.m). You might work on the project based on either the R version or the MATLAB version. 

