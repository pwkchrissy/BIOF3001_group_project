# 2024-09-30
# MPOX HK
# Estimation of Rt

rm(list = ls());

library(EpiEstim)
library(ggplot2)
library(MASS)
library(openxlsx)

dateZero = as.Date('2023-06-01',format = '%Y-%m-%d');

# read the excel file into R
hkLineList = read.xlsx('2023_09_23_mpox_CHP_press_release.xlsx');
hkLineList$date_symptom_onset_rev = hkLineList$date_symptom_onset-as.numeric(dateZero);

# convert the dates into date format instead of numerical format (optional)
for (col in c("date_press_release","date_clinic_visit","date_symptom_onset","date_exposure_start","date_exposure_end")) {
  hkLineList[[col]] <- format(as.Date(hkLineList[[col]], origin = "1899-12-30"), "%d-%m-%Y")
}

# UK study: https://www.bmj.com/content/379/bmj-2022-073153
# Incubation period: 7.6 (95% CrI 6.5-9.9) or 7.8 (6.5-9.2)
# Serial interval: 8.0 (95% CrI 6.5-9.8) or 9.5 (7.4-12.3)

# read the epi curve for mpox (obtained using Matlab codes)
hkEpiCurve = read.csv('output/deconvolInfection.csv',header = TRUE);
hkEpiCurve$all_cases = hkEpiCurve$all_cases;
hkEpiCurve$date = dateZero+(1:length(hkEpiCurve$all_cases));

# obtain the minimum and maximum dates of the epi curve
minDate = min(hkEpiCurve$date,na.rm = TRUE);
maxDate = max(hkEpiCurve$date,na.rm = TRUE);

# adjust the epi curve for future analysis (make it integers)
hkEpiCurve$all_cases[hkEpiCurve$all_cases<0.2] =  0; 
hkEpiCurve$all_cases[hkEpiCurve$all_cases>=0.2&hkEpiCurve$all_cases<1] =  1;
rec = data.frame(dates = maxDate+(-(length(hkEpiCurve$date)-1):0),
                 all_cases = round(hkEpiCurve$all_cases));
minDate = min(rec$dates[rec$all_cases>=1]);
rec = rec[rec$dates>=minDate,];

# 6. Estimate Rt
# call function
## fixing the random seeds
MCMC_seed  = 1;
overall_seed  = 2;
mcmc_control = make_mcmc_control(seed = MCMC_seed,burnin = 1000,thin = 10);
# fitting a Gamma distribution for the SI
dist = "G";
maxT = length(rec$dates)
t_start = seq(2,(maxT-6));
t_end = t_start+6;
config <- make_config(list(mean_si = 9.7, 
                           std_si = 11.2,
                           mcmc_control = mcmc_control,
                           seed = overall_seed, 
                           t_start = t_start,
                           t_end=t_end,
                           n1 = 1000, 
                           n2 = 100));

# Excluding imported cases, start estimation
inputIncRev = rec[,c('dates','all_cases')];
colnames(inputIncRev) = c('dates','I');
minDate = rec$dates[1];
res_si_from_data <- estimate_R(inputIncRev,
                               method = "parametric_si",
                               config = config);

# record Rt estimated 
tempRec = res_si_from_data;
outputRec = tempRec$R;
# add dates
outputRec$t_start = outputRec$t_start+minDate-1;
outputRec$t_end = outputRec$t_end+minDate-1;
# create the file storing the Rt estimates
write.xlsx(outputRec,paste0('output/mpox_Rt_local_transmissibility.xlsx'),sheetName = 'Sheet1',rowNames = FALSE);
# do the same for the serial interval estimates
write.xlsx(tempRec$SI.Moments,'output/mpox_Rt_serial_interval.xlsx',sheetName = 'Sheet1',rowNames = FALSE)

