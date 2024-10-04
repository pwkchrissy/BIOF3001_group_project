rm(list = ls());

# assume that the mean and SD of the incubation period are 7.6 and 5.8
# respectively (gamma distributed)
incubationMean = 7.6;
incubationSD = 5.8;

# obtain the scale and shape parameters of incubation period
gammaScale = incubationSD*incubationSD/incubationMean;
gammaShape = incubationMean/gammaScale;

# obtain the convolution "pdf" (in terms of days) 
pdfConvolFG  = pgamma(1:31,shape=gammaShape,scale=gammaScale)-pgamma(0:30,shape=gammaShape,scale=gammaScale);
pdfConvolFG = c(pdfConvolFG,1-sum(pdfConvolFG));

# reformat the data for later deconvolution
dateZero = as.Date('2023-06-01',format = '%Y-%m-%d');
mpoxLinelist = read.xlsx('2023_09_23_mpox_CHP_press_release.xlsx'); # read the file
mpoxLinelist$date_symptom_onset_rev = convertToDate(mpoxLinelist$date_symptom_onset, "1899-12-30")-dateZero; # convert to day 1, 2, 3 etc
mpoxLinelist$date_symptom_onset_rev = as.numeric(mpoxLinelist$date_symptom_onset_rev); # convert to numeric form

days_sequence <- seq(0, max(mpoxLinelist$date_symptom_onset_rev)+1);
counts_per_day <- table(factor(mpoxLinelist$date_symptom_onset_rev, levels = days_sequence));
counts_per_day <- as.data.frame(counts_per_day);
# fine tune to make the date one day earlier
counts_per_day$Var1 = 1:(max(mpoxLinelist$date_symptom_onset_rev)+2);
counts_per_day <- counts_per_day[-nrow(counts_per_day),];

colnames(counts_per_day) = c('date','new_cases_onset')
write.csv(counts_per_day,'output/mpox_epi_curve_onset_R.csv',row.names = FALSE);

# obtain the time series of the infection dates of all cases
caseReport = counts_per_day$new_cases_onset;
caseReport[caseReport==0] = 1e-9; # to avoid calculation errors

source("DeconvolutionIncidence1.R")
deconvolInfection2 = DeconvolutionIncidence1(caseReport,pdfConvolFG);
deconvolInfection2 = as.data.frame(deconvolInfection2);
colnames(deconvolInfection2)= c('all_cases');

# plot both the time series of infection dates and reported dates 
library(ggplot2)

caseReport_table = as.data.frame(cbind(1:length(caseReport), caseReport));
colnames(caseReport_table) = c('date', 'number');
deconvol_table = as.data.frame(cbind(1:length(caseReport), deconvolInfection2));
deconvol_table = head(deconvol_table, -5);
colnames(deconvol_table) = c('date', 'number');
  
ggplot()+ geom_line(data= caseReport_table, aes(x=date, y=number),color = "red")+
  geom_line(data= deconvol_table, aes(x=date, y=number), color = 'blue')+theme_classic();

# save the result 
write.csv(deconvolInfection2, 'output/deconvolInfection_R.csv',,row.names = FALSE)


