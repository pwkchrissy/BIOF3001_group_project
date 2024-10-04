
clear all

% assume that the mean and SD of the incubation period are 7.6 and 5.8
% respectively (gamma distributed)
incubationMean = 7.6;
incubationSD = 5.8;

% obtain the scale and shape parameters of incubation period
gammaScale = incubationSD*incubationSD/incubationMean;
gammaShape = incubationMean/gammaScale;

% obtain the convolution "pdf" (in terms of days) 
pdfConvolFG  = gamcdf(1:31,gammaShape,gammaScale)-gamcdf(0:30,gammaShape,gammaScale);
pdfConvolFG = [pdfConvolFG,1-sum(pdfConvolFG)];

% reformat the data for later deconvolution
dateZero = datenum('2023-06-01','yyyy-mm-dd');
mpoxLinelist = readtable('2023_09_23_mpox_CHP_press_release.xlsx'); % read the file
mpoxLinelist.date_symptom_onset_rev = datenum(mpoxLinelist.date_symptom_onset)-dateZero; % convert to day 1, 2, 3 etc
caseReportHistc = histcounts(mpoxLinelist.date_symptom_onset_rev,0:(max(mpoxLinelist.date_symptom_onset_rev)+1));
caseReportCount = array2table([1:(max(mpoxLinelist.date_symptom_onset_rev)+1);caseReportHistc]',"VariableNames",{'date','new_cases_onset'});
writetable(caseReportCount,'output/mpox_epi_curve_onset.csv');

% obtain the time series of the infection dates of all cases
caseReport = caseReportCount.new_cases_onset;
caseReport(caseReport==0) = 1e-9;
deconvolInfection2 = DeconvolutionIncidence1(caseReport,pdfConvolFG);

% plot both the time series of infection dates and reported dates 
clf;
figure(2)
plot(caseReport, color = 'red')
hold on
plot(deconvolInfection2(1:(end-5)), color = 'blue')

writetable(array2table([...
    deconvolInfection2],...
    'VariableNames',{'all_cases'}),'output/deconvolInfection.csv');