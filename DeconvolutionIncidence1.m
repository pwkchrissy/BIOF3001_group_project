%%%To deconvolute a hospitalization incidence curve into a case incidence
%%%curve
%%%Input 1: daily hospitalization number (HospitalIncidence)
%%%Input 2: PDF for hospitalization delay (starting from 0 day, 1 day
%%%increment) (HospitalDelay)
%%%Output: Deconvoluted incidence curve (starting from TypicalHospitalDelay
%%% days shifted back from the start of the HospitalIncidence curve)

function IncidenceCurve = DeconvolutionIncidence1(HospitalIncidence, HospitalDelay) 

%%%%%%%%%%%%%%


%sizes
Temp = size(HospitalIncidence);
NumOfDays = max(Temp);
if (NumOfDays == 0)
    disp('Error in deconvolution: Empty hospitalization incidence vector');
    pause;
end
Temp = size(HospitalDelay);
NumOfDelay = max(Temp);
if (NumOfDelay == 0)
    disp('Error in deconvolution: Empty hospitalization delay vector');
    pause;
end
  
HospitalDelayCDF = HospitalDelay;
for i = 1 : length(HospitalDelay)
    HospitalDelayCDF(i) = sum(HospitalDelay(1:i));
end

%time to hospitalization with highest prob
TypicalHospitalDelay = find (HospitalDelay == max(HospitalDelay));

%%%The initial guess for infection incidence (up to a multiple)
IncidenceCurve = HospitalIncidence;

%%%transition probability
p = zeros(NumOfDays);
for j = 1 : NumOfDays
    for i = j : NumOfDays
        if (i-j < NumOfDelay)
           p(i,j) = HospitalDelay(i-j+1);
        end
    end
end

%%%observation prob
q = zeros(NumOfDays,1);
for j = 1: NumOfDays
    q(j) = sum(p(:,j));
    if ((q(j)>1) || (q(j)==0))
        disp('Error in deconvolution: q_j');
        q(j)
       
        pause;
    end
end

p_hat = zeros(NumOfDays);
for j = 1: NumOfDays
    for i = j : NumOfDays
        p_hat(i,j) = p(i,j)/q(j);
    end
end

%%%initialize the intermediate parameters
DelayTemp = zeros(NumOfDays,1);

%%%routine
for k = 1 : 10
    for i = 1 : NumOfDays
    DelayTemp(i) = dot(p_hat(i,1:NumOfDays),IncidenceCurve');
    end

    for j = 2 : NumOfDays
        itemp = 0;
        for i = 2 : NumOfDays
            itemp = itemp + (p(i,j)*HospitalIncidence(i)/DelayTemp(i));
        end
    IncidenceCurve(j) = IncidenceCurve(j)/q(j)*itemp;
    end 
end

for i = 0 : length(HospitalDelay) -1
    IncidenceCurve(end-i) = IncidenceCurve(end-i)/HospitalDelayCDF(i+1);
end

return