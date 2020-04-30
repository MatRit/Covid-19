clear; close all; clc;
%[] \

% Paramters
Place='Berlin_March';
alpha=0.01:0.01:0.20; %Share of infections needing ICU
Kmax=1:14; % Average time in hospital
L=1:10; % average lag between reported infection and ICU admission
ForecastHorizon=30;

Data=xlsread(['Data_' Place '.xlsx']);

set(0,'defaulttextfontsize',14)
set(0,'defaultaxesfontsize',14)
set(0,'defaultlinelinewidth',2)

Date=datenum(Data(:,1:3));
T=size(Date,1);
Infections=[datenum(Data(:,1:3)) Data(:,4)];
ICUPatients=[datenum(Data(:,1:3)) Data(:,5)];

RelChangeInfections=NaN(T,2);
RelChangeInfections(:,1)=Date;
RelChangeInfections(2:end,2)=(Infections(2:end,2)-Infections(1:end-1,2))./Infections(1:end-1,2);

AbsChangeInfections=NaN(T,2);
AbsChangeInfections(:,1)=Date;
AbsChangeInfections(2:end,2)=Infections(2:end,2)-Infections(1:end-1,2);

Row100=find(Infections(:,2)>=100, 1 );

%% Plot number of infections and ICU patients
figure
plot(Infections(:,1),Infections(:,2))
hold on
plot(ICUPatients(:,1),ICUPatients(:,2))
datetick('x',6,'keeplimits')
legend('Reported infections','ICU patients','Location','Best')
title(Place)


%% Plot absolute and relative change after at least 100 infections
figure
bar(Date(Row100:end,1),AbsChangeInfections(Row100:end,2))
hold on
ylabel('Absolute change of number of infections')
yyaxis right
plot(Date(Row100:end,1),RelChangeInfections(Row100:end,2))
datetick('x',6,'keeplimits')
legend('Absolute change','Relative change','Location','Best')
title(Place)
ylabel('Relative change of number of infections')


%% Fit a model for number of ICU patients
%ModelComparison=NaN;
s=1;
for a=alpha
    for K=Kmax
        for l=L
        EstimatedNewICUPatients=zeros(T,2);
        EstimatedNewICUPatients(:,1)=Date;
        EstimatedICUPatients=zeros(T,2);
        EstimatedICUPatients(:,1)=Date;
        for t=1:T
        if t-max(l+1,K+1)>0 
        EstimatedNewICUPatients(t,2)=a*AbsChangeInfections(t-l,2);
        EstimatedICUPatients(t,2)=sum(EstimatedNewICUPatients(max(1,t-K+1):t,2));
        end
        end
        ModelComparison(s,:)= [a K l sqrt(immse(EstimatedICUPatients(EstimatedICUPatients(:,2)>0,2),ICUPatients(EstimatedICUPatients(:,2)>0,2))) corr(EstimatedICUPatients(EstimatedICUPatients(:,2)>0,2),ICUPatients(EstimatedICUPatients(:,2)>0,2))^2]; %#ok<SAGROW>
        s=s+1;
        end
    end
end

ModelComparison=sortrows(ModelComparison,4);

%% Best model according to RMSE
disp(['Best Model out of ' num2str(size(ModelComparison,1)) ':'  ])
disp(ModelComparison(1,:))

abest=ModelComparison(1,1);
Kbest=ModelComparison(1,2);
lbest=ModelComparison(1,3);

for a=abest
    for K=Kbest
        for l=lbest
        EstimatedNewICUPatients=zeros(T,2);
        EstimatedNewICUPatients(:,1)=Date;
        EstimatedICUPatients=zeros(T,2);
        EstimatedICUPatients(:,1)=Date;
        for t=1:T
        if t-max(l+1,K+1)>0 
        EstimatedNewICUPatients(t,2)=a*AbsChangeInfections(t-l,2);
        EstimatedICUPatients(t,2)=sum(EstimatedNewICUPatients(max(1,t-K+1):t,2));
        end
        end
        end
    end
end

figure
plot(ICUPatients(EstimatedICUPatients(:,2)>0,1),ICUPatients(EstimatedICUPatients(:,2)>0,2))
hold on
plot(EstimatedICUPatients(EstimatedICUPatients(:,2)>0,1),EstimatedICUPatients(EstimatedICUPatients(:,2)>0,2))
datetick('x',6,'keeplimits')
legend('Reported ICU patients',['Predicted ICU patients (\alpha=' num2str(ModelComparison(1,1)) ', K=' num2str(Kbest) ', l=' num2str(lbest) ')' ],'Location','Best')
title(Place)

%% Forecast infections

PredDate=Date(end,1)+1:Date(end,1)+ForecastHorizon;

AvgChangeabs7=(Infections(end,2)-Infections(end-7,2))/7; % Average absolute change in the last 7 days
AvgChangerel7=(Infections(end,2)/Infections(end-7,2))^(1/7); % Averagle relative change in the last 7 days


PredInfectionslin=[Infections; NaN(ForecastHorizon,2)];
PredInfectionslin(T+1:end,1)=PredDate;
PredInfectionsexp=PredInfectionslin; 



for t=T+1:T+ForecastHorizon
PredInfectionslin(t,2)=PredInfectionslin(t-1,2)+AvgChangeabs7;
PredInfectionsexp(t,2)=PredInfectionsexp(t-1,2)*AvgChangerel7;
end

figure
plot(PredInfectionslin(:,1),PredInfectionslin(:,2))
hold on
plot(PredInfectionsexp(:,1),PredInfectionsexp(:,2))
plot(Infections(:,1),Infections(:,2),'k')
datetick('x',6,'keeplimits')
legend(['Linear growth, \Delta= ' num2str(AvgChangeabs7,3)],['Exponential growth, r= ' num2str((AvgChangerel7-1)*100,2)  '%' ],'Location','Best')
title(Place)
ylabel('Predicted infections')

%% Forecast ICU patients

PredNewICUPatientslin=[EstimatedNewICUPatients; NaN(ForecastHorizon,2)];
PredNewICUPatientslin(T+1:end,1)=PredDate;
PredNewICUPatientsexp=PredNewICUPatientslin;

PredICUPatientslin=[ICUPatients; NaN(ForecastHorizon,2)];
PredICUPatientslin(T+1:end,1)=PredDate;
PredICUPatientsexp=PredICUPatientslin;


for a=abest
    for K=Kbest
        for l=lbest
        for t=T+1:T+ForecastHorizon
        PredNewICUPatientslin(t,2)=a*(PredInfectionslin(t-l,2)-PredInfectionslin(t-l-1,2));
        PredICUPatientslin(t,2)=sum(PredNewICUPatientslin(max(1,t-K+1):t,2));
        PredNewICUPatientsexp(t,2)=a*(PredInfectionsexp(t-l,2)-PredInfectionsexp(t-l-1,2));
        PredICUPatientsexp(t,2)=sum(PredNewICUPatientsexp(max(1,t-K+1):t,2));

        end
        end
    end
end

figure
plot(PredICUPatientslin(:,1),PredICUPatientslin(:,2))
hold on
plot(PredICUPatientsexp(:,1),PredICUPatientsexp(:,2))
plot(ICUPatients(:,1),ICUPatients(:,2),'k')
datetick('x',6,'keeplimits')
legend(['Linear growth, \Delta= ' num2str(AvgChangeabs7,3)],['Exponential growth, r= ' num2str((AvgChangerel7-1)*100,2)  '%' ],'Location','Best')
title(Place)
ylabel('Predicted ICU patients')


