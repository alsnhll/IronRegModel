% SYNOPSIS: simulates change in iron levels with various intakes
%           ('interventions') with various initial conditions of chronic
%           anemia (depletions in Hb levels and iron stores)
%
% INPUT:      Current code contains parameters for Indian females, with a weight
%       of 55 kg, and median healthy Hb of 13 g/dL and 0.7 g of non-Hb
%       iron. Anthropomorphic information is used to convert [Hb] to total
%       grams of iron in Hb.
%
%
% OUTPUT: Plot 4 figures, 1) timecourse of change in [Hb], 2) timecourse of
%           iron absorption rate and 3) timecourse of other body iron, 4)
%           erythropoiesis rate
%
%
% Other functions called: 
%           ironsolve.m     contains differential equations
%           ode45.m         (MATLAB function) numerically integrates equations
%           absp.m          calculates the absorption rate
%           eryth.m         calculates the erythropoeisis rate
%
% Difference from iron_intervene_2 : has depletion of body stores
%
% Written by Alison Hill, alhill@fas.harvard.edu, last updated Sept 21 2010

function iron_intervene_2

weight=55; %female weight in kg
PV=weight*0.2*0.2; % healthy plasma volume
BV=PV/(1-0.38); %blood volume, 0.38 is healthy hematocrit
conv=285/(10*BV); %285 is conversion of g Fe to g Hb, BV is blood volume

e1=0.00106; %(g/day) daily menstrual excretion, as long as body Fe > 0. (0.001)
e2=0.00060; %(g/day) daily other excretion, as long as body Fe > 0. (0.001)
d=0.0055; %0.0055(/day) rate of Hb turnover, goes back to body iron,=1/death rate=ln(2)/half life, half life=127 days
h0=(d*(13/conv)+e1)/0.7;  %0.013; %(/day) erythropoiesis rate (of body iron going to HB iron). 
                            %Calculated to get s.s. with normal [Fe]'s

Tend=30; %time to run simulation for,  in months

L_range=[10+10 20+10 50+10]*10^(-3); % (g/day) daily intake.

%I (g) initial iron concentrations 
% to get expected steady state iron stores for a given chronic anemic [Hb],
% use differential equations
L=L_range(1); [T1,Y1] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,[((0.5*13/conv)*d+e1)./eryth(0.5*13/conv,h0,conv); 0.5*13/conv]);
L=L_range(2); [T2,Y2] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,[((0.5*13/conv)*d+e1)./eryth(0.5*13/conv,h0,conv); 0.5*13/conv]);
L=L_range(3); [T3,Y3] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,[((0.5*13/conv)*d+e1)./eryth(0.5*13/conv,h0,conv); 0.5*13/conv]);

L=L_range(1); [T5,Y5] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,[((0.8*13/conv)*d+e1)./eryth(0.5*13/conv,h0,conv); 0.8*13/conv]);
L=L_range(2); [T6,Y6] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,[((0.8*13/conv)*d+e1)./eryth(0.5*13/conv,h0,conv); 0.8*13/conv]);
L=L_range(3); [T7,Y7] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,[((0.8*13/conv)*d+e1)./eryth(0.5*13/conv,h0,conv); 0.8*13/conv]);


gcf=figure(7);
set(gcf,'DefaultAxesColorOrder',jet(3))

subplot(2,1,1)
plot(T1/30,Y1(:,2)*conv,T2/30,Y2(:,2)*conv,T3/30,Y3(:,2)*conv,T5/30,Y5(:,2)*conv,'--',T6/30,Y6(:,2)*conv,'--',T7/30,Y7(:,2)*conv,'--')
xlabel('time(months)')
ylabel('Hemoglobin (g/dL)')
legend(sprintf('initially 50%% depleted \n 10 mg/day'),'20 mg/day','50 mg/day',sprintf('initially 20%% depleted \n 10 mg/day'),'20 mg/day','50 mg/day')
%legend('5 mg/day,50% loss','10 mg/day,50% loss','20 mg/day,50% loss','50 mg/day,50% loss','5 mg/day,20% loss','9/day,20% loss','20 mg/day,20% loss','50 mg/day,20% loss')
%ylim([11 15])

%{
subplot(2,2,2)
plot(T1/30,100*absp(Y1(:,2),conv),T2/30,100*absp(Y2(:,2),conv),T3/30,100*absp(Y3(:,2),conv),T4/30,100*absp(Y4(:,2),conv),T5/30,100*absp(Y5(:,2),conv),'--',T6/30,100*absp(Y6(:,2),conv),'--',T7/30,100*absp(Y7(:,2),conv),'--',T8/30,100*absp(Y8(:,2),conv),'--')
xlabel('time (months)')
ylabel('Fe absorption (%)')
%legend('5 mg/day,50% loss','9/day,50% loss','20 mg/day,50% loss','50 mg/day,50% loss','5 mg/day,20% loss','9/day,20% loss','20 mg/day,20% loss','50 mg/day,20% loss')
% xlim([0 Tend])
% ylim([2 10])
%}
subplot(2,1,2)
plot(T1/30,Y1(:,1),T2/30,Y2(:,1),T3/30,Y3(:,1),T5/30,Y5(:,1),'--',T6/30,Y6(:,1),'--',T7/30,Y7(:,1),'--')
xlabel('time(months)')
ylabel('Other body Fe (g)')

%legend('5 mg/day,50% loss','9/day,50% loss','20 mg/day,50% loss','50 mg/day,50% loss','5 mg/day,20% loss','9/day,20% loss','20 mg/day,20% loss','50 mg/day,20% loss')
%ylim([11 15])
%{

subplot(2,2,4)
plot(T1/30,1/h0*eryth(Y1(:,2),h0,conv),T2/30,1/h0*eryth(Y2(:,2),h0,conv),T3/30,1/h0*eryth(Y3(:,2),h0,conv),T4/30,1/h0*eryth(Y4(:,2),h0,conv),T5/30,1/h0*eryth(Y5(:,2),h0,conv),'--',T6/30,1/h0*eryth(Y6(:,2),h0,conv),'--',T7/30,1/h0*eryth(Y7(:,2),h0,conv),'--',T8/30,1/h0*eryth(Y8(:,2),h0,conv),'--')
xlabel('time (months)')
ylabel('erythropoesis rate (X baseline) (/day)')
legend(sprintf('initially 50%% depleted \n 5 mg/day'),'10 mg/day','20 mg/day','50 mg/day',sprintf('initially 20%% depleted \n 5 mg/day'),'10 mg/day','20 mg/day','50 mg/day')
%}


end
