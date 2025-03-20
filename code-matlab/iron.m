function iron


% SYNOPSIS: Master function to simulate the differential equations and plot
%           the results for a two compartment model of iron regulation, including 
%           increased erythropoiesis and intestinal absorption in response to Hb levels
%
% INPUT:      Current code contains parameters for Indian females, with a weight
%       of 55 kg, and median healthy Hb of 13 g/dL and 0.7 g of non-Hb
%       iron. Anthropomorphic information is used to convert [Hb] to total
%       grams of iron in Hb.
%
%       Initial conditions are currently set to represent a 10% blood loss,
%       for example, a blood donation. 
%
% OUTPUT: Plot 4 figures, 1) timecourse of change in [Hb] and other body
%       iron (OBI), 2) timecourse of erythropoiesis absorption rate, 3) dose
%       response curve for iron absorption versus [Hb], 4) dose response curve
%       for erythropoiesis rate versus [Hb]
%
% Other functions called: 
%           ironsolve.m     contains differential equations
%           ode45.m         (MATLAB function) numerically integrates equations
%           absp.m          calculates the absorption rate
%           eryth.m         calculates the erythropoeisis rate
%
% Written by Alison Hill, alhill@fas.harvard.edu, last updated Sept 21 2010


weight=55; %female weight in kg
PV=weight*0.2*0.2; % healthy plasma volume
BV=PV/(1-0.38); %blood volume, 0.38 is healthy hematocrit
conv=285/(10*BV) %285 is conversion of g Fe to g Hb, BV is blood volume


e1=0.00106; %(g/day) daily menstrual excretion, as long as body Fe > 0. (0.001)
e2=0.00060; %(g/day) daily other excretion, as long as body Fe > 0. (0.001)
d=0.0055; %0.0055(/day) rate of Hb turnover, goes back to body iron,=1/death rate=ln(2)/half life, half life=127 days
h0=(d*(13/conv)+e1)/0.7  %0.013; %(/day) erythropoiesis rate (of body iron going to HB iron). 
                            %Calculated to get s.s. with normal [Fe]'s

L=55*10^(-3); % (g/day) daily intake.

Tend=10; %time to run simulation for,  in months

I=[0.7 0.9*13/conv]; %(g) initial iron concentrations, 0.9 is for a 10% blood donation

[T,Y] = ode45(@(t,x)ironsolve(t,x,L,e1,e2,d,h0,conv),[0 Tend]*30,I); % run the simulation

figure(1)

subplot(1,2,1);

[AX,H1,H2]=plotyy(T/30,Y(:,2)*conv,T/30,Y(:,1));

set(get(AX(2),'Ylabel'),'String','Other body Fe (g)---  ')
set(get(AX(1),'Ylabel'),'String','Hemoglobin (g.dL ^{-1})')
set(AX(2),'Ylim',[0.6 0.70])
set(AX(2),'YTick',[0.6:0.05:0.70])
set(AX(1),'Ylim',[11 13])
set(AX(1),'YTick',[11:1:13])
xlabel('time (months)')

set(H1,'Color','k')
set(H2,'Color','k')

set(AX(1),'YColor','k')
set(AX(2),'YColor','k')

subplot(1,2,2)
[AX,H1,H2]=plotyy(T/30,100*absp(Y(:,2),conv),T/30,1/h0*eryth(Y(:,2),h0,conv));
set(get(AX(1),'Ylabel'),'String','Fe absorption (%)')
set(get(AX(2),'Ylabel'),'String','erythropoeisis rate (X baseline) (.d^{-1})---  ')
set(AX(1),'Ylim',[0 8])
set(AX(1),'YTick',[0:2:8])
set(AX(2),'Ylim',[1 3])
set(AX(2),'YTick',[1:0.5:3])
xlabel('time (months)')

set(H1,'Color','k')
set(H2,'Color','k')

set(AX(1),'YColor','k')
set(AX(2),'YColor','k')

figure(2)
x=0:0.01:4;
plot(x*conv,100*absp(x,conv),13,100*absp(13/conv,conv),'.r')
xlabel('Hemoglobin (g/dL)')
ylabel('Fe absorption (%)')
legend('response curve','normal value')
xlim([0 20])

figure(3)
plot(x*conv,1/h0*eryth(x,h0,conv),13,1,'.r')
xlabel('Hemoglobin (g/dL)')
ylabel('erythropoesis rate (X baseline) (/day)')
legend('rate','normal')
xlim([0 20])


end
