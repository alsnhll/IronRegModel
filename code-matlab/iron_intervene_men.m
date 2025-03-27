% SYNOPSIS: Simulates change in iron levels with various interventions
%           to reduce menstruation for various initial conditions of chronic
%           anemia (depletions in Hb levels and iron stores)
%
% INPUT:    Current code contains parameters for Indian females, with a weight
%           of 55 kg, and median healthy Hb of 13 g/dL and 0.7 g of non-Hb
%           iron. Anthropomorphic information is used to convert [Hb] to total
%           grams of iron in Hb.
%
% OUTPUT:   Plot 2 figures, 1) timecourse of change in [Hb], and 2) timecourse 
%           of other body iron 
%
% Other functions called: 
%           ironsolve.m     contains differential equations
%           ode45.m         (MATLAB function) numerically integrates equations
%           absp.m          calculates the absorption rate
%           eryth.m         calculates the erythropoeisis rate
%
% Written by Alison Hill, alhill@fas.harvard.edu, last updated Sept 21 2010
% Modified by ET, March 20 2025

function iron_intervene_men

%%%%%%%%%%%%%%%%% Input parameters %%%%%%%%%%%%%%%%%%%%%%%%
weight = 55;                % [kg] female weight
e1 = 0.00060;               % [g/day] baseline daily menstrual excretion, as long as body Fe > 0. (0.001)
e2 = 0.00106;               % [g/day] baseline daily other excretion, as long as body Fe > 0. (0.001)
Tend = 30;                  % [months] time to run simulation
hb0 = 13;                   % [g/dL] "Healthy" Hb levels
dep_hb = [0, 20, 50];       % Percent depletion for iron in Hb
red_men = [0, 20, 50];      % Percent reduction in Fe lost to menstruation

%%%%%%%%%%%%%%%%% Conversion %%%%%%%%%%%%%%%%%%%%%%%%
PV=weight*0.2*0.2;          % healthy plasma volume
BV=PV/(1-0.38);             % blood volume, 0.38 is healthy hematocrit
conv=285/(10*BV);           % 285 is conversion of g Fe to g Hb, BV is blood volume

%%%%%%%%%%%%%%%%% Steady-State Scenario %%%%%%%%%%%%%%%%%%%%%%%%
d=0.0055;                   % [/day] rate of Hb turnover, goes back to body iron,=1/death rate=ln(2)/half life, half life=127 days
h0=(d*(hb0/conv)+e1)/0.7;   % 0.013; % [/day] erythropoiesis rate (of body iron going to HB iron). 
                            % Calculated to get s.s. with normal [Fe]'s

% Expected steady state iron stores for a given chronic anemic [Hb]
OBI0 = (((1-dep_hb/100)*hb0/conv)*d+e1)./eryth((1-dep_hb/100)*hb0/conv,h0,conv); % Fe in OBI based on SS anemia value
Int = (e1+e2)./absp((1-dep_hb/100)*hb0/conv,conv);     % [g/day] daily intake to maintain SS anemia value

%%%%%%%%%%%%%%%%% Intervention of Menstruation Reduction %%%%%%%%%%%%%%%%%%%%%%%%
results = cell(length(dep_hb), length(red_men));

for j = 1:length(dep_hb)
    for k=1:length(red_men)
        [results{j, k}.T, results{j, k}.Y] = ode45(@(t, x) ironsolve(t, x, Int(j), e1, e2*(1-red_men(k)/100), d, h0, conv), [0 Tend]*30, [OBI0(j); (1-dep_hb(j)/100)*13/conv]);
    end
end

%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(dep_hb)
    fprintf('Scenario %d:\n', i);
    fprintf('Percent depletion of Iron from Healthy Hb: %.1f %%\n', dep_hb(i));
    fprintf('Fe in OBI: %.4f g\n', OBI0(i));
    fprintf('S.S. Intake values: %.1f mg/day\n', Int(i)*1000);
    fprintf('\n');
end

gcf = figure(1);
set(gcf, 'DefaultAxesColorOrder', [
    0.2, 0.4, 0.6;      % Light blue
    0.6, 0.8, 0.2;      % Light green
    0.8, 0.1, 0.1;      % Dark red
    0.5, 0, 0.5;        % Purple
    1, 0.5, 0;          % Orange
    1, 0.75, 0.8;       % Pink
    0, 0.5, 0.5;        % Teal
    0.6, 0.3, 0.1;      % Brown
    0.5, 0.5, 0.5;      % Gray
    ])

for i=1:length(dep_hb)
    for j=1:length(red_men)
        plot(results{i, j}.T/30, results{i, j}.Y(:, 2)*conv,'DisplayName', sprintf('%.1f%% anemic, %.0f%% reduction in menstruation', dep_hb(i), red_men(j)))
        hold on
    end
    hold on
end
hold off

xlabel('Time (months)')
ylabel('Hemoglobin (g/dL)')
legend('Location', 'best')
ylim([6 14])

gcf = figure(2);
set(gcf, 'DefaultAxesColorOrder', [
    0.2, 0.4, 0.6;      % Light blue
    0.6, 0.8, 0.2;      % Light green
    0.8, 0.1, 0.1;      % Dark red
    0.5, 0, 0.5;        % Purple
    1, 0.5, 0;          % Orange
    1, 0.75, 0.8;       % Pink
    0, 0.5, 0.5;        % Teal
    0.6, 0.3, 0.1;      % Brown
    0.5, 0.5, 0.5;      % Gray
    ])

for i=1:length(dep_hb)
    for j=1:length(red_men)
        plot(results{i, j}.T/30, results{i, j}.Y(:, 1),'DisplayName', sprintf('%.1f%% anemic, %.0f%% reduction in menstruation', dep_hb(i), red_men(j)))
        hold on
    end
 
    hold on
end
hold off

xlabel('Time (months)')
ylabel('Other Body Fe (g)')
legend('Location', 'best')
ylim([0 1])

end