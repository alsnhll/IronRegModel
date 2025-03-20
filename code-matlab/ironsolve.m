
%SYNOPSIS: Differential equations describing a two compartment
%       model of iron regulation, including increased erythropoiesis and
%       intestinal absorption in response to Hb levels
%
% Called by: ode45(@(t,x)ironsolve(t,x,L,e0,d,h0,conv) in iron.m, iron_let.m, 
%           iron_ss.m, iron_intake.m, iron_intervene.m
%
% INPUT:     t  current time (not used in current equations but could be
%               changed), independant variable
%            y  dependent variable, vector of iron concentrations
%            L  iron intake, g/day
%            e0 iron excretion, g/day
%            d  death rate of RBCs, /day
%            h0 baseline rate of erythropoiesis, /day 
%            conv conversion between plasma iron (g) and [Hb] (g/day),
%            [plasma Fe]=[Hb]/conv
%
% OUTPUT:   dy, vector of derivatives (ie dy/dt)
%
% Other functions called: 
%           absp.m          calculates the absorption rate
%           eryth.m         calculates the erythropoeisis rate
%           exc.m           step function so that excretion is constant
%           until iron levels go to zero, then excretion goes to zero
%
% Written by Alison Hill, alhill@fas.harvard.edu, last updated Sept 21 2010


function dy = ironsolve(t,y,L,e1,e2,d,h0,conv)
        dy = zeros(2,1);    % a column vector
        
        dy(1) = L*absp(y(2),conv)-exc(y(1),e2)-eryth(y(2),h0,conv)*y(1)+d*y(2);
        
        dy(2) = eryth(y(2),h0,conv)*y(1)-d*y(2)-exc(y(2),e1);
        
    end