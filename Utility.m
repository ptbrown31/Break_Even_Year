%Algorithm for calculating the social utility of each state under ceratin
%action (a3)

function [ U ]  =  Utility( S4, a4, t4 )
%Social Utility for a given state St
%   state = St
%   action = a
%   time = t2
global optlrsav sai1 sai2 sai3 theta1 theta2 L elasmu pai a0
global optimize_for_total_utility abatement_cost_on Burke_damage_on alpha utility_function
global damage_scalar mitigation_scalar

if t4==1
    a4=a0;
end

%   Damage cost function
    % Previous Burke damage emulation
    
    Damage = 1;
    
    if Burke_damage_on == 0
    
        %Damage = 1-(sai1 * S4(2) + sai2 * (S4(2) ^ sai3)); % DICE-2013R
        Damage = damage_scalar*(1-(sai1 * S4(2) + sai2 * (S4(2) ^ sai3))); % DICE-2013R
        
    end
        
% the TFP and depreciation rate damage pathways are already included in gross output S4(7)

% Abatement cost = AbatedEmission * Carbonprice / theta2

if abatement_cost_on == 1
    
    %Abate = (pai(t4) ^ (1 - theta2)) * theta1(t4) * (a4 ^ theta2);
    Abate = mitigation_scalar*((pai(t4) ^ (1 - theta2)) * theta1(t4) * (a4 ^ theta2));
    
end

if abatement_cost_on == 0
    Abate = 0;
end

%   Net output after damage and abatement
Q = (Damage -Abate) * S4(7);

%   Consumption
C = (1 - optlrsav) * Q;

%   Consumption per capita
c = C / L(t4) * 1000;

%   Utility per capita
% u = 1 + c ^ (1 - alpha) / (1 - alpha); This is used in DICE-2007

if utility_function == 1
    u = (c^(1-elasmu)-1)/(1-elasmu)-1;
end

if utility_function == 0
    u = (c ^ (1 - alpha)) / (1 - alpha);
end

if utility_function == 2 %maximize per capita consumption
    u = c;
end

%what if we are just maximizing consumption per capita
%u = c;

%   Social utility

    if optimize_for_total_utility == 1;

        U = u * L(t4);

    end

    if optimize_for_total_utility == 0;

        U = u;

    end
end

