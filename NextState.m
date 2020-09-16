function [ S2]  =  NextState( S1, a1, t1 )
%Finds the next state given the current state and actions
%   S1: Initial state (length 7) at time t
%   act: action vector
%       act(1): saving rate
%       act(2): GHG reduction rate
%   S2: Next states
%       S2(1:14): Next state

global tstep gama sai1 sai2 sai3 deltak optlrsav theta1 theta2 damage_dep_rate_coeff_on_T damage_TFP_growth_coeff_on_T
global Sig Eland b11 b12 b21 b22 b23 b32 b33 pai a0
global etha1 etha3 etha4 deltarf Fex L A deltaT elasmu Ag
global optimize_for_total_utility abatement_cost_on Burke_damage_on utility_function alpha
global damage_scalar mitigation_scalar


%zero carbon capital DICE config looks like this
    %collect the output
    % S
    % 1 - regular capital X
    % 2 - GMST X
    % 3 - lower ocean temperature X
    % 4 - atmospheric CO2 X
    % 5 - upper ocean CO2 X
    % 6 - lower ocean CO2 X
    % 7 - gross output X
    % 8 - zero-carbon emissions captial (leave NaN)
    % 9 - total capital (same as 1)
    % 10 - abatement cost
    % 11 - epsilon (leave NaN)
    % 12 - Industrial CO2 emissions
    % 13 - CO2 Radiative Forcing
    % 14 - Climate Damage
    % 15 - 
    % 16 - total utility
    % 17 - per capita utility
    % 18 - consumption
    % 19 - consumption per capita
    % 20 - net output
    % 21 - deltaK
    % 22 - TFP growth
    % 23 - TFP    
    % 24 - no damage gross output

S2 = zeros(24, 1);

if t1==1
    a1= a0;
end

% stop negative emissions when GMST reaches zero

if S1(2) < 0
    
    a1 = 1;
    
end

%Industrial Emission
IE = Sig(t1) * (1 - a1) * S1(7);

S2(12) = IE;

%Total Emission
E = IE + Eland(t1);

%CO2 doubling coefficient
etha2 = deltarf / deltaT;

%Lower Ocean Concentration
S2(6) = b23 / 100 * S1(5) + b33 / 100 * S1(6);

%Upper Ocean Concentration
S2(5) = b12 / 100 * S1(4) + b22 / 100 * S1(5)+ b32 / 100 * S1(6);

%Atmospheric Concentration
S2(4) = (5 / 3.666) * E + b11 / 100 * S1(4) + b21 / 100 * S1(5);

%Lower Ocean Temperature
S2(3) = S1(3) + etha4 * (S1(2) - S1(3));

%Radiative Forcings
F = deltarf * ((log(S2(4)) - log(588)) / log(2)) + Fex(t1 + 1);

S2(13) = F;

%Atmospheric Temperature
S2(2) = S1(2) + etha1 * (F - etha2 * S1(2) - etha3 * (S1(2) - S1(3)));

%Climate Damage
    %default DICE:
    
    Damage = 1;
    
        if Burke_damage_on == 0

            %Damage = 1 - (sai1 * S1(2) + sai2 * S1(2) ^ sai3);  % 
            Damage = damage_scalar*(1 - (sai1 * S1(2) + sai2 * S1(2) ^ sai3));  % 
            
            S2(14) = Damage;
            
        end
    
% Damage on depcreciation rate

if Burke_damage_on == 1

      deltak_damaged = 0.108.*S1(2) - damage_dep_rate_coeff_on_T.*S1(2).^2;
      if deltak_damaged < deltak; deltak_damaged = deltak; end

      S2(21) = deltak_damaged;
        
      deltak_undamaged = deltak;
    
    % Damage on TFP growth rate
    
        Ag_damaged = Ag(t1) - damage_TFP_growth_coeff_on_T*S1(2);
                    
        A_damaged = S1(23)/(1-Ag_damaged);  %S1(23) is the previous timesteps damaged TFP

        S2(23) = A_damaged;
        
        A_undamaged = A(t1+1);
    
end

%Abatement cost = AbatedEmission * Carbonprice / theta2

if abatement_cost_on == 1
    
    %Abate = pai(t1) ^ (1 - theta2) * theta1(t1) * a1 ^ theta2;
    Abate = mitigation_scalar*(pai(t1) ^ (1 - theta2) * theta1(t1) * a1 ^ theta2);
    
end

if abatement_cost_on == 0

    Abate = 0;
    
end

%store abatement cost
S2(10) = Abate;

%Economic output net of damage and abatement costs
Q = (Damage - Abate) * S1(7);

S2(20) = Q;

if Burke_damage_on == 1

    %Capital at time t1+1
    S2(1) = tstep * optlrsav * Q + (1 - deltak_damaged) ^ tstep * S1(1);

    %Gross output at time t1+1
    S2(7) = A_damaged * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);
    
    % ------ %
    
    %calculate fraction of damage
    Q_undamaged = (1 - Abate) * S1(24);
    
    capital_undamaged = tstep * optlrsav * Q_undamaged + (1 - deltak_undamaged) ^ tstep * S1(1);
    gross_output_undamaged = A_undamaged * (capital_undamaged ^ gama) * (L(t1 + 1)/1000)^(1 - gama);
    
    %damaged gross output
    S2(7) = A_damaged * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);
    
    %fraction of damaged gross output to imagninary undamaged output
    Damage = S2(7)./gross_output_undamaged;
    
    S2(14) = Damage;
    
    %Gross output undamaged at time t1+1
    
    S2(24) = A_undamaged * (capital_undamaged ^ gama) * (L(t1 + 1)/1000)^(1 - gama);
    
end
        
if Burke_damage_on == 0 || Burke_damage_on == 2

    %Capital at time t1+1
    S2(1) = tstep * optlrsav * Q + (1 - deltak) ^ tstep * S1(1);

    %Gross output at time t1+1
    S2(7) = A(t1+1) * (S2(1) ^ gama) * (L(t1 + 1)/1000)^(1 - gama);
    
end


% per capital consumption and utility

    %   Consumption
    C = (1 - optlrsav) * Q;
    
    S2(18) = C;

    %   Consumption per capita
    c = C / L(t1+1) * 1000;
    
    S2(19) = c;

    %   Utility per capita
    
    if utility_function == 1
        u = (c^(1-elasmu)-1)/(1-elasmu)-1;
    end

    if utility_function == 0
        u = (c ^ (1 - alpha)) / (1 - alpha);
    end
    
    S2(17) = u;
    
    %   Social utility
        if optimize_for_total_utility == 1

            U = u * L(t1+1);
            
        end
    
        if optimize_for_total_utility == 0

            U = u;

        end

    S2(16) = U;

end

