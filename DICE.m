% clc
% clear;
Initialset;
% Initial year 2010; time step 5 yr

global  S S0 K0 Tat0 Tlo0 Mat0 Mlo0 Mup0 pbacktime pai theta2
global  L A T gama fval aopt deltak
global optimize_on

% Initialization
S = zeros(24, T);
y0 = A(1)*(K0^gama)*(L(1)/1000)^(1-gama);
S0 = [K0, Tat0, Tlo0, Mat0, Mlo0, Mup0, y0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, deltak, NaN, A(1), y0];

myoptions = optimset('Display','iter','FunValCheck','on','algorithm','sqp','MaxFunEvals',20000);
act0 = zeros(T, 1)+0.5;
actlower = zeros(T,1);
actupper = ones(T,1);

actlower(1)=0;
actupper(1)=0.001;

actupper(2)=0.039;

% % Upper limit on control rate for negative emissions after 2150          / 1.2  /
% actupper(1:5,1) = linspace(0,1.27,5);
% actupper(6:T,1) = actupper(6:T,1).*1.27;

if optimize_on == 1
    
    [aopt, fval] = fmincon(@DICE_fun,act0,[],[],[],[],actlower,actupper, [], myoptions);
    
end

if optimize_on == 0
    
    fval = DICE_fun(aopt);

end