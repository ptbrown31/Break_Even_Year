global tstep gama sai1 sai2 sai3 deltak optlrsav pd theta1 pbacktime
global Sig0  Sigg Sig E0 Eland Mat0 Mlo0 Mup0
global b11 b12 b21 b22 b23 b32 b33 F2010 F2100 Tat0 Tlo0 
global etha1 etha3 etha4 deltarf Fex L0 Lg0 LA L Ag A A0 pai ro
global R K0 lambda a0 cpricebase actbase alpha

% OUTPUT AND CAPITAL ACCUMULATION
%   T = Total Time horizon
time = 1:T;
tstep = 5; %  (5 years per period)
fosslim = 6000;      %  Maximum cumulative extraction fossil fuels (GtC) /6000/

% OUTPUT AND CAPITAL ACCUMULATION
%   capital share
gama = 0.3;
%   damage coefficient on temperature
sai1 = 0.0;
%   damage coefficient on temperature squared
sai2 = 0.00236; % DICE-2013R; 0.0028388 for DICE-2007
%   Exponent on damages
%sai3 = 2.0;

%   rate of depreciation
deltak = 0.1; 

%   Optimal long-run savings rate used for transversality
optlrsav = (deltak + 0.004)/(deltak + 0.004*elasmu + prstp)*gama;
%   Initial capital stock ($ trillion 2010 USD)
K0 = 223; %

% EMISSIONS

% Industrial emissions 2010 (GtCO2 per year)
e0=35.85;
% Initial emissions control rate for base case 2010
miu0=0.03;
% Initial world gross output (trill 2010 USD)
q0=105.5;
%    Initial Sig
Sig0 = e0/(q0*(1-miu0)); % DICE-2013R
%    Initial growth rate of Sig (percent per decade)
%    Growth rate of Sig (percent per decade)
Sigg(1) = Sigg0;
%    Sig (industrial CO2 emissions/output -- MTC/$1000)
Sig = zeros(1,T);
Sig(1) = Sig0;
for i = 2:T
    Sigg(i) = Sigg(i-1)*((1+Siggd)^tstep);
    Sig(i) = Sig(i-1)*exp(Sigg(i)*tstep);
end

%    Initial carbon emissions from land use change (GTCO2 per year)
E0 = 2.6; % 
%    Decline rate of land emissions (per period)
deland=0.115;
%    Carbon emissions from land use change (GTCO2 per year)
Eland = zeros(1,T);
for i=1:T % DICE-2013R
    Eland(i) = E0*(1-deland)^(time(i)-1);
end

%   Abatement cost function coefficient
theta1 = zeros(1,T+5); % theta1 = (pb*Sig/theta2).*((pr-1+exp(-pd*(time-1)))/pr); % DICE-2007
%   cost of backstop 2010 $ per t C in 2010
pb = 550; 

%   Backstop price
pbacktime = zeros(1,T);
for i=1:T % DICE-2013R
    pbacktime(i) = pb*(1-pd)^(time(i)-1);
    theta1(i) = pbacktime(i)* (Sig(i) )/theta2/1000;
end

% CONCENTRATIONS

%    Initial atmospheric concentration of CO2 (GTC) in 2010
Mat0 = 851;
%    Initial concentration of CO2 in biosphere/shallow oceans (GTC) in 2010
Mlo0 = 1527; 
%    Initial concentration of CO2 in deep oceans (GTC) in 2010
Mup0 = 10010; 
%    Carbon cycle transition coefficients (percent per decade)
%       atmosphere to atmosphere (b11)
b11 = 91.2; 
%       biosphere/shallow oceans to atmosphere (b21)
b21 = 3.83; 
%        atmosphere to biosphere/shallow oceans (b12)
b12 = 8.80; 
%       biosphere/shallow oceans to biosphere/shallow oceans (b22)
b22 = 95.92; 
%        deep oceans to biosphere/shallow oceans (b32)
b32 = 0.03375;
%       biosphere/shallow oceans to deep oceans (b23)
b23 = 0.25; 
%       deep oceans to deep oceans (b33)
b33 = 99.97; 

% TEMPERATURE

%   2010 forcings other ghg
F2010 = 0.5;
%   2100 forcings other ghg
F2100 = 1.0;
%   Initial atmospheric temperature (deg. C above 1900)
Tat0 = 0.85;
%   Initial temperature of deep oceans (deg. C above 1900)
Tlo0 = 0.0068;
%   Speed of adjustment parameter for atmospheric temperature
etha1 = 0.098; 
%   FCO22x Forcings of equilibrium CO2 doubling (Wm-2)
deltarf = 3.6813;
%   Coefficient of heat loss from atmosphere to oceans
etha3 = 0.088; 
%   Coefficient of heat gain by deep oceans
etha4 = 0.025;
%   Exogenous forcing (Watts per square meter)
Fex = zeros(1,T);
for i=1:T % DICE-2013R
    if i<19
        Fex(i) = F2010+ (1/18)*(F2100-F2010)*(time(i)-1);
    else
        Fex(i) = F2100;
    end
end

% POPULATION

%   Initial population (millions) 2010
L0 = 7403;
%   Growth rate to calibrate to 2050 pop projection
Lg0 = 0.134;
%   Population (millions)
L = zeros(1,T);
L(1) = L0;
for i=1:(T-1) % DICE-2013R
    L(i+1) = L(i)*(LA/L(i))^Lg0;
end

% PRODUCTIVITY

%   Initial level of total factor productivity
A0 = 5.115;
%   Rate of growth of productivity (percent per decade)
Ag = zeros(1,T);
%   Total factor productivity
A = zeros(1,T);
A(1) = A0;
for i = 1:(T-1)
    Ag(i)=Ag0*exp(-Agd*tstep*(time(i)-1));
    A(i+1) = A(i)/(1-Ag(i));
end

% PARTICIPATION

%   PARTFRACT (participation rate)
pai = ones(1,T+5);

% WELFARE

%   Social rate of time preference (percent per year)
ro = prstp*100*ones(1,T);
%   Social time preference factor
R = zeros(1,T);
R(1) = 1;
for i = 2:T
    R(i) = R(i-1)/(1+ro(i-1)/100)^tstep;
end
%alpha = 2;
alpha = elasmu;
lambda = 1/(1+ro(1)/100)^tstep;

% Initialization

% Initial Action
a0 = 0.03;

% ** Abatement cost
% Period before which no emissions controls base
tnopol = 45;
% Initial base carbon price (2010$ per tCO2)
cprice0 = 1.0; 
% Growth rate of base carbon price per year
gcprice = 0.02;
% *Base Case    Control Rate
actbase = zeros(1,T);
% *Base Case      Carbon Price
cpricebase = zeros(1,T);
for i = 1:T
    cpricebase(i) = cprice0*(1+gcprice)^(tstep*(time(i)-1));
    actbase(i) = pai(i) * exp( (log(cpricebase(i)/pbacktime(i))) / (theta2-1) );
end