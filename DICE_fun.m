function [ welfare ] = DICE_fun( act )

global S lambda T S0

% ** Scaling and inessential parameters
%         scale1      Multiplicative scaling coefficient               /0.016408662 /
%         scale2      Additive scaling coefficient                     /-3855.106895/ ;
scale1 = 0.016408662;
scale2 = -3855.106895;
tstep = 5;

S = zeros(24, T);
S(:, 1) = S0;
welfare = 0;

for t = 1:1:(T-1)
    util = Utility(S(:, t), act(t), t);
    welfare = welfare + lambda^(t-1) * util;
    S(:, t + 1) = NextState(S(:, t), act(t), t);
end

welfare = welfare*scale1*tstep + scale2;

welfare = -1 *welfare;

end

