%=========================================================================
%
%   Approximate asymptotic critical values for PP coefficient test
%
%=========================================================================
clear all
clc

reps = 100000;
t    = 1000;
s    = zeros(reps,1);

for j=1:reps;

  b    = cumsum(randn(t,1))/sqrt(t);
  s(j) = 0.5*(b(t)^2-1)/mean(trimr(b,0,1).^2);

end

disp(' Critical Values')
disp('      1%        5%        10% ')
disp(quantile(s,[0.01 0.05 0.1]));
