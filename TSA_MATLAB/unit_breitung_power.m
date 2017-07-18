%=========================================================================
%
%   Empirical power of the Breitung test
%
%=========================================================================
clear all
clc

% Parameters
t    = 1000;
cv   = [0,-5,-10,-15,-20];
reps = 100000;
tdf  = zeros(reps,length(cv));
rho  = zeros(reps,length(cv));

for i = 1:length(cv)
    
  RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) );

  for j = 1:reps 
      
    y        = recserar(randn(t,1),randn(1,1),1+cv(i)/t);
    phihat   = trimr(y,0,1)\trimr(y,1,0);
    vhat     = trimr(y,1,0)-phihat*trimr(y,0,1);
    tdf(j,i) = (phihat-1)/sqrt((vhat'*vhat/length(vhat))/(trimr(y,0,1)'*trimr(y,0,1)));
    rho(j,i) = sum(cumsum(y).^2)/(sum(y.^2)*(t^2));
    
  end
end

disp('    c          tdf        rho') 
disp('-----------------------------------------------')

disp([ cv' mean(tdf<quantile(tdf(:,1),0.05))' mean(rho < quantile(rho(:,1),0.05))' ]);

