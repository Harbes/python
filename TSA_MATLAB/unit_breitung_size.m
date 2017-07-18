%=========================================================================
%
%   Monte Carlo experiment on the Breitung test
%
%=========================================================================
clear all
clc

reps = 10000;
Tv   = [100,200];
arv  = [0,0.3,0.6,0.9];
mav  = [-0.8,-0.4,0,0.4,0.8];
rej  = zeros(length(arv),length(mav));
rho  = zeros(reps,1);

disp( '     T          arv      mav       Rej Freq ')
disp('---------------------------------------------')
for i = 1:length(Tv)
    
  for j = 1:length(arv)
      
    for m = 1:length(mav)
        
      RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )
      
      for k =1:reps
          
          v      = randn(Tv(i)+1,1);
          y      = cumsum(recserar(trimr(v,1,0)+mav(m)*trimr(v,0,1),v(1),arv(j)));
          rho(k) = sum(cumsum(y).^2)/(sum(y.^2)*(length(y)^2));
          
      end  
      
      rej(j,m) = mean(rho < 0.02);
      disp( [Tv(i) arv(j) mav(m) rej(j,m)])  
      
    end  
    
  end
  
end

    


