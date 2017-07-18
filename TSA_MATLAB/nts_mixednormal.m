%=========================================================================
%
%   Program to simulate a mixed normal distribution and 
%   compare it to the standard normal distribution
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',1234) )

t    = 4000;
reps = 20000;
b    = zeros(reps,1);

for i = 1:reps

  % Independent random numbers used for y2 and u
  y2 = cumsum(randn(t,1));           
  u  = randn(t,1);

  % regress y2 on u and scale the slope estimate by t
  b(i) = t*(y2\u);                     

end

disp(['Mean     = ',num2str(mean(b) ) ]); 
disp(['Std.dev. = ',num2str(std(b)) ]);
disp(['Skewmess = ',num2str(mean((b-mean(b)).^3)./std(b)^3) ]);
disp(['Kurtosis = ',num2str(mean((b-mean(b)).^4)./std(b)^4) ]);

minx = -10; 
maxx = 10;
x    = seqa( minx,(maxx-minx)/200,201 );

ftrue = normpdf( (x-mean(b))/std(b) )/std(b);
fhat  = ksdensity( b,x );

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;

%--------------------------------------------------------%
plot(x,ftrue,'-k', ...
     x,fhat,'--k');
box off
ylim([0 0.3]);
xlim( [-10 10]);
set(gca,'xtick',[-10 -5 0 5 10]);
set(gca,'ytick',[0.0 0.1 0.2 0.3]);
xlabel('$m$');
ylabel('$f(m)$');

% Print the tex file to the relevant directory
%laprint(1,'mixednormal','options','factory');





