%=========================================================================
%
%   Program to reproduce Diebold and Yilmaz (2009) spillover model
%=========================================================================
clear all
clc

% Read the data
% stock market returns and volatility (4-Dec-1996 to 23-Nov-2007)
% Order of countries
%                         1.  US
%                         2.  UK
%                         3.  France
%                         4.  Germany
%                         5.  Hong Kong
%                         6.  Japan
%                         7.  Australia
%                         8.  Indonesia
%                         9.  S. Korea
%                         10. Malaysia
%                         11. Philippines
%                         12. Singapore
%                         13. Taiwan
%                         14. Thailand
%                         15. Argentina
%                         16. Brazil
%                         17. Chile
%                         18. Mexico
%                         19. Turkey

load('diebold_yilmaz.mat') 

% Choose the data
y = returns;        
z = bsxfun(@rdivide,bsxfun(@minus,y,mean(y)),std(y));

disp('    Mean      Median    Max       Min       Std.dev.  Skew      Kurt')
disp( [ mean(y)' median(y)' max(y)' min(y)' std(y)' mean(z.^3)' mean(z.^4)'] )

% Estimate a VAR with lag p=2  
x = [ones(length(y)-2,1)   trimr(y,1,1)   trimr(y,0,2) ];
b = x\trimr(y,2,0);

mu  = trimr(b,0,38)';  % Vector of intercepts                                              
phi1 = trimr(b,1,19)'; % Lag 1 parameter estimates   
phi2 = trimr(b,20,0)'; % Lag 2 parameter estimates     
  
% Generate VMA (non-orthogonalized) for horizons 1 to 10 
si1  = eye(size(y,2));
si2  = phi1;
si3  = phi1*si2 + phi2;
si4  = phi1*si3 + phi2*si2;
si5  = phi1*si4 + phi2*si3;
si6  = phi1*si5 + phi2*si4;
si7  = phi1*si6 + phi2*si5;
si8  = phi1*si7 + phi2*si6;
si9  = phi1*si8 + phi2*si7;
si10  = phi1*si9 + phi2*si8;

% Generate VMA (orthogonalized) for horizons 1 to 10  
v  = trimr(y,2,0) - x*b;   % VAR residuals      
vc = v'*v/length(v);
d  = diag(vc);
s  = chol(vc)';

ir1  = si1*s;
ir2  = si2*s;
ir3  = si3*s;
ir4  = si4*s;
ir5  = si5*s;
ir6  = si6*s;
ir7  = si7*s;
ir8  = si8*s;
ir9  = si9*s;
ir10 = si10*s;

% Compute variance decompositions for horizons 1 to 10    

vd1  = ir1.^2;
vd1  = 100*bsxfun(@rdivide,vd1,sum(vd1,2));

vd2  = ir1.^2 + ir2.^2;
vd2  = 100*bsxfun(@rdivide,vd2,sum(vd2,2));

vd3  = ir1.^2 + ir2.^2 + ir3.^2;
vd3  = 100*bsxfun(@rdivide,vd3,sum(vd3'));

vd4  = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2;
vd4  = 100*bsxfun(@rdivide,vd4,sum(vd4'));

vd5  = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2 + ir5.^2;
vd5  = 100*bsxfun(@rdivide,vd5,sum(vd5'));

vd6  = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2 + ir5.^2 + ir6.^2;
vd6  = 100*bsxfun(@rdivide,vd6,sum(vd6'));

vd7  = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2 + ir5.^2 + ir6.^2 + ir7.^2;
vd7  = 100*bsxfun(@rdivide,vd7,sum(vd7'));

vd8  = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2 + ir5.^2 + ir6.^2 + ir7.^2 + ir8.^2;
vd8  = 100*bsxfun(@rdivide,vd8,sum(vd8'));

vd9  = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2 + ir5.^2 + ir6.^2 + ir7.^2 + ir8.^2 + ir9.^2;
vd9  = 100*bsxfun(@rdivide,vd9,sum(vd9'));

vd10 = ir1.^2 + ir2.^2 + ir3.^2 + ir4.^2 + ir5.^2 + ir6.^2 + ir7.^2 + ir8.^2 + ir9.^2 + ir10.^2;
vd10 = 100*bsxfun(@rdivide,vd10,sum(vd10,2));

str = [  'US     '  ;    
         'UK     '  ;    
         'FRA    '  ;
         'GER    '  ;
         'HKG    '  ;
         'JPN    '  ;
         'AUS    '  ;
         'IDN    '  ;    
         'KOR    '  ;    
         'MYS    '  ;    
         'PHL    '  ;    
         'SGP    '  ;    
         'TAI    '  ;    
         'THA    '  ;    
         'ARG    '  ;    
         'BRA    '  ;    
         'CHL    '  ;    
         'MEX    '  ;    
         'TUR    '   ];

disp('Variance decomposition at period 10')
disp(vd10); 

tmp=sum(vd10,2)-diag(vd10);
disp('Contribution From Others')
disp( [str num2str(tmp) ]);

tmp=(sum(vd10)-diag(vd10)')';
disp('Contribution To Others')
disp( [str num2str(tmp) ]);

