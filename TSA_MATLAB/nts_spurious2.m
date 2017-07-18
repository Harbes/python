%=========================================================================
%
%   Program to demonstrate the spurious regression problem 
%   using least squares regression
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) )

t = 100;                            
r = 10000;                          



% Generate sampling distributions and plot the histograms
% I(0) versus I(0)      
b = zeros(r,2);

for i = 1:r
    
    y1 = trimr(recserar(randn(t+100,1),0.0,0.0),100,0);
    y2 = trimr(recserar(randn(t+100,1),0.0,0.0),100,0);
        
    b(i,:) = [ones(t,1) y1]\y2;

end
figure(1)
hist(b(:,2),21);
title(' I(0) vs. I(0)' );

% I(1) versus I(1)       
for i = 1:r
    
    y1 = trimr(recserar(randn(t+100,1),0.0,1.0),100,0);
    y2 = trimr(recserar(randn(t+100,1),0.0,1.0),100,0);
    
    b(i,:) = [ones(t,1) y1]\y2;
end

figure(2)
hist(b(:,2),21);
title(' I(1) vs. I(1) ' );

% I(1) versus I(2)       
for i = 1:r 
    
    y1 = trimr(recserar(randn(t+100,1),0.0,1.0),100,0);
    y2 = trimr(recserar(randn(t+100,1),[0.0;0.0],[2.0;-1.0]),100,0);
        
    b(i,:) = [ones(t,1) y1]\y2;
end
figure(3)
hist(b(:,2),21);
title(' I(1) vs. I(2) ' );

% I(2) versus I(2)       
for i = 1:r
    
    y1 = trimr(recserar(randn(t+100,1),[0.0;0.0],[2.0;-1.0]),100,0);
    y2 = trimr(recserar(randn(t+100,1),[0.0;0.0],[2.0;-1.0]),100,0);

    b(i,:) = [ones(t,1) y1]\y2;
end
figure(4)
hist(b(:,2),21);
title(' I(2) vs. I(2) ' );

