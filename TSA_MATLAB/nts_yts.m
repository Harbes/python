%=========================================================================
%
%   Program to construct the step function of Y[ts]
%
%=========================================================================
clear all;
clc;

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123) )

% Generate simulated values for t = 5
t  = 5;                              
y0 = 0;
vt5 = randn( t,1 );
yt5 = zeros( t,1 );

for j = 1:t
        
    yt5(j) = y0 + sum( vt5(1:j) );
    
end

yts5 = zeros( 100,1);

k = 1;
for s = 0.01:0.01:1
    
    p = floor(s*t);
    if p == 0;
        yts5(k) = y0;
    else
        yts5(k) = y0 + sum( vt5(1:p) );
    end
    
    k =k+1;
end

% Generate simulated values for t = 40
t  = 40;                              
y0 = 0;
vt = randn( t,1 );
yt = zeros( t,1 );

for j = 1:t
        
    yt(j) = y0 + sum( vt(1:j) );
    
end

yts = zeros( 100,1);

k = 1;
for s = 0.01:0.01:1
    
    p = floor(s*t);
    if p == 0;
        yts(k) = y0;
    else
        yts(k) = y0 + sum( vt(1:p) );
    end
    
    k =k+1;
end


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
% Panel (a)
subplot(2,2,1)
plot(1:1:5,yt5,'.k');
title('(a) Discrete Representation T=5');
ylabel('$y_t$' )
xlabel('t')
box off
%ylim([1 3 ]);
xlim( [0 5]);
set(gca,'xtick',[0 1 2 3 4 5]);
     
%--------------------------------------------------------%
% Panel (b)
subplot(2,2,2)
plot(1:1:40,yt,'.k');
title('(b) Discrete Representation T=40');
ylabel('$y_t$' )
xlabel('t')
box off
xlim( [0 40]);
set(gca,'xtick',[0 10 20 30 40]);
hold off
  
%--------------------------------------------------------%
% Panel (c)
subplot(2,2,3)
plot(0.01:0.01:1,yts5,'-k');
title('(c) Continuous Representation T=5');
ylabel('$y_{[Ts]}$' )
xlabel('t')
box off
%ylim([1 3])
xlim( [0 1]);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);

%--------------------------------------------------------%
% Panel (d)
subplot(2,2,4)
plot(0.01:0.01:1,yts,'-k');
title('(d) Continuous Representation T=40');
ylabel('$y_{[Ts]}$' )
xlabel('t')
box off
xlim( [0 1]);
set(gca,'xtick',[0 0.2 0.4 0.6 0.8 1]);
 
% % Print the tex file to the relevant directory
%laprint(1,'yts','options','factory');
 


