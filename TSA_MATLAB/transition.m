%=========================================================================
%
%   Plot transition functions for the TAR models
%
%=========================================================================

clear all
clc

% Parameters
gamma = 5;
c     = 0.0;


% Take x to be between -1 and 1
x = -2:0.01:2;


G1 = zeros( length(x),1 );
G2 = zeros( length(x),1 );
G3 = zeros( length(x),1 );
G4 = zeros( length(x),1 );

for i = 1:length(x)

        if x(i) >= 0;
            G1(i) = 1.0;
        end
        G2(i) = normcdf( gamma*(x(i)-c) ); 
        G3(i) =1/(1+exp(-gamma*(x(i)-c)));       
        G4(i) =1-exp(-gamma*(x(i)-c)^2);
        
        
        
        
end
        

%**********************************************************************
%***
%***     Generate graph
%***
%**********************************************************************

% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
fig1 = figure(1);
clf;

plot(x,G2,':k',x,G3,'-k',x,G1,'--k',x,G4,'-.k','LineWidth',0.75);
ylabel('$w_t$');
xlabel('$y_{t-d}$')
axis([-2 2 -0.1 1.1 ])
set(gca,'YTick',[0.0 0.2 0.4 0.6 0.8 1.0]);
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'});
set(gca,'XTick',-2:1:2);
legend box off
legend('STAR','LSTAR','SETAR','ESTAR','Location','SouthEast');
legend2latex(fig1); 

% Print the tex file to the relevant directory
laprint(1,'threshold','options','factory');

  
