%=========================================================================
%
%   Program to identify properties of US macro variables
%
%=========================================================================
function stsm_roots( )

    clear all
    clc
%%
    % Read the data: quarterly US data from Jan-1959 to Dec-1998
    load simsdata.mat
 
    % Define varaibles
    r    = ytdata(:,1);        
    lex  = log( ytdata(:,2) );
    lcp  = log( ytdata(:,3) );
    lm   = log( ytdata(:,4) );
    lp   = log( ytdata(:,5) );
    lo   = log( ytdata(:,6) );
    sdum = ytdata(:,7:17);
    
    % Construct variables for use in VAR: from 1960:1 to 1998:12    
    % interest rate and the annual % growth rates in money, price and output
    yvar = [r  lm  lp  lo ];             
    y    = [trimr(yvar(:,1),12,0)   100*(trimr(yvar(:,2:4),12,0) - trimr(yvar(:,2:4),0,12)) ];

    % Compute the roots of the AR(4) model for the interest rate
    y     = r;
    X     = [ones(length(y)-4,1) trimr(y,3,1) trimr(y,2,2) trimr(y,1,3) trimr(y,0,4)];    
    theta = X\trimr(y,4,0); 
    
    c  = [-theta(5); -theta(4); -theta(3); -theta(2); 1 ];
    rt = roots( c );
 %%   
    format short
    disp('AR(4) model for interest rates')
    disp('------------------------------')
    disp('  Roots               Absolute Value of Roots')
    disp( [ rt abs( rt ) ] )
    disp('  ')
    
    clear c rt;
    
    % Compute the roots of the AR(2) model for real output growth
    y     = 100*(trimr(lo,12,0) - trimr(lo,0,12));
    theta = [ones(length(y)-2,1) trimr(y,1,1) trimr(y,0,2) ]\trimr(y,2,0);
    
    c = [ -theta(3); -theta(2); 1 ];
    rt = roots( c );
    
    format short
    disp('AR(2) model for real output growth')
    disp('----------------------------------')
    disp('  Roots       Absolute Value of Roots')
    disp( [ rt abs( rt ) ] )
    disp('  ')

    clear c rt;
        
    % Compute the roots of the AR(2) model for log output 
    y     = log( ytdata(:,6) );
    theta = [ones(length(y)-2,1) trimr(y,1,1) trimr(y,0,2) ]\trimr(y,2,0);
    
    c = [ -theta(3); -theta(2); 1 ];
    rt = roots( c );
    
    format short
    disp('AR(2) model for log of real output')
    disp('----------------------------------')
    disp('  Roots       Absolute Value of Roots')
    disp( [ rt abs( rt ) ] )
    disp('  ')
    
    clear c rt;
    
    % Compute the roots of the VAR(2) model  
    yvar  = [ lm lp ];
    y     = trimr(yvar,12,0) - trimr(yvar,0,12);
    X     = [ones(length(y)-2,1) trimr(y,1,1) trimr(y,0,2)];
    theta = X\trimr(y,2,0); 
        
    mu   = trimr(theta,0,4)';
    phi1 = trimr(theta,1,2)';
    phi2 = trimr(theta,3,0)';
      
    [~,rt] = polyeig(eye(2),-phi1,-phi2);
    
    format short
    disp('VAR(2) model of money and inflation')
    disp('----------------------------------')
    disp('  Roots       Absolute Value of Roots')
    disp( [ rt abs( rt ) ] )
    disp('  ')
    
    clear rt
    
    % Compute the roots of the VAR(2) model containing  
    % interest rate, money growth rate, inflation and real GDP growth rate    
    yvar  = [ r  lm  lp  lo ];  
    y     = [trimr(yvar(:,1),12,0)  100*(trimr(yvar(:,2:4),12,0) - trimr(yvar(:,2:4),0,12)) ];
    X     = [ ones(length(y)-2,1) trimr(y,1,1) trimr(y,0,2) ];      
    theta = X\trimr(y,2,0);

    mu   = trimr(theta,0,8)';                                            
    phi1 = trimr(theta,1,4)';                                               
    phi2  = trimr(theta,5,0)';     
    
    [~,rt] = polyeig(eye(4),-phi1,-phi2);

    format short
    disp('VAR(2) model of interest rate, money, inflation and real gdp')
    disp('------------------------------------------------------------')
    disp('  Roots       Absolute Value of Roots')
    disp( [ rt abs( rt ) ] )
    disp('  ')

    clear rt
    
    % Compute the roots of the VAR(4) model containing  
    % interest rate, money growth rate, inflation and real GDP growth rate    
    yvar  = [ r  lm  lp  lo ];  
    y     = [trimr(yvar(:,1),12,0)  100*(trimr(yvar(:,2:4),12,0) - trimr(yvar(:,2:4),0,12)) ];
    X     = [ ones(length(y)-4,1) trimr(y,3,1) trimr(y,2,2) trimr(y,1,3) trimr(y,0,4)];      
    theta = X\trimr(y,4,0);

    mu   = trimr(theta,0,16)';                                                 
    phi1 = trimr(theta,1,12)';                                                 
    phi2 = trimr(theta,5,8)';                                                 
    phi3 = trimr(theta,9,4)';                                                  
    phi4 = trimr(theta,13,0)';                                                 
   
    [~,rt] = polyeig(eye(4),-phi1,-phi2,-phi3,-phi4);

    format short
    disp('VAR(2) model of interest rate, money, inflation and real gdp')
    disp('------------------------------------------------------------')
    disp('  Roots       Absolute Value of Roots')
    disp( [ rt abs( rt ) ] )
    disp('  ')

end

