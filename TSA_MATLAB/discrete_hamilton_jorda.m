%=========================================================================
%
%   Program to estimate the Hamilton and Jorda (2002)
%   ordered probit model of US monetary policy.
%
%=========================================================================
function discrete_hamilton_jorda( )

    clear all
    clc

    % Read the data: US weekly yields   
    % 1st week of February 1984 to last week of April 2001
    load hamilton_jorda.mat
    
    event   = data(:,1);
    target  = data(:,2);
    change  = data(:,3);
    bin     = data(:,4);
    spread6 = data(:,5);


	for i = 2:length(event)
		if event(i,1) == 1
			change(i,1) = change(i,1);
        else
			change(i,1) = change(i-1,1);
        end
    end
        
    % Convert data to event time
    capt      = length(event);
    eventdata = zeros(capt,2);
    ybin      = zeros(capt,1);

    lagchange = 0;
    i = 1;
    for t = 1:capt
        if event(t,1) == 1;
            eventdata(i,1) =  lagchange;
            eventdata(i,2) = spread6(t,1);  
            ybin(i,1)      = bin(t,1);

                if t > 1
                    lagchange = target(t,1) - target(t-1,1);  

                end
                i = i+1;
        end
    end
    summ = sum(event);
    ybin = ybin(1:summ); 
    x    = eventdata(1:summ,:);

    % Create dummy variables for each interest rate change
    d1 = double(ybin == -0.50);
    d2 = double(ybin == -0.25);
    d3 = double(ybin ==  0.00);
    d4 = double(ybin ==  0.25);
    d5 = double(ybin ==  0.50);

    d  = [ d1 d2 d3 d4 d5];
        
    % Choose event days from 1984 to 1997
    x = x(1:102,:);         
    d = d(1:102,:);

    % Estimate the ordered probit model
    ops    = optimset('LargeScale','off','Display','iter');
    theta0 = [
                2.5449149
                0.54142729
               -1.8948826
               -0.42001991 
               -0.0052515480
                1.5173916 
              ];

    [theta1,l1,~,~,~,h] = fminunc(@(b) lprobit(b,x,d),theta0,ops);
    
    disp(['Unrestricted log-likelihood function =     ',num2str(l1) ]);
    disp(['T x unrestricted log-likelihood function = ',num2str(t*l1)]);

    disp(' Parameter estimates ');
    disp( theta1 );

end
%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
%  Unrestricted Probit negative log-likelihood function 
%-------------------------------------------------------------------------
function lf = lprobit(b,x,d)

    % Cut off points
    c = b(3:6);
    
    % Regression part excluding the intercepts
    xb = x*b(1:2);
    
    % Cut off points
    f1 = normcdf( c(1) - xb,0,1);
	f2 = normcdf( c(2) - xb,0,1) - f1;
	f3 = normcdf( c(3) - xb,0,1) - f1 - f2;
	f4 = normcdf( c(4) - xb,0,1) - f1 - f2 - f3;
	f5 = 1                       - f1 - f2 - f3 - f4;
    f  = [ f1  f2  f3  f4  f5 ];   
    
    % Log-likelihood function
    tp = bsxfun(@times,d,log(f));
    lf = -mean( sum(tp,2) );     

end 
