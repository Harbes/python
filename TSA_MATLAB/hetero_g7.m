% ========================================================================
%
%   Program to estimate heteroskedastic models of the business cycle 
%
% ========================================================================
function hetero_g7(  )

    % Set country
    country = 1;    % 1 for Canada
                    % 2 for France
                    % 3 for Germany
                    % 4 for Italy
                    % 5 for Japan
                    % 6 for UK
                    % 7 for US

    
    % Get the data 
    data = getData( country );
    y  = data(:,1);                 % Growth rate of GDP
    x  = data(:,2);                 % Spread
     

    % Mean variables: lagged gdp growth and lagged interest spread
    x = trimr([y  x],0,1);       

    % Variance variables: lagged gdp growth and lagged interest spread
    w = x;
    
    y = trimr(y,1,0);
    t = length(y);

    % Estimate unrestriced model
    start = 0.1*ones(6,1);   
    ops   = optimset('LargeScale', 'off', 'Display', 'final');
    
    [bhat1,lf1] =  fminunc(@(p) neglog(p,y,x,w),start,ops);  
    lf1 = -lf1;

    % Estimate restriced model
    start = 0.1*ones(4,1);
        
    [bhat0,lf0] =  fminunc(@(p) neglog0(p,y,x),start,ops);  
    lf0 = -lf0;

    % LR test
    lr  = -2*t*(lf0 - lf1);
    dof = length(bhat1) - length(bhat0);

    disp(' ')
    disp(['Log-likelihood function (unrestricted) = ',num2str(lf1) ]);
    disp(['Log-likelihood function (restricted)   = ',num2str(lf0) ]);
    disp(['LR statistic                           = ',num2str(lr) ]);
    disp(['p-value                                = ',num2str(1-chi2cdf(lr,dof)) ]);

end

%
%--------------------------- Functions -----------------------------------
% 
%-------------------------------------------------------------------------
% Unrestricted likelihood function 
%-------------------------------------------------------------------------
function lf = neglog(p,y,x,w)

    m  = p(1) + p(2)*x(:,1) + p(3)*x(:,2);
    s2 = exp(p(4) + p(5)*w(:,1) + p(6)*w(:,2));           
    lt = -0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m).^2)./s2;        

	lf = - mean( lt );

end
%-------------------------------------------------------------------------
% Restricted likelihood function 
%-------------------------------------------------------------------------
function lf = neglog0(p,y,x)

    m  = p(1) + p(2)*x(:,1) + p(3)*x(:,2);
    s2 = exp(p(4));           
    lt = -0.5*log(2*pi) - 0.5*log(s2) - 0.5*((y - m).^2)./s2;        

	lf = - mean( lt );

end
%-------------------------------------------------------------------------
% Gets the data for a specfic country and change dates in specific cases
%-------------------------------------------------------------------------
function y = getData( flag )

    if flag == 1        % Canada
    
        y = xlsread('G7Data.xlsx','B8:C162');
        y = trimr(y,0,19); % Change Canada dates: June 1961 to December 1999     

    elseif flag == 2    % France    
    
        y = xlsread('G7Data.xlsx','D43:E157');
    
    elseif flag == 3    % Germany
    
        y = xlsread('G7Data.xlsx','F4:G162');
    
    elseif flag == 4    % Italy
    
        y = xlsread('G7Data.xlsx','H48:I162');
    
    elseif flag == 5    % Japan
    
        y = xlsread('G7Data.xlsx','J43:K160');
    
    elseif flag == 6    % UK
    
        y = xlsread('G7Data.xlsx','L4:M162');
        y = trimr(y,0,20); % Change UK dates: June 1960 to December 1999 
    
    else
        
        y = xlsread('G7Data.xlsx','N4:O162');
        y = trimr(y,0,20);  % Change US dates: June 1960 to December 1999
    
    end
end
