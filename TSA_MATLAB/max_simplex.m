%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	Compute the steps of simplex and find the step thatis used 
%%  	in first iteration.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ch3_7
clear all;
clc;

% Read in the data for y
y = [ 3.5; 1.0; 1.5; ];
t = length( y );       	% Define the sample size

th  = [1 ; 3 ];     % Initialize the estimates

neglog = lnl_neg(th,y,t);

if(neglog(1) > neglog(2))
    flipud (th);
    flipud(neglog);    
end

th_avg  = mean(th);

%   The three steps that may be followed

%   Reflect
alpha = 0.5;
th_r  = th_avg + alpha*(th_avg - th(2)); 

%   Expand
beta = 1.1;
th_e = th_avg + beta*(th_avg - th(2));

%   Contract
gamma = 0.5;
th_c = th_avg + gamma*(th_avg - th(2));

% Decide which step is followed

llen_r = lnl_neg(th_r,y,t);

if(llen_r < neglog(2))
    
    if(llen_r < neglog(1))
        fprintf('Expanding the simplex \n');
        llen_e = lnl_neg(th_e,y,t);
        if (llen_e < llen_r)
            fprintf('Theta2 is replaced by theta_e \n');
            th(2) = th_e;
        else
            fprintf('Theta2 is replaced by theta_r \n');
            th(2) = th_r;
        end  
    else
        fprintf('Reflecting the simplex \n');
        fprintf('Theta2 is replaced by theta_r \n');
        th(2) = th_r;
    end  
else
    llen_c = lnl_neg(th_c,y,t);
    if(llen_c < neglog(2))
        fprintf('Contracting the simplex \n');
        fprintf('Theta2 is replaced by theta_c \n');
        th(2) = th_c;
    else
         fprintf('Shrinking the simplex \n');
         th = (th + th_avg) / 2;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Negative log likelihood evaluator

function lnl_neg = lnl_neg(theta,y,t)
lnl_neg = t*log(theta) + sum(y)./theta;
end


