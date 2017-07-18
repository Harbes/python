%=========================================================================
%
%   Program to generate simulated time series plots for 
%   alternative specifications of the vecm
%
%=========================================================================
clear all
clc

RandStream.setDefaultStream( RandStream('mt19937ar','seed',123457) );

t = 200;            

experiment = 4;                                                  

if experiment == 1

    beta  = [1;-1];                 % Normalized cointegrating vector                              
    beta0 = 2;                      % Long-run constant                                            
    beta1 = 0;                      % Long-run time trend                                          

    alpha = [-0.1;0.1];             % Error-correction parameters                                  

    alpha0 = null(alpha*alpha');
    delta0 = alpha0*0;              % Short-run intercepts 
    delta1 = alpha0*0;              % Short-run time trend using orthogonal complement of alpha    

elseif experiment == 2

    beta = [1;-1];                  % Normalized cointegrating vector                              
    beta0 = 2;                      % Long-run constant                                            
    beta1 = 0;                      % Long-run time trend                                          

    alpha = [-0.0;0.01];           % Error-correction parameters                                  

    alpha0 = null(alpha*alpha');
    delta0 = alpha0*0;              % Short-run intercepts using orthogonal complement of alpha    
    delta1 = alpha0*0;              % Short-run time trend using orthogonal complement of alpha    

elseif experiment == 3

    beta = [1;-1];                  % Normalized cointegrating vector                             **/
    beta0 = 2;                      % Long-run constant                                           **/
    beta1 = 0;                      % Long-run time trend                                         **/

    alpha = [-0.1;0.1];             % Error-correction parameters                                 **/

    alpha0 = null(alpha*alpha');
    delta0 = alpha0*2;              % Short-run intercepts using orthogonal complement of alpha   **/
    delta1 = alpha0*0;              % Short-run time trend using orthogonal complement of alpha   **/

else

    beta = [1;-1];                  % Normalized cointegrating vector                             **/
    beta0 = 2;                      % Long-run constant                                           **/
    beta1 = 0.1;                    % Long-run time trend                                         **/

    alpha = [-0.1;0.1];             % Error-correction parameters                                 **/

    alpha0 = null(alpha*alpha');
    delta0 = alpha0*0;              % Short-run intercepts using orthogonal complement of alpha   **/
    delta1 = alpha0*0;              % Short-run time trend using orthogonal complement of alpha   **/

end

% Simulate the model   

y = zeros(t,2);         
v = randn(t,2);

for j = 2:t
    
    u = beta0 + beta1*j + y(j-1,:)*beta;
    y(j,:) = y(j-1,:) + (delta0 + delta1*j + alpha*u)' + v(j,:);

end

% Graph simulated data 
figure(1)
clf
plot(1:t,y)  
legend('y1','y2','Location','SouthEast')


