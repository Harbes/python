%=========================================================================
%
%   Bivariate normal distribution properties
%
%=========================================================================

clear all
clc

y1 = -4:0.2:4;
y2 = -4:0.2:4;

%********************************************************************
%***
%***     Generate graph
%***
%********************************************************************
    
% Switch off TeX interpreter and clear figure
set(0,'defaulttextinterpreter','none');
figure(1);
clf;


mu  = [0 0]';
Sig = [1 0.6; 0 0.6];

f = zeros( length(y1),length(y2) );
for i=1:length(y1)
    for j=1:length(y2)
    
        f(i,j)=(1/(2*pi*det(Sig)^1/2)) ...
            *exp((-1/2)*([y1(i) y2(j)]-mu')*inv(Sig)*([y1(i);y2(j)]-mu));
    end
end

subplot(2,2,1)
mesh(y1,y2,f)
title('$\rho = 0.6$')
xlabel('$y_1$')
ylabel('$y_2$')
zlim([0 0.4]);
zlabel('$f(y_1,y_2)$')


subplot(2,2,3)
contour(y1,y2,f)
xlabel('$y_1$')
ylabel('$y_2$')


mu  = [0 0]';
Sig = [1 0; 0 1];

f = zeros( length(y1),length(y2) );
for i=1:length(y1)
    for j=1:length(y2)
    
        f(i,j)=(1/(2*pi*det(Sig)^1/2)) ...
            *exp((-1/2)*([y1(i) y2(j)]-mu')*inv(Sig)*([y1(i);y2(j)]-mu));
    end
end

subplot(2,2,2)
mesh(y1,y2,f)
title('$\rho=0.0$')
xlabel('$y_1$')
ylabel('$y_2$')
zlim([0 0.4])
zlabel('$f(y_1,y_2)$')


subplot(2,2,4)
contour(y1,y2,f)
xlabel('$y_1$')
ylabel('$y_2$')


laprint(1,'bivariatenormal','options','factory');

