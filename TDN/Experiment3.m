%% Experiment 3: 1D test of the Davidenko-Branin method

% Problem parameters
f = @(x) 0.5*x + sin(x);
fp= @(x) 0.5   + cos(x);

f = @(x) -f(x);
fp = @(x) -fp(x);

% Initialization
% x0 = 5;
x = linspace(-10,10,1001);
xN0 = x;
xB0 = x;
h = 1;
% plot(x,f(x),'k',xN0,f(xN0),'ro')
% axis([x(1),x(end),-10,10])
% pause(0.5)
for k = 1:100
    xN1 = xN0 - h*f(xN0)./(fp(xN0));
    xN0 = xN1;
    
    xB1 = xB0 - h*f(xB0)./abs(fp(xB0));
    xB0 = xB1;
    
%     plot(x,f(x),'k',x1,f(x1),'ro')
%     axis([x(1),x(end),-10,10])
%     pause(h)
end

indN = abs(xN1) > 1e-8;
indB = abs(xB1) > 1e-8;

yyaxis left
plot(x(indN),2*indN(indN),'b.',x(indB),indB(indB),'r.')
xlabel('x_0')
ylabel('Failure intervals')
axis([x(1),x(end),0,3])

yyaxis right
plot(x,f(x),'k')
ylabel('f(x_0)')
legend('Newton','Davidenko-Branin','f(x_0)')
title('Convergence for f(x) = -x/2 - sin(x)')

% This method is designed to get over 'humps' found while converging to
% zero using Newton's method. However, this requires that you have some
% knowledge of where the root is with regards to the 'hump': if the root is
% on the other side of the hump then this works fine, but if the root is on
% this side of the hump then you're going in the wrong direction. In many
% circumstances you will have no a priori knowledge of where the root is.