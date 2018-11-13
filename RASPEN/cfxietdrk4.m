%% Continuous fixed point iteration tests
% Casual tests of the continuous fixed point iteration method using ETDRK4
% bug: not converging to correct value

% Complex plane
x = linspace(17,19,201);
y = linspace(-1,1,201);
[xx,yy] = meshgrid(x,y);
zz = xx + 1i*yy;

% g = @(x) exp(-x);
% g = @(x) (1 + x)./(1 + exp(x));
% g = @(x) x + 1 - x.*exp(x);
g = @(x) x - (x.*exp(x) - 1)./((x+1).*exp(x));
G = 0;
h = 0.1;
T = round(200/h);

for i = 0:T
    
    % integrate and iterate
    k1 = g(zz);
    k2 = g(zz + h*k1/2)*exp(h/2);
    k3 = g(zz + h*k2/2)*exp(h/2);
    k4 = g(zz + h*k3)  *exp(h);
    zz = exp(-h)*(zz + h*(k1 + 2*k2 + 2*k3 + k4)/6);
        
    % plot
    figure(1)
%     subplot(2,1,1)
%     contourf(xx,yy,real(zz))
%     subplot(2,1,2)
%     contourf(xx,yy,imag(zz))
    contourf(xx,yy,log10(abs(zz - lambertw(0,1))))
    colorbar
    pause(eps)
end