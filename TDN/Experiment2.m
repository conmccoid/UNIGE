%% Experiment 2: pathing methods
% (renaming the algorithm)
% We look at a number of examples containing a variety of features that may
% complicate certain pathing methods.

%% Choosing between two roots

clear
clc

f = @(x) x.^2 - 1;
fp = @(x) 2*x;
x = -2:0.0001:2;

dt = 1/2^12;

p = -0.001;
a = f(p);
t = 0;
iter = 0;

while t<1
    
    k1 = -a./fp(p);
    k2 = -a./fp(p + dt*k1/2);
    k3 = -a./fp(p + dt*k2/2);
    k4 = -a./fp(p + dt*k3);
    
    p = p + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    t = t+dt;
    iter = iter+1;
    
    if mod(iter,100)==0
        plot(x,f(x),'b',p,f(p),'ro',p,a*(1-t),'k*')
        pause(0.01)
    end
    
end

%% Between two local extrema
% For this particular path and function the initial guess must be greater
% than 1 in absolute value.

clear
clc
f = @(x) x.^4 - 2*x.^2 - 1;
fp = @(x) 4*x.^3 - 4*x;
x = -2:0.0001:5;

dt = 1/2^13;

p = 5;
a = f(p);
t = 0;
iter = 0;

plot(x,f(x),'b',p,f(p),'ro',p,a*(1-t),'k*')
hold on

while t<1
    
    k1 = -a./fp(p);
    k2 = -a./fp(p + dt*k1/2);
    k3 = -a./fp(p + dt*k2/2);
    k4 = -a./fp(p + dt*k3);
    
    p = p + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    t = t+dt;
    iter = iter+1;
    
    if mod(iter,1000)==0
        plot(p,f(p),'ro',p,a*(1-t),'k*')
    end
    
end

hold off

%% Function with singularity: 1 - 1/x
% Can we circle around it into the complex plane?

clear
clc

x = -2:0.01:2;
y = x;
[xx,yy] = meshgrid(x,y);
rr = xx + 1i*yy;

f = @(x) 1 - 1./x;
fp = @(x) 1./x.^2;

% contour(x,y,real(f(xx+1i*yy)),'r')
% hold on
% contour(x,y,imag(f(xx+1i*yy)),'b')
% hold off

p = -1;
a = f(p);
t = 0;
iter = 0;
dt = 1/2^8;

g = @(t) a*(1-t) + 1i*sin(pi*t);
gp = @(t) -a + 1i*pi*cos(pi*t);

% g = @(t) a*(1-t);
% gp = @(t) -a;

figure(1)
plot(real(f(p)),imag(f(p)),'ro',real(g(t)),imag(g(t)),'k*',...
    real(p),imag(p),'b^')
hold on

while t<1
    
    k1 = gp(t)./fp(p);
    k2 = gp(t + dt/2)./fp(p + dt*k1/2);
    k3 = gp(t + dt/2)./fp(p + dt*k2/2);
    k4 = gp(t + dt)./fp(p + dt*k3);
    
    p = p + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    t = t+dt;
    iter = iter+1;
    
    if mod(iter,2^5)==0
        figure(1)
        plot(real(f(p)),imag(f(p)),'ro',real(g(t)),imag(g(t)),'k*',...
            real(p),imag(p),'b^')
    end
    
end

figure(1)
axis([-2,2.5,-0.5,1.5])
xlabel('Re(r)')
ylabel('Im(r)')
legend('f(p)','g(t)','p')
hold off

figure
plot(real(f(rr)),imag(f(rr)))