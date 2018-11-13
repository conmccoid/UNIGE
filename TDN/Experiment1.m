%% Experiment 1 for time dependent Newton (TDN)
% We try out TDN using some choices of the g(t) function and the example
% xe^x - 1 = 0, which converges for initial guesses in [0.5,1] (and
% elsewhere).

clear
clc

f = @(x) x.*exp(x) - 1;
fp = @(x) (x+1).*exp(x);
% f = @(x) log(x) + x;
% fp = @(x) 1 + 1./x;
x = -5:0.0001:1;

%% g(t) = a exp(-t)
% This should result in standard Newton's method

dt = 0.001;

p = 1;
a = f(p);
t = 0;

while t < 10
    p = p - dt*a./fp(p);
    t = t+dt;
    a = a*exp(-dt);
    
    plot(x,f(x),'b',p,f(p),'ro',p,a*sin(1/t),'k*')
    pause(dt)
end

%% g(t) = a sin(1/t)

dt = 0.001;

p = 1;
a = f(p);
t = 2/pi;

while t < 10
    p = p - dt*a*cos(1/t)./(t^2 * fp(p));
    t = t + dt;
    
    plot(x,f(x),'b',p,f(p),'ro',p,a*sin(1/t),'k*')
    pause(dt)
end

%% g(t) = a sin(t)

dt = 0.001;

p = 0.75;
a = f(p);
t = pi/2;

while t<3*pi/2
    p = p + dt*a*cos(t)./fp(p);
    t = t+dt;
    
    plot(x,f(x),'b',p,f(p),'ro',p,a*sin(t),'k*')
    pause(dt)
end

%% g(t) = a * (1-t)

dt = 0.001;

p = -0.5;
a = f(p);
t = 0;

while t<1
    p = p - dt*a./fp(p);
    t = t+dt;
    
    plot(x,f(x),'b',p,f(p),'ro',p,a*(1-t),'k*')
    pause(dt)
end

%% g(t) = a * (1 - t) with ode45

f = @(x) x.^2 - 1;
fp = @(x) 2*x;

NN = 0:8; NN = 2.^NN;
error = zeros(size(NN));

for i = 1:length(NN)
    N = NN(i);

    p = 0.5;
    a = f(p);
    t0 = 0;
    tf = 1;
    tspan = linspace(t0,tf,N+1);

    func = @(t,x) -a./fp(x);
%     [tout, yout] = ode45(func,tspan,p); p = yout(end);
%     for k = 0:N-1
%         p = p + func(k/N,p)/N;
%     end
    p = rk4(func,tspan,p);

%     error(i) = abs(p - lambertw(0,1));
    error(i) = abs(p - 1);
end

p = polyfit(log10(1./NN),log10(error),1);
P = polyval(p,log10(1./NN));
plot(log10(1./NN),log10(error),'bo--',log10(1./NN),P,'r-')
xlabel('log10(h)')
ylabel('log10(Error)')
legend('Data',[ num2str(10^p(2)) ' h^{' num2str(p(1)) '}'])