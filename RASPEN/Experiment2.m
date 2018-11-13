%% Experiment 2: Kepler's equation
% As in experiment 1, we look to precondition Newton's method by finding a
% function that shares a root with the function of concern. The function of
% concern here is Kepler's equation: x - 0.8 sin(x) - 2 pi / 10.

%% Newton's method
clear
clc

limit = 1e1;
offset = 0;
pp = linspace(-limit,limit,401) + offset;

iter = zeros(size(pp));
r = iter;
itermax = 1e2;
tol = 1e-8;

for i = 1:length(pp)
    
    p = pp(i);
    iteration = 0;
    error = 1;
    
    while error > tol && iteration < itermax
        
        iteration = iteration + 1;
        
        p_new = p - (p - 0.8*sin(p) - 2*pi/10)./(1 - 0.8*cos(p));
        error = abs(p_new - p);
        p = p_new;
        
    end
    
    if iteration < itermax && isfinite(p)
        iter(i) = iteration;
        r(i) = p;
    else
        iter(i) = NaN;
        r(i) = NaN;
    end
    
end

figure(1)
plot(pp,iter,'b.--')
xlabel('Initial guess')
ylabel('Number of iterations')

figure(2)
plot(pp,r,'b.')
xlabel('Initial guess')
ylabel('Root')

% This converges perfectly fine on its own, there seems no need for
% preconditioning. Moreover, finding a different function that shares a
% root that works better may prove challenging.

%% Complicating the function: m-th power
% We look for more complicated functions that share a root with this
% function, to see where Newton's method might fail. We start by taking the
% m-th power:
% (x - 0.8 sin(x) - c)^m = 0

m = 2;

iter = zeros(size(pp));
r1 = iter;

for i = 1:length(pp)
    
    p = pp(i);
    iteration = 0;
    error = 1;
    
    while error > tol && iteration < itermax
        
        iteration = iteration + 1;
        g1 = p - 0.8*sin(p) - 2*pi/10;
        p_new = p - g1^m ./ (m * (g1^(m-1)) * (1 - 0.8*cos(p)));
        
        error = abs(p_new-p);
        p = p_new;
        
    end
    
    if iteration < itermax && isfinite(p)
        iter(i) = iteration;
        r1(i) = p;
    else
        iter(i) = NaN;
        r1(i) = NaN;
    end
    
end

figure(1)
hold on
plot(pp,iter,'r.--')
hold off

figure(2)
hold on
plot(pp,r1,'r.')
hold off

%% Complicating the function: log(f(x) + 1)

iter = zeros(size(pp));
r1 = iter;

for i = 1:length(pp)
    
    p = pp(i);
    iteration = 0;
    error = 1;
    
    while error > tol && iteration < itermax
        
        iteration = iteration + 1;
        g1 = p - 0.8*sin(p) - 2*pi/10;
        g2 = 1 - 0.8*cos(p);
        p_new = p - (g1+1).*log(g1 + 1)./g2;
        
        error = abs(p_new-p);
        p = p_new;
        
    end
    
    if iteration < itermax && isfinite(p)
        iter(i) = iteration;
        r1(i) = p;
    else
        iter(i) = NaN;
        r1(i) = NaN;
    end
    
end

figure(1)
hold on
plot(pp,iter,'k.--')
hold off

figure(2)
hold on
plot(pp,r1,'k.')
hold off

%% Complicating the function: exp(-1/f(x))

iter = zeros(size(pp));
r1 = iter;

for i = 1:length(pp)
    
    p = pp(i);
    iteration = 0;
    error = 1;
    
    while error > tol && iteration < itermax
        
        iteration = iteration + 1;
        g1 = p - 0.8*sin(p) - 2*pi/10;
        g2 = 1 - 0.8*cos(p);
        p_new = p - sign(g1)*abs(g1)^2/g2;
        
        error = abs(p_new-p);
        p = p_new;
        
%         plot(p,'ro')
%         pause(0.1)
        
    end
    
    if iteration < itermax && isfinite(p)
        iter(i) = iteration;
        r1(i) = p;
    else
        iter(i) = NaN;
        r1(i) = NaN;
    end
    
end

figure(1)
hold on
plot(pp,iter,'m.--')
legend('Newton','f^2','log(f+1)','exp(-1/f)')
hold off

figure(2)
hold on
plot(pp,r1,'m.')
legend('Newton','f^2','log(f+1)','exp(-1/f)')
hold off