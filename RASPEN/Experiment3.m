%% Experiment 3
% f(x) = x^4 + 4x^3 - 2x^2 + 1 = 0

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
        
        p_new = p - (p.^4 + 4*p.^3 - 2*p.^2 + 1)./(4*p.^3 + 12*p.^2 - 4*p);
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

figure(2)
plot(pp,r,'b.')

%% Preconditioned with g(x) = -(-4x^3 + 2x^2 -1)^(1/4)

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
        
        g1 = -4*p.^3 + 2*p.^2 - 1;
        g2 = -12*p.^2 + 4*p;
        p_new = p - (p + g1.^(1/4))./(1 + g2./(4*g1.^(3/4)));
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
hold on
plot(pp,iter,'r.--')
hold off

figure(2)
hold on
plot(pp,r,'r.')
hold off

% This does not converge well at all.

%% Preconditioned with g(x) = -sqrt(0.5x^4 + 2x^3 + 0.5)

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
        
        g1 = 0.5*p^4 + 2*p^3 + 0.5;
        p_new = p - (p + sqrt(g1))./(1 + 0.5*(2*p^3 + 6*p^2)./sqrt(g1));
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
hold on
plot(pp,iter,'k.--')
hold off

figure(2)
hold on
plot(pp,r,'k.')
hold off

% This is a rather nice preconditioner, except near the cutoff between
% roots.

%% Preconditioned with g(x) = 0.5x^3 + 2x^2 + 1/2x

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
        
        g1 = 0.5*p^4 + 2*p^3 + 0.5;
        p_new = p - (p -0.5*p.^3 - 2*p.^2 - 0.5./p)./(1 ...
            - 1.5*p.^2 - 4*p + 0.5./p.^2);
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
hold on
plot(pp,iter,'m.--')
hold off

figure(2)
hold on
plot(pp,r,'m.')
hold off

% Not great, possibly better but not likely.

%% Preconditioned with g(x) = (-0.25x^4 + 0.5x^2 - 0.25)^(1/3)

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
        
        g1 = -0.25*p.^4 + 0.5*p.^2 - 0.25;
        p_new = p - (p - g1.^(1/3))./(1 - (-p.^3 + p)./(3*(g1^(2/3))));
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
hold on
plot(pp,iter,'g.--')
xlabel('Initial guess')
ylabel('Number of iterations')
legend('Newton','()^(1/4)','()^(1/2)','(..+1/x)','()^(1/3)')
hold off

figure(2)
hold on
plot(pp,r,'g.')
xlabel('Initial guess')
ylabel('Root')
legend('Newton','()^(1/4)','()^(1/2)','(..+1/x)','()^(1/3)')
hold off