%% Preconditioning Newton's method using a fixed point iteration
% For a starting example, consider the LambertW function x*exp(x)-1=0.
% This has an infinite number of roots throughout the complex plane.
% We intend to solve it by constructing first a fixed point iteration
% g(x) = x and then using Newton's method on F(x) = g(x) - x.

%% Non-preconditioned
% f(x) = x * exp(x) - 1 = 0

clear
clc

pp = logspace(-8,20,100);
iter = zeros(size(pp));
itermax = 1e8;
tol = 1e-8;

r = iter;
for n = 1:100

    p = pp(n);
    error = 1;
    iteration = 0;

    while error > tol && iteration < itermax
        iteration = iteration + 1;

        exp_p = exp(p);
        p_new = p - (p * exp_p - 1)./( exp_p + p * exp_p );

        error = norm(p_new - p);
        p = p_new;
    end

    if iteration < itermax && isfinite(p)
        disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
        iter(n) = iteration;
        r(n) = p;
    else
        disp(['Failure to converge within ' num2str(itermax) ' iterations'])
        iter(n) = NaN;
        r(n) = NaN;
    end
    
end

figure(1)
semilogx(pp,iter,'b.--')
xlabel('Initial guess')
ylabel('Number of iterations to convergence')

figure(2)
loglog(pp,abs(r - lambertw(1)),'r.--')
xlabel('Initial guess')
ylabel('Root found')

%% g(x) = x^2 * exp(x)
% g'(x) = 2*x *exp(x) + x^2 * exp(x)
% F(x) = x^2 * exp(x) - x; F'(x) = 2*x * exp(x) + x^2 * exp(x) - 1;
% Newton's iteration: p_n+1 = p_n - F(p_n)/F'(p_n)

clear
clc

pp = logspace(-2,3,100);
iter = zeros(size(pp));
itermax = 1e8;
tol = 1e-8;


r = iter;
for n = 1:100

    p = pp(n);
    error = 1;
    iteration = 0;

    while error > tol && iteration < itermax
        iteration = iteration + 1;
        
        exp_p = exp(p);
        p_new = p - (p.^2 .* exp_p - p)./( 2*p .* exp_p + ...
            p.^2 .* exp_p - 1);

        error = norm(p_new - p);
        p = p_new;
    end

    if iteration < itermax && isfinite(p)
        disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
        iter(n) = iteration;
        r(n) = p;
    else
        disp(['Failure to converge within ' num2str(itermax) ' iterations'])
        iter(n) = NaN;
        r(n) = NaN;
    end
    
end

figure(1)
semilogx(pp,iter,'b.--')
xlabel('Initial guess')
ylabel('Number of iterations to convergence')

figure(2)
loglog(pp,abs(r - lambertw(1)),'r.--')
xlabel('Initial guess')
ylabel('Root found')

%% g(x) = exp(-x)
% g'(x) = -exp(-x)
% F(x) = x - exp(-x); F'(x) = 1 + exp(-x);
% Newton's iteration: p_n+1 = p_n - F(p_n)/F'(p_n)

clear
clc

pp = logspace(-8,8,100);
iter = zeros(size(pp));
itermax = 1e8;
tol = 1e-8;


r = iter;
for n = 1:100

    p = pp(n);
    error = 1;
    iteration = 0;

    while error > tol && iteration < itermax
        iteration = iteration + 1;
        
        exp_p = exp(-p);
        p_new = p - (p - exp_p)./(1 + exp_p);

        error = norm(p_new - p);
        p = p_new;
    end

    if iteration < itermax && isfinite(p)
        disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
        iter(n) = iteration;
        r(n) = p;
    else
        disp(['Failure to converge within ' num2str(itermax) ' iterations'])
        iter(n) = NaN;
        r(n) = NaN;
    end
    
end

figure(1)
semilogx(pp,iter,'b.--')
xlabel('Initial guess')
ylabel('Number of iterations to convergence')

figure(2)
loglog(pp,abs(r - lambertw(1)),'r.--')
xlabel('Initial guess')
ylabel('Root found')

%% g(x) = -log(x)
% g'(x) = -1/x;
% F(x) = -log(x) - x; F'(x) = -1/x - 1;

clear
clc

pp = logspace(-8,20,100);
% pp = 1;
iter = zeros(size(pp));
itermax = 1e2;
tol = 1e-8;

r = iter;
for n = 1:length(pp)

    p = pp(n);
    error = 1;
    iteration = 0;
%     err_list = [];

    while error > tol && iteration < itermax
        iteration = iteration + 1;
        
        logp = log(p);
        p_new = p - (-logp - p)./(-1/p - 1);
        
        error = norm(p_new - p);
%         err_list = [ err_list error ];
        p = p_new;
        
    end

    if iteration < itermax && isfinite(p)
        disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
        iter(n) = iteration;
        r(n) = p;
    else
        disp(['Failure to converge within ' num2str(itermax) ' iterations'])
        iter(n) = NaN;
        r(n) = NaN;
    end
    
end

figure(1)
semilogx(pp,iter,'b.--')
xlabel('Initial guess')
ylabel('Number of iterations to convergence')

figure(2)
loglog(pp,abs(r - lambertw(1)),'r.--')
xlabel('Initial guess')
ylabel('Root found')

for n = 1:100
    
    iteration = 0;
    
    fp = pp(n);
    error = 1;

    while error > tol && iteration < itermax
        iteration = iteration + 1;

        fp_new = -log(fp);

        error = norm(fp_new - fp);
        fp = fp_new;
        
    end

    if iteration < itermax && isfinite(fp)
        disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
        iter(n) = iteration;
        r(n) = fp;
    else
        disp(['Failure to converge within ' num2str(itermax) ' iterations'])
        iter(n) = NaN;
        r(n) = NaN;
    end
    
end

figure(3)
semilogx(pp,iter,'b.--')
xlabel('Initial guess')
ylabel('Number of iterations to convergence')

figure(4)
loglog(pp,abs(r - lambertw(1)),'r.--')
xlabel('Initial guess')
ylabel('Root found')

% Using this fixed point significantly preconditions the problem.
% It should be noted that for no values of initial guess does the
% corresponding fixed point iteration converge. It always ends up in a
% two-cycle, neither point corresponding to a root of the Lambert W
% function.

%% 2D comparison of non-preconditioned and -log(x) preconditioned problem
% We now look over the complex plane for the entire basin of attraction. We
% take a fairly fine grid of the complex plane and check each point as an
% initial guess. We will begin by considering only the region sufficiently
% close to the first few roots: from -2 to 2 along the real axis and from
% -30 to 30 along the imaginary axis

clear
clc

limit = 3;
preal = linspace(-limit,limit,401);
pimag = preal;

branch = 0;

N = length(preal); M = length(pimag);

iter_Newton = zeros(N,M);
iter_Precon = iter_Newton;
r_Newton = iter_Newton;
r_Precon = iter_Newton;
iter_Precon2 = iter_Newton;
r_Precon2 = iter_Newton;

itermax = 100;
tol = 1e-8;

for i = 1:M
    for j = 1:N
        
        p = preal(j) + 1i*pimag(i);
        error = 1;
        iteration = 0;
        
        while error > tol && iteration < itermax
            
            iteration = iteration + 1;
            
%             exp_p = exp(p);
%             p_new = p - (p * exp_p - 1)./( exp_p + p * exp_p );
            p_new = ( p^2 + exp(-p) )./( p+1 );

            error = norm(p_new - p);
            p = p_new;
            
        end
        
        if iteration < itermax && isfinite(p)
%             disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
            iter_Newton(i,j) = iteration;
            r_Newton(i,j) = p;
        else
%             disp(['Failure to converge within ' num2str(itermax) ' iterations'])
            iter_Newton(i,j) = NaN;
            r_Newton(i,j) = NaN;
        end
        
        p = preal(j) + 1i*pimag(i);
        error = 1;
        iteration = 0;
        
        while error > tol && iteration < itermax
            
            iteration = iteration + 1;
            
            p_new = p - (-log(p) - p + 1i*branch*2*pi)./(-1/p - 1);
        
            error = norm(p_new - p);
            p = p_new;
            
        end
        
        if iteration < itermax && isfinite(p)
%             disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
            iter_Precon(i,j) = iteration;
            r_Precon(i,j) = p;
        else
%             disp(['Failure to converge within ' num2str(itermax) ' iterations'])
            iter_Precon(i,j) = NaN;
            r_Precon(i,j) = NaN;
        end
        
        % we attach additionally the preconditioning with x = exp(-x)
        p = preal(j) + 1i*pimag(i);
        error = 1;
        iteration = 0;
        
        while error > tol && iteration < itermax
            
            iteration = iteration + 1;
            
            p_new = p - (p - exp(-p))./(1 + exp(-p));
        
            error = norm(p_new - p);
            p = p_new;
            
        end
        
        if iteration < itermax && isfinite(p)
%             disp(['Convergence in ' num2str(iteration) ' iteration(s)'])
            iter_Precon2(i,j) = iteration;
            r_Precon2(i,j) = p;
        else
%             disp(['Failure to converge within ' num2str(itermax) ' iterations'])
            iter_Precon2(i,j) = NaN;
            r_Precon2(i,j) = NaN;
        end
    
    end
end

figure(1)
contourf(preal,pimag,iter_Newton)
hold on
plot(real(r_Newton),imag(r_Newton),'r.')
hold off
xlabel('Re(p_0)')
ylabel('Im(p_0)')
colorbar
axis([-limit,limit,-limit,limit])
title('Newton''s method basin of attraction')

figure(2)
contourf(preal,pimag,iter_Precon)
hold on
plot(real(r_Precon),imag(r_Precon),'r.')
hold off
xlabel('Re(p_0)')
ylabel('Im(p_0)')
colorbar
axis([-limit,limit,-limit,limit])
title('log(x) preconditioner basin of attraction')

% The basin of attraction of Newton's method for this problem has fractal 
% boundaries. The preconditioned Newton's method completely removes these,
% extending the basin to essentially the entire complex plane. However, it
% is necessary to pick which branch of the function we are interested in
% before using the preconditioned method. Also, the points -1 and 0 are not
% ideal for initial guesses.
%%
figure(3)
contourf(preal,pimag,real(r_Precon-lambertw(branch,1)))
title('Accuracy of log(x) preconditioner: real')

figure(4)
contourf(preal,pimag,imag(r_Precon-lambertw(branch,1)))
title('Accuracy of log(x) preconditioner: imag')

% same results for exp(-x) preconditioner
figure(5)
contourf(preal,pimag,iter_Precon2)
hold on
plot(real(r_Precon2),imag(r_Precon2),'r.')
hold off
xlabel('Re(p_0)')
ylabel('Im(p_0)')
colorbar
axis([-limit,limit,-limit,limit])
title('exp(-x) preconditioner basin of attraction')

figure(6)
contourf(preal,pimag,real(r_Precon2-lambertw(branch,1)))
title('Accuracy of exp(-x) preconditioner: real')

figure(7)
contourf(preal,pimag,imag(r_Precon2-lambertw(branch,1)))
title('Accuracy of exp(-x) preconditioner: imag')

%% Testing the progression of the errors for the preconditioned method
% We now examine the errors at each step, seeing how they progress through
% 12 iterations.

clear
clc

limit = 1;
preal = linspace(-limit,limit,101);
pimag = preal;
[ppreal,ppimag] = meshgrid(preal,pimag');
p = ppreal + 1i*ppimag + lambertw(0,1);

branch = 0;

for iteration = 1:1
    
    p_new = p - (log(p) + p -1i*branch*2*pi)./(1 + 1./p);
%     p_new = ( p.^2 + exp(-p) )./(p+1);
    
    error = p_new - lambertw(branch,1);
    p = p_new;
    
    figure(1)
    contourf(preal,pimag,real(error))
    colorbar
    
    figure(2)
    contourf(preal,pimag,log(abs(imag(error))))
    colorbar
    
    figure(3)
    contourf(preal,pimag,real(p))
    colorbar
    
    figure(4)
    contourf(preal,pimag,imag(p))
    colorbar
    
    pause
end

figure(5)
contour(ppreal,ppimag,real(p-lambertw(0,1)),'b')
hold on
contour(ppreal,ppimag,imag(p-lambertw(0,1)),'r')
hold off

%% Testing basin of attraction
% We propose that the basin of attraction is the region for which
% f''(p)/f(p) < 2/|p* - p|. We test this for the current experiment.

limit = 1e0;
pp = linspace(-limit,limit,401);
rr = lambertw(0,1);

iter = zeros(size(pp));
r = iter;
tol = 1e-8;
itermax = 100;

for i = 1:length(pp)
    
    p = pp(i);
    error = 1;
    iteration = 0;

    while error > tol && iteration < itermax
        iteration = iteration + 1;
        
        % Newton's method
%         p_new = p - (p - exp(-p))./( 1 + p );
        
        % Precon. with log(x)
        p_new = p - (log(p) + p)./(1/p + 1);
        
        % Precon. with exp(-x)
%         p_new = p - (p - exp(-p))./(1 + exp(-p));

        error = norm(p_new - p);
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

yyaxis left
plot(pp,iter,'r')

yyaxis right
% plot(pp,abs( abs((pp+2).*(pp-rr)./(pp+1)) )<2,'b')
plot(pp,abs( (pp-rr)./(pp.*(pp+1)) )<2,'b')
% plot(pp,abs( (pp-rr)./(exp(pp)+1) )<2,'b')