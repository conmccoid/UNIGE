%% Fractal form
% Produces a picture of the fractal basin of attraction of Newton's method
% for the Lambert W function.

clear
clc

% Define the complex space we examine
limit = 10;
preal = linspace(-limit,limit,401);
pimag = preal;

N = length(preal); M = length(pimag);

% Initialize storage
iter = zeros(N,M);
r = iter;

% Set end conditions
itermax = 100;
tol = 1e-8;

for i = 1:M
    for j = 1:N
        
        % Pick out one of the points of concern
        p = preal(j) + 1i*pimag(i);
        
        error = 1;
        iteration = 0;
        
        while error > tol && iteration < itermax
            
            iteration = iteration + 1;
            
            % Newton's iteration
%             exp_p = exp(p);
%             p_new = p - (p * exp_p - 1)./( exp_p + p * exp_p );

            % Condensed Newton's iteration
            p_new = ( p^2 + exp(-p) )./( p+1 );

            error = norm(p_new - p);
            p = p_new;
            
        end
        
        % For plotting purposes, pick out successes and remove failures
        if iteration < itermax && isfinite(p)
            iter(i,j) = iteration;
            r(i,j) = p;
        else
            iter(i,j) = NaN;
            r(i,j) = NaN;
        end
    
    end
end

figure(1)
contourf(preal,pimag,iter)
hold on
plot(real(r),imag(r),'r.')  % Where are the roots?
hold off
xlabel('Re(p_0)')
ylabel('Im(p_0)')
colorbar
axis([-limit,limit,-limit,limit])
title('Newton''s method basin of attraction')