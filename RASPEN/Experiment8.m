%% Experiment 8 - seeking to break the method
% Having shown examples where preconditioning by Newton severly improves
% the convergence, we'd now like to find an example where the precondition
% fails. We do so by first finding a function G(y) for which Newton does
% poorly, or has a very small basin of attraction. We then find a PDE that
% would produce such a G(y) as part of the preconditioning.

%% Finding a function G(y)

% Grid
y = linspace(-1,1,1001);
a = 0.25;
b = 0.5;
G = @(y,c) -a*y - (b/c) * sin(c*y);
Gp= @(y,c) -a   - b*cos(c*y);

c = 10;
tol = 1e-6;
itermax = 100;
save = zeros(size(y));

for k = 1:length(y)
    y0 = y(k);
    error= 1;
    iter = 1;
    while error > tol && iter < itermax
%         y1 = y0 - (G(y0,c) - y0)./(Gp(y0,c) - 1);
%         y1 = y0 - G(y0,c)./Gp(y0,c);
        y1 = G(y0,c)+y0;
        error = abs(y1 - y0);
        y0 = y1;
        iter = iter+1;
        
%         plot(y,G(y,c),'k',y0,G(y0,c),'ro')
%         axis([-1,1,-1,1])
%         pause(0.1)
    end
    
    if error < tol && iter < itermax
        save(k) = iter;
    else
        save(k) = NaN;
    end

end

plot(y,save,'b*--')