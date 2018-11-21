%% Experiment 7 - Repeat experiment 6 for viscous Burgers

% Problem parameters
e = 0.1;

% Grid
nx = 1001; dx = 1/(nx-1);
a  = -1/3;
b  = 1/3;
x  = linspace(-1,1,nx)';
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);
l1 = length(x1); l2 = length(x2);

% Diff. mat.
I1 = eye(l1);
D1 = ones(l1-2,1);
D1 = (nx-1)^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[0,1,2],l1-2,l1);
D1 = [I1(1,:) ; e*D1 ; I1(end,:) ];
I2 = eye(l2);
D2 = ones(l2-2,1);
D2 = (nx-1)^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[0,1,2],l2-2,l2);
D2 = [I2(1,:) ; e*D2 ; I2(end,:) ];

d1 = ones(l1-2,1);
d1 = (nx-1)/2 * [ -d1, zeros(l1-2,1), d1 ];
d1 = spdiags(d1,[0,1,2],l1-2,l1);
d1 = [ zeros(1,l1) ; d1 ; zeros(1,l1) ];
d2 = ones(l2-2,1);
d2 = (nx-1)/2 * [ -d2, zeros(l2-2,1), d2 ];
d2 = spdiags(d2,[0,1,2],l2-2,l2);
d2 = [ zeros(1,l2) ; d2 ; zeros(1,l2) ];

% Jacobian
J1 = ones(l1-2,1);
J1 = [J1, -2*J1, J1];
J1 = e * J1 * (nx-1)^2;
J1 = spdiags(J1,[-1,0,1],l1-2,l1-2);
J2 = ones(l2-2,1);
J2 = [J2, -2*J2, J2];
J2 = e * J2 * (nx-1)^2;
J2 = spdiags(J2,[-1,0,1],l2-2,l2-2);

% Initialization
k = 0;
tol = 1e-8;
itermax = 100;
itersave = zeros(101,1);
testu2 = linspace(-1,1,101);
for u2b0 = testu2
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = -x1; u1(1) = 1;
    u2 = -x2; u2(end) = -1;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        u1(end) = u2b;
        for j = 1:1
            f1 = -D1 * u1 + d1 * (u1.^2 / 2);
            JJ = spdiags([u1(2:end-1) , -u1(3:end)+u1(1:end-2) , -u1(2:end-1)],[-1,0,1],l1-2,l1);
            JJ = JJ(:,2:end-1);
            JJ = J1 + JJ*(nx-1)/2;
%             u1(2:end-1) = u1(2:end-1) + (JJ \ f1(2:end-1));
            JJ = JJ \ eye(l1-2); S1 = JJ(x1(2:end-1)==a,end); U1 = u1(end-1);
            u1(2:end-1) = u1(2:end-1) + JJ * f1(2:end-1);
        end
        u1a = u1(x1==a);

        % Step 2: solve u in second domain
        u2(1) = u1a;
        for j = 1:1
            f2 = -D2 * u2 + d2 * (u2.^2/2);
            JJ = spdiags([u2(2:end-1) , -u2(3:end)+u2(1:end-2) , -u2(2:end-1)],[-1,0,1],l2-2,l2);
            JJ = JJ(:,2:end-1);
            JJ = J2 + JJ*(nx-1)/2;
%             u2(2:end-1) = u2(2:end-1) + (JJ \ f2(2:end-1));
            JJ = JJ \ eye(l2-2); S2 = JJ(x2(2:end-1)==b,1); U2 = u2(2);
            u2(2:end-1) = u2(2:end-1) + JJ * f2(2:end-1);
        end
        u2b = u2(x2==b);

%         % Preconditioning with Newton
%         % Step 3: solve g in first domain
%         dF1 = diag(d1 * u1);
%         g1 = (D1 - dF1 - diag(u1)*d1) \ [ zeros(l1-1,1) ; 1 ];
%         
% 
%         % Step 4: solve g in second domain
%         dF2 = diag(d2 * u2);
%         g2 = (D2 - dF2 - diag(u2)*d2) \ [ g1(x1==a) ; zeros(l2-1,1) ];
%         g2b = g2(x2==b);
% 
%         % Step 5: update u2b
%         u2b = ( g2b*u2b - u2(x2==b) )./ ( g2b - 1 ); u1(end) = u2b;

        % New preconditioning, for Newton instead of fixed point as
        % nonlinear solver
%         g = (e*(nx-1)^2 + (2*U1 - u1(end-1))*(nx-1)/2) * ( (u2(2) - U2)*(nx-1)/2 - e*(nx-1)^2 ) * S1 * S2;
%         u2b = u2bold - (u2b - u2bold)/(g - 1);

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
    end
    
    if iter < itermax && error < tol
        itersave(k) = iter;
    else
        itersave(k) = NaN;
    end
    
    plot(x1,u1,x2,u2)
    axis([-1,1,-1,1])
    pause(eps)
        
end

plot(testu2,itersave)
xlabel('\gamma')
ylabel('Number of iterations to convergence')
title('Newton precond. on transmission condition')

%% Showing order of convergence
% Again we improve rate of convergence from linear to quadratic, but only
% for e >= 1e-1.

% Problem parameters
e = 0.1;

% Grid
nx = 1001; dx = 1/(nx-1);
L = 0.2;
a  = -L;
b  = L;
x  = linspace(-1,1,nx)';
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);
l1 = length(x1); l2 = length(x2);

% Diff. mat.
I1 = eye(length(x1));
D1 = ones(length(x1)-2,1);
D1 = (nx-1)^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[0,1,2],l1-2,l1);
D1 = [I1(1,:) ; e*D1 ; I1(end,:) ];
I2 = eye(length(x2));
D2 = ones(length(x2)-2,1);
D2 = (nx-1)^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[0,1,2],l2-2,l2);
D2 = [I2(1,:) ; e*D2 ; I2(end,:) ];

d1 = ones(l1-2,1);
d1 = (nx-1)/2 * [ -d1, zeros(l1-2,1), d1 ];
d1 = spdiags(d1,[0,1,2],l1-2,l1);
d1 = [ zeros(1,l1) ; d1 ; zeros(1,l1) ];
d2 = ones(l2-2,1);
d2 = (nx-1)/2 * [ -d2, zeros(l2-2,1), d2 ];
d2 = spdiags(d2,[0,1,2],l2-2,l2);
d2 = [ zeros(1,l2) ; d2 ; zeros(1,l2) ];

% Jacobian
J1 = ones(l1-2,1);
J1 = [J1, -2*J1, J1];
J1 = e * J1 * (nx-1)^2;
J1 = spdiags(J1,[-1,0,1],l1-2,l1-2);
J2 = ones(l2-2,1);
J2 = [J2, -2*J2, J2];
J2 = e * J2 * (nx-1)^2;
J2 = spdiags(J2,[-1,0,1],l2-2,l2-2);

% Initialization
u0 = -x;
itermax = 200;
tol = 1e-8;
error0_u = zeros(itermax+1,1);
error0_u(1) = 1;
error1_u = error0_u;
error2_u = error0_u;
error0_g = error0_u;
error1_g = error0_u;
error2_g = error0_u;

% Standard additive Schwarz
u2b = -b;
u2bold = -b;
u1 = u0(x<=b); u1(1) = 1;
u2 = u0(x>=a); u2(end) = -1;
uold = u0;
for iter = 1:itermax

    u1(end) = u2b;
        % Step 1: solve u in first domain
        f1 = -D1 * u1 + u1 .* ( d1 * u1 );
        JJ = spdiags([u1(2:end-1) , -u1(3:end)+u1(1:end-2) , -u1(2:end-1)],[-1,0,1],l1-2,l1);
        JJ = JJ(:,2:end-1);
        JJ = J1 + JJ*(nx-1)/2;
        JJ = JJ \ eye(l1-2); S1 = JJ(x1(2:end-1)==a,end);
        u1(2:end-1) = u1(2:end-1) + JJ * f1(2:end-1);
        
        u1a = u1(x1==a);
        u2(1) = u1a;
        
        % Step 2: solve u in second domain
        f2 = -D2 * u2 + u2 .* ( d2 * u2 );
        JJ = spdiags([u2(2:end-1) , -u2(3:end)+u2(1:end-2) , -u2(2:end-1)],[-1,0,1],l2-2,l2);
        JJ = JJ(:,2:end-1);
        JJ = J2 + JJ*(nx-1)/2;
        JJ = JJ \ eye(l2-2); S2 = JJ(x2(2:end-1)==b,1);
        u2(2:end-1) = u2(2:end-1) + JJ * f2(2:end-1);
    u2b = u2(x2==b);
    
    u = [ u1(x1<a) ; 0.5*u1(x1>=a) + 0.5*u2(x2<=b) ; u2(x2>b) ];
    
    error0_u(iter) = norm(u - uold);
    error0_g(iter) = abs(u2b - u2bold);
    uold = u;
    u2bold = u2b;

end

% Preconditioning with Newton (u then bc.)
u2b = -b;
u2bold = -b;
iter = 0;
u1 = u0(x<=b); u1(1) = 1;
u2 = u0(x>=a); u2(end) = -1;
u1old = u1;
u2old = u2;
nonlinsolves = 30;
while iter < itermax
    
    u1(end) = u2b;
    errorA = 1;
    iterA = iter;
%     while errorA > tol && iterA < itermax
    while iterA < iter + nonlinsolves && iterA < itermax
        iterA = iterA+1;
        
        % Step 1: solve u in first domain
        f1 = -D1 * u1 + u1 .* ( d1 * u1 );
        JJ = spdiags([u1(2:end-1) , -u1(3:end)+u1(1:end-2) , -u1(2:end-1)],[-1,0,1],l1-2,l1);
        JJ = JJ(:,2:end-1);
        JJ = J1 + JJ*(nx-1)/2;
        JJ = JJ \ eye(l1-2); S1 = JJ(x1(2:end-1)==a,end);
        u1(2:end-1) = u1(2:end-1) + JJ * f1(2:end-1);
        
        errorA = norm(u1 - u1old);
        u1old = u1;
        
        error1_u(iterA) = errorA;
    end
        
        u1a = u1(x1==a);
        u2(1) = u1a;
        
    errorB = 1;
    iterB = iter;
%     while errorB > tol && iterB < itermax
    while iterB < iter + nonlinsolves && iterB < itermax
        iterB = iterB+1;
        
        % Step 2: solve u in second domain
        f2 = -D2 * u2 + u2 .* ( d2 * u2 );
        JJ = spdiags([u2(2:end-1) , -u2(3:end)+u2(1:end-2) , -u2(2:end-1)],[-1,0,1],l2-2,l2);
        JJ = JJ(:,2:end-1);
        JJ = J2 + JJ*(nx-1)/2;
        JJ = JJ \ eye(l2-2); S2 = JJ(x2(2:end-1)==b,1);
        u2(2:end-1) = u2(2:end-1) + JJ * f2(2:end-1);
        
        errorB = norm(u2 - u2old);
        u2old = u2;
        
        error1_u(iterB) = error1_u(iterB) + errorB;
    end
    iter = min(iterA,iterB);
    u2b = u2(x2==b);
            
    % Preconditioning with Newton
    % Step 3: solve g in first domain
    dF1 = diag(d1 * u1);
    g1 = (D1 - dF1 - diag(u1)*d1) \ [ zeros(l1-1,1) ; 1 ];

    % Step 4: solve g in second domain
    dF2 = diag(d2 * u2);
    g2 = (D2 - dF2 - diag(u2)*d2) \ [ g1(x1==a) ; zeros(l2-1,1) ];

    % Step 5: update u2b
    g2b = g2(x2==b);
    u2b = ( g2b*u2bold - u2b )./ ( g2b - 1 );

    error1_g(iter+1) = abs(u2b - u2bold);
    u2bold= u2b;

end

% Preconditioning with Newton (bc. then u)
error2_g = zeros(4*itermax,1);
u2b = -b;
u2bold = -b;
u1 = u0(x<=b); u1(1) = 1;
u2 = u0(x>=a); u2(end) = -1;
u1old = u1;
u2old = u2;
uold = u0;
for iter = 1:4*itermax
        
    u1 = u1old;
    u1(end) = u2b;

    f1 = -D1 * u1 + u1 .* ( d1 * u1 );
    JJ = spdiags([u1(2:end-1) , -u1(3:end)+u1(1:end-2) , -u1(2:end-1)],[-1,0,1],l1-2,l1);
    JJ = JJ(:,2:end-1);
    JJ = J1 + JJ*(nx-1)/2;
    JJ = JJ \ eye(l1-2); S1 = JJ(x1(2:end-1)==a,end);
    u1(2:end-1) = u1(2:end-1) + JJ * f1(2:end-1);

    u1a = u1(x1==a);
    u2  = u2old;
    u2(1) = u1a;

    f2 = -D2 * u2 + u2 .* ( d2 * u2 );
    JJ = spdiags([u2(2:end-1) , -u2(3:end)+u2(1:end-2) , -u2(2:end-1)],[-1,0,1],l2-2,l2);
    JJ = JJ(:,2:end-1);
    JJ = J2 + JJ*(nx-1)/2;
    JJ = JJ \ eye(l2-2); S2 = JJ(x2(2:end-1)==b,1);
    u2(2:end-1) = u2(2:end-1) + JJ * f2(2:end-1);

    u2b = u2(x2==b);

    % New preconditioning, for Newton instead of fixed point as
    % nonlinear solver
    g = -( -e*(nx-1)^2 + u1(end-1)*(nx-1)/2 ) * ( e*(nx-1)^2 + u2(2)*(nx-1)/2 ) * S1 * S2;
    u2b = ( g*u2bold - u2b )./ (g - 1);

    error2_g(iter) = abs(u2b - u2bold);
    u2bold = u2b;

%     if error2_g(iter) < tol
    if mod(iter,4)==0

        u = [ u1(x1<a) ; 0.5*u1(x1>=a) + 0.5*u2(x2<=b) ; u2(x2>b) ];

        error2_u(iter/4) = norm(u - uold);
        u1old = u1;
        u2old = u2;
        uold = u;

    end

end

subplot(1,2,1)
semilogy(1:(itermax+1),error0_u,'b*--',1:(itermax+1),error1_u,'ro--',1:(itermax+1),error2_u,'k.--')
axis([0, 200, 1e-18, 1])
xlabel('Iteration')
ylabel('2-norm of change in u')
legend('AS','ASPN','PNAS')

subplot(1,2,2)
semilogy(1:(itermax+1),error0_g,'b*--',1:(itermax+1),error1_g,'ro--',1:itermax,error2_g(1:4:end),'k.--')
axis([0, 200, 1e-18, 1])
xlabel('Iteration')
ylabel('Change in u_2(\beta)')
legend('AS','ASPN','PNAS')