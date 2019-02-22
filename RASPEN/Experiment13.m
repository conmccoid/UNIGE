%% Experiment 13 - searching for a period doubling example
% As part of a series of counterexamples, I am searching for a period
% doubling bifurcation in the Newton-Raphson iterations of a nonlinear BVP
% using alternating Schwarz and Dirichlet BCs.

% So far no tested functions have displayed any kind of period doubling.
% These include: exponentials, logarithms, polynomials, power functions

%% Implementing test

% Problem parameters
C = 0.15;
F = @(x) tan(C*x);
Fp= @(x) C*sec(C*x).^2;
P = 1;

% Grid
nx = 1001;
a  = -0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
dx = 1./abs(x(2) - x(1));
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);
l1 = length(x1); l2 = length(x2);

% Diff. mat.
X1 = [ [x1(2:end) ; 0], x1, [0 ; x1(1:end-1)] ];
X2 = [ [x2(2:end) ; 0], x2, [0 ; x2(1:end-1)] ];

I1 = eye(l1);
D1 = ones(l1,1);
D1 = dx^2 * [ D1, -2*D1, D1 ];% - dx/2 * [ -D1, zeros(l1,1), D1 ].*X1;
D1 = spdiags(D1,[-1,0,1],l1,l1);% + spdiags(ones(l1,1),0,l1,l1);
I2 = eye(l2);
D2 = ones(l2,1);
D2 = dx^2 * [ D2, -2*D2, D2 ];% - dx/2 * [ -D2, zeros(l2,1), D2 ].*X2;
D2 = spdiags(D2,[-1,0,1],l2,l2);% + spdiags(ones(l2,1),0,l2,l2);

% Initialization
k = 0;
tol = 1e-8;
itermax = P+1;
G = zeros(1,101);
Gp= G;
itersaveNewton = G;
itersaveReg = G;
testu2 = linspace(-10,10,101);
% testu2=0.03;
nonlinsolves = 50;
for u2b0 = testu2
% u2b0 = 0;
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = 0*ones(size(x1)); u1(1)  = 0;
    u2 = 0*ones(size(x2)); u2(end)= 0;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - F(u1(2:end-1));
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        u1a = u1(x1==a);
        
        % Step 2: solve u in second domain
        u2(1) = u1a;
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2(2:end-1,:)*u2 - F(u2(2:end-1));
            u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
        end
        u2b = u2(x2==b);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(Fp(u1),0,l1,l1);
        g1 = dF1(2:end-1,2:end-1) \ ( -dF1(2:end-1,end) );
        g1 = [0 ; g1 ; 1 ];
        g1a= g1(x1==a);

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(Fp(u2),0,l2,l2);
        g2 = dF2(2:end-1,2:end-1) \ ( -g1a*dF2(2:end-1,1) );
        g2 = [ g1a ; g2 ; 0];
        g2b= g2(x2==b);

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)/(g2b - 1);
        
        if iter==P
            Gp(k) = u2b;
        end

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
    end
    
    if iter < itermax && error < tol
        itersaveNewton(k) = iter;
    else
        itersaveNewton(k) = NaN;
    end
        
end

k=0;
for u2b0 = testu2
% u2b0 = 0;
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = 0*ones(size(x1)); u1(1)  = 0;
    u2 = 0*ones(size(x2)); u2(end)= 0;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - F(u1(2:end-1));
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        u1a = u1(x1==a);
        
        % Step 2: solve u in second domain
        u2(1) = u1a;
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2(2:end-1,:)*u2 - F(u2(2:end-1));
            u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
        end
        u2b = u2(x2==b);
        
        if iter==1
            G(k) = u2b;
        end

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
    end
    
    if iter < itermax && error < tol
        itersaveReg(k) = iter;
    else
        itersaveReg(k) = NaN;
    end
        
end

figure(2)
plot(testu2,G,'r',testu2,Gp,'k',testu2,testu2,testu2,-testu2,'linewidth',2)
xlabel('\gamma')
ylabel('G(\gamma)')
legend('FP','NR')
axis([-10,10,-10,10])
set(gca,'fontsize',26,'linewidth',2)