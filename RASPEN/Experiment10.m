%% Experiment 10 - Cycle for ASPN
% Using u''(x) - sin(Cu(x)) = 0 as our problem, we enter into the theorized
% cycle(s) for DD preconditioned with Newton-Raphson.
% For C<=3 there are no cycles. For C>=4 the cycles are all unstable
% (divergent? what is the term here?). I have found C=3.6 to have stable
% cycles for initial guess in the region of 1.6.

% Problem parameters
C = 3.6;

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

% Initialization for ASPN
tol = 1e-16;
itermax = 100;
% u2b0 = -0.606336103110884206; % closest to the fixed point for C=8
u2b0 = -1.6;
nonlinsolves = 50;
u2b    = u2b0;
u2bold = u2b0;
u1 = 0*ones(size(x1)); u1(1)  = 0;
u2 = 0*ones(size(x2)); u2(end)= 0;

U2b    = u2b0;
U2bold = u2b0;
error  = 1;
iter   = 1;
U1 = 0*ones(size(x1)); U1(1)  = 0;
U2 = 0*ones(size(x2)); U2(end)= 0;

% ASPN
while error > tol && iter < itermax

    % Step 1: solve u in first domain
    u1(end) = u2b;
    for i = 1:nonlinsolves
        J1 = D1 - spdiags(C*cos(C*u1),0,l1,l1);
        F1 = D1(2:end-1,:)*u1 - sin(C*u1(2:end-1));
        u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
    end
    u1a = u1(x1==a);

    % Step 2: solve u in second domain
    u2(1) = u1a;
    for i = 1:nonlinsolves
        J2 = D2 - spdiags(C*cos(C*u2),0,l2,l2);
        F2 = D2(2:end-1,:)*u2 - sin(C*u2(2:end-1));
        u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
    end
    u2b = u2(x2==b);

    % Preconditioning with Newton
    % Step 3: solve g in first domain
    dF1= D1 - spdiags(C*cos(C*u1),0,l1,l1);
    g1 = dF1(2:end-1,2:end-1) \ ( -dF1(2:end-1,end) );
    g1 = [0 ; g1 ; 1 ];
    g1a= g1(x1==a);

    % Step 4: solve g in second domain
    dF2= D2 - spdiags(C*cos(C*u2),0,l2,l2);
    g2 = dF2(2:end-1,2:end-1) \ ( -g1a*dF2(2:end-1,1) );
    g2 = [ g1a ; g2 ; 0];
    g2b= g2(x2==b);

    % Step 5: update u2b
    u2b = u2bold - (u2b - u2bold)/(g2b - 1);
    u2bold = u2b;
    
    % Step 1: solve u in first domain
    U1(end) = U2b;
    for i = 1:nonlinsolves
        J1 = D1 - spdiags(C*cos(C*U1),0,l1,l1);
        F1 = D1(2:end-1,:)*U1 - sin(C*U1(2:end-1));
        U1(2:end-1) = U1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
    end
    U1a = U1(x1==a);

    % Step 2: solve u in second domain
    U2(1) = U1a;
    for i = 1:nonlinsolves
        J2 = D2 - spdiags(C*cos(C*U2),0,l2,l2);
        F2 = D2(2:end-1,:)*U2 - sin(C*U2(2:end-1));
        U2(2:end-1) = U2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
    end
    U2b = U2(x2==b);

    error = abs(U2b - U2bold);
    U2bold= U2b;
    iter  = iter+1;
    
    figure(2)
    subplot(1,2,1)
    plot(x1,u1,x2,u2,'linewidth',2)
    xlabel('x')
    ylabel('u(x)')
    axis([-1,1,-2,2])
    title('NR')
    set(gca,'fontsize',26,'linewidth',2)
    
    subplot(1,2,2)
    plot(x1,U1,x2,U2,'linewidth',2)
    xlabel('x')
    ylabel('u(x)')
    axis([-1,1,-2,2])
    title('FP')
    set(gca,'fontsize',26,'linewidth',2)
    pause(0.5)

end