%% Experiment 6 - Applying the new Newton precond. to SAM
% Same problem as experiment 5

% Problem parameters
f = @(x) x.*(x-1);
df = @(x) 2*x - 1;

% Grid
nx = 101; dx = 1/(nx-1);
a  = 1/3;
b  = 2/3;
x  = linspace(0,1,nx)';
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);

% Diff. mat.
I1 = eye(length(x1));
D1 = ones(length(x1)-2,1);
D1 = (nx-1)^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[0,1,2],length(x1)-2,length(x1));
D1 = [I1(1,:) ; D1 ; I1(end,:) ];
[L1,U1] = lu(D1);
I2 = eye(length(x2));
D2 = ones(length(x2)-2,1);
D2 = (nx-1)^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[0,1,2],length(x2)-2,length(x2));
D2 = [I2(1,:) ; D2 ; I2(end,:) ];
[L2,U2] = lu(D2);

% Initialization
k = 0;
tol = 1e-8;
itermax = 100;
itersave = zeros(1001,1);
testu2 = linspace(0,1,1001);
for u2b0 = testu2
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    while error > tol && iter < itermax
        u1 = zeros(size(x1));
        u2 = zeros(size(x2));

        % Step 1: solve u in first domain
        for j = 1:20
            f1 = f(u1); f1(1) = 0; f1(end) = u2b;
            u1 = U1 \ (L1 \ f1);
        end
        u1a = u1(x1==a);

        % Step 2: solve u in second domain
        for j = 1:20
            f2 = f(u2); f2(end) = 1; f2(1) = u1a;
            u2 = U2 \ (L2 \ f2);
        end
%         u2b = u2(x2==b);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1 = diag(df(u1));
        dF1([1;end],:) = zeros(2,length(x1));
        g1 = (D1 - dF1) \ [ zeros(length(x1)-1,1) ; 1 ];

        % Step 4: solve g in second domain
        dF2 = diag(df(u2));
        dF2([1;end],:) = zeros(2,length(x2));
        g2 = (D2 - dF2) \ [ g1(x1==a) ; zeros(length(x2)-1,1) ];
        g2b = g2(x2==b);

        % Step 5: update u2b
        u2b = ( g2b*u2b - u2(x2==b) )./ ( g2b - 1 );

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
%         semilogy(iter-1,error,'b*')
%         hold on
%         pause(0.1)
    end
    
    if iter < itermax && error < tol
        itersave(k) = iter;
    else
        itersave(k) = NaN;
    end
end

plot(testu2,itersave)
xlabel('\gamma')
ylabel('Number of iterations to convergence')
title('Newton precond. on transmission condition')

%% Showing order of convergence

% Problem parameters
f = @(x) x.*(x-1);
df = @(x) 2*x - 1;

% Grid
nx = 101; dx = 1/(nx-1);
a  = 1/3;
b  = 2/3;
x  = linspace(0,1,nx)';
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);

% Diff. mat.
I1 = eye(length(x1));
D1 = ones(length(x1)-2,1);
D1 = (nx-1)^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[0,1,2],length(x1)-2,length(x1));
D1 = [I1(1,:) ; D1 ; I1(end,:) ];
[L1,U1] = lu(D1);
I2 = eye(length(x2));
D2 = ones(length(x2)-2,1);
D2 = (nx-1)^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[0,1,2],length(x2)-2,length(x2));
D2 = [I2(1,:) ; D2 ; I2(end,:) ];
[L2,U2] = lu(D2);

% Initialization
itermax = 100;
tol = 1e-8;
error0 = zeros(itermax,1);
error0(1) = 1;
error1 = error0;

% Standard additive Schwarz
u2b = 0;
u2bold = 0;
iter = 0;
while error0(iter+1) > tol && iter < itermax
    iter = iter+1;
    u1 = zeros(size(x1));
    u2 = zeros(size(x2));

    % Step 1: solve u in first domain
    for j = 1:20
        f1 = f(u1); f1(1) = 0; f1(end) = u2b;
        u1 = U1 \ (L1 \ f1);
    end
    u1a = u1(x1==a);

    % Step 2: solve u in second domain
    for j = 1:20
        f2 = f(u2); f2(end) = 1; f2(1) = u1a;
        u2 = U2 \ (L2 \ f2);
    end
    
    % Step 3: update u2b
    u2b = u2(x2==b);

    error0(iter+1) = abs(u2b - u2bold);
    u2bold= u2b;

end

% Preconditioning with Newton
u2b = 0;
u2bold = 0;
iter = 0;
while error1(iter+1) > tol && iter < itermax
    iter = iter+1;
    u1 = zeros(size(x1));
    u2 = zeros(size(x2));

    % Step 1: solve u in first domain
    for j = 1:20
        f1 = f(u1); f1(1) = 0; f1(end) = u2b;
        u1 = U1 \ (L1 \ f1);
    end
    u1a = u1(x1==a);

    % Step 2: solve u in second domain
    for j = 1:20
        f2 = f(u2); f2(end) = 1; f2(1) = u1a;
        u2 = U2 \ (L2 \ f2);
    end

    % Step 3: solve g in first domain
    dF1 = diag(df(u1));
    dF1([1;end],:) = zeros(2,length(x1));
    g1 = (D1 - dF1) \ [ zeros(length(x1)-1,1) ; 1 ];

    % Step 4: solve g in second domain
    dF2 = diag(df(u2));
    dF2([1;end],:) = zeros(2,length(x2));
    g2 = (D2 - dF2) \ [ g1(x1==a) ; zeros(length(x2)-1,1) ];
    g2b = g2(x2==b);

    % Step 5: update u2b
    u2b = ( g2b*u2b - u2(x2==b) )./ ( g2b - 1 );

    error1(iter+1) = abs(u2b - u2bold);
    u2bold= u2b;

end

semilogy(1:itermax,error0,'b*--',1:itermax,error1,'r*--')
xlabel('Iteration')
ylabel('Change in u_2(\beta)')
legend('AS','ASPN')