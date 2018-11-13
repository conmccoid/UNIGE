%% Experiment 5: Schwarz method as a fixed point iteration

%% u_xx = u(1-u) on [0,1], u(0) = 0, u(1) = 1

% Problem parameters
f = @(x) x.*(1-x);

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
u1 = zeros(size(x1));
u2 = zeros(size(x2)); u2(x2==b) = 0.5;

for k = 1:20
    u2b(k) = u2(x2==b);
    u1a(k) = u1(x1==a);
    for j = 1:20
        f1 = f(u1); f1(1) = 0; f1(end) = u2b(k);
        u1 = U1 \ (L1 \ f1);
        f2 = f(u2); f2(end) = 1; f2(1) = u1a(k);
        u2 = U2 \ (L2 \ f2);
    end
    
    plot(x1,u1,'r--',x2,u2,'b--')
    pause
end

u2b(k+1) = u2(x2==b);
u1a(k+1) = u1(x1==a);
semilogy(1:21,abs(u2b - b),'*--')

%% The process as viewed as functions acting from R to R

% Problem parameters
f = @(x) exp(x);

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
u2b = linspace(0,1,101);
% u1a = u2b;

for k = 1:101
    u1 = zeros(size(x1));
    u2 = zeros(size(x2));
    
%     % Searching for two-cycle
%     for j = 1:20
%         f1 = f(u1); f1(1) = 0; f1(end) = u2b(k);
%         u1 = U1 \ (L1 \ f1);
%         f2 = f(u2); f2(end) = 1; f2(1) = u1a(k);
%         u2 = U2 \ (L2 \ f2);
%     end

    % Searching for fixed point
    for j = 1:20
        f1 = f(u1); f1(1) = 0; f1(end) = u2b(k);
        u1 = U1 \ (L1 \ f1);
    end
    u1a(k) = u1(x1==a);
    for j = 1:20
        f2 = f(u2); f2(end) = 1; f2(1) = u1a(k);
        u2 = U2 \ (L2 \ f2);
    end        
    u1_save(k) = u1(x1==a);
    u2_save(k) = u2(x2==b);
    
    plot(x1,u1,'r--',x2,u2,'b--')
    pause(0.1)
end

% plot(u2b,u1_save,'r*',u1a,u2_save,'bo')
plot(u2b,u2_save,'bo',u2b,u2b,'r-')