%% Experiment 11 - Testing G(y) for Burgers
% Using the equation u''(x) - 2 u u' = 0 with homogeneous BCs (so u=0) we
% repeat the experiment seeking to establish which types of G functions are
% allowable.
% Note that since u=0 is the fixed point the problem reduces to G''=0 in
% the g-space. This suggests we'll get exactly the same regions we found
% for the simple u''=0 setup.
% We'll use a symmetric domain: a = -b = -1, alpha = -beta = -0.2.

% Problem appears to be singular, especially for choices of D1, gamma < 0.

%% D1 as a function of gamma and C1

% Parameters
a = -1;
b =  1;
al=-0.2;
be= 0.2;

% Grid
y = linspace(-1,1,1001);
C1= linspace(0.01 ,(pi/(2*(be-a)))^2,1001);
C1= C1';

% D1
D1= bsxfun(@minus,y,C1.*sec(sqrt(C1)*(be-a)).^2);
D1= bsxfun(@rdivide,D1,sqrt(C1).*tan(sqrt(C1)*(be-a)));
contourf(y,C1,D1,-10:10)

%% Mapping G(y)

% Grid
nx = 1001;
a  =-0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
h  = x(2) - x(1);
x1 = x(x<=b); b = x1(end); l1 = length(x1);
x2 = x(x>=a); a = x2(1);   l2 = length(x2);

% Differentiation matrices
I1 = eye(l1);
I2 = eye(l2);

D11= ones(l1,1)/(2*h);
D11= [-D11 , zeros(l1,1) , D11];
D11= spdiags(D11,[-1,0,1],l1,l1);
D12= ones(l2,1)/(2*h);
D12= [-D12 , zeros(l2,1) , D12];
D12= spdiags(D12,[-1,0,1],l2,l2);

D21= ones(l1,1)/h^2;
D21= [D21, -2*D21, D21];
D21= spdiags(D21,[-1,0,1],l1,l1);
D22= ones(l2,1)/h^2;
D22= [D22, -2*D22, D22];
D22= spdiags(D22,[-1,0,1],l2,l2);

% Problem parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
% c1 = 1; d1 =-5/6 + 1e-9;
% c2 = 1; d2 =-5/4 + 1e-8;
c1 = 1; d1 = 1;
c2 = 1; d2 = 1;

A  = 3.6;
BCL= 0;
BCR= 0;

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);

% Initialization
yy = -2:0.01:2;
nonlinsolves = 50;
G  = zeros(size(yy));
Gp = G;
N  = G;

% Map
for k = 1:length(yy)
    u1 = ones(size(x1));
    u2 = ones(size(x2));
    
    % 1st domain
    for i=1:nonlinsolves
        F1 = D21*u1 - D11*(u1.^2);
        J1 = D21 - 2*spdiags(u1,0,l1,l1)*D11 - 2*spdiags(D11*u1,0,l1,l1);
        J1 = [J1BC(1,:);J1(2:end-1,:);J1BC(2,:)];
        F1(1) = J1BC(1,:)*u1 - BCL; F1(end) = J1BC(2,:)*u1 - yy(k);
        u1 = u1 - J1 \ F1;
    end
    u1a= c2*D11(x1==a,:)*u1 + d2*u1(x1==a);
    g1 = J1 \ [ zeros(l1-1,1); 1];
    g1a= c2*D11(x1==a,:)*g1 + d2*g1(x1==a);
    
    % 2nd domain
    for i=1:nonlinsolves
        F2 = D22*u2 - D12*(u2.^2);
        J2 = D22 - 2*spdiags(u2,0,l2,l2)*D12 - 2*spdiags(D12*u2,0,l2,l2);
        J2 = [J2BC(1,:);J2(2:end-1,:);J2BC(2,:)];
        F2(1) = J2BC(1,:)*u2 - u1a; F2(end) = J2BC(2,:)*u2 - BCR;
        u2 = u2 - J2 \ F2;
    end
    u2b= c1*D12(x2==b,:)*u2 + d1*u2(x2==b);
    g2 = J2 \ [ g1a; zeros(l2-1,1)];
    g2b= c1*D12(x2==b,:)*g2 + d1*g2(x2==b);
    
    G(k) = u2b;
    Gp(k)= g2b;
    N(k) = yy(k) - (u2b - yy(k))/(g2b - 1);
    
end

figure(1)
plot(yy,G,'b',yy,N,'k',yy,yy,'--',yy,-yy,'--')%,'Linewidth',2)
axis([-2,2,-2,2])
axis square
xlabel('\gamma')
ylabel('G(\gamma)')
legend('G(\gamma)','NR')
% set(gca,'linewidth',2,'fontsize',26)

xx = -2:0.01:2;
C  = -10:0.5:10;
yC = bsxfun(@times,C',sqrt(abs(xx))); yC = bsxfun(@plus,xx,yC);

figure(2)
plot(yy,G,'b',xx,yC,'k')%,'linewidth',2)
axis([-2,2,-2,2])
axis square
xlabel('\gamma')
ylabel('G(\gamma)')
% set(gca,'linewidth',2,'fontsize',26)