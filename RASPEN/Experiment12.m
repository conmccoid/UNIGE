%% Experiment 12 - Counterexamples
% We collect a series of counterexamples to be used to dispute convergence
% conjectures for FP and NR in transmission conditions of alternating
% Schwarz methods.

%% u'' - Csin(u)=0 with homogeneous Dirichlet BCs

% Intervals
a = -1;
b =  1;
al=-0.2;
be= 0.2;

% Parameter(s)
C = linspace(0,2,1001);
% C = sqrt(0.9) * pi / 2;
% C = 1;
% C1=linspace(-2,2,101);
% C2=C1';
p = linspace(-10,10,101)';

% Possible regions
% Gp= ( (sinh(C*al).*cosh(C*a) - cosh(C*al).*sinh(C*a)) .* ...
%     (sinh(C*be).*cosh(C*b) - cosh(C*be).*sinh(C*b)) ) ./ ...
%     ( (sinh(C*be).*cosh(C*a) - cosh(C*be).*sinh(C*a)) .* ...
%     (sinh(C*al).*cosh(C*b) - cosh(C*al).*sinh(C*b)) );
% Gp= ( (sin(C*al).*cos(C*a) - cos(C*al).*sin(C*a)) .* ...
%     (sin(C*be).*cos(C*b) - cos(C*be).*sin(C*b)) ) ./ ...
%     ( (sin(C*be).*cos(C*a) - cos(C*be).*sin(C*a)) .* ...
%     (sin(C*al).*cos(C*b) - cos(C*al).*sin(C*b)) );
% semilogy(C,Gp)
% xlabel('C^{1/2}')
% ylabel('G''(\gamma)')

% Gp= bsxfun(@times, C*cosh(C*(al-a)) + C2*sinh(C*(al-a)), C*cosh(C*(be-b)) + C1*sinh(C*(be-b)));
% Gp= bsxfun(@rdivide, Gp, C*cosh(C*(be-a)) + C1*sinh(C*(be-a)));
% Gp= bsxfun(@rdivide, Gp, C*cosh(C*(al-b)) + C2*sinh(C*(al-b)));
k1= bsxfun(@times, p, sinh(C*(al-a)));
k2= bsxfun(@times, p, sinh(C*(be-b)));
k3= bsxfun(@times, p, sinh(C*(be-a)));
k4= bsxfun(@times, p, sinh(C*(al-b)));
Gp= bsxfun(@times, C.*cosh(C*(al-a)) - k1, C.*cosh(C*(be-b)) + k2);
Gp= bsxfun(@rdivide, Gp, C.*cosh(C*(be-a)) + k3);
Gp= bsxfun(@rdivide, Gp, C.*cosh(C*(al-b)) - k4);

contourf(C,p,Gp,[-1,0,1])
xlabel('C')
ylabel('p')
% contourf(C1,C2,Gp,[-1,0,1])
% xlabel('C_1')
% ylabel('C_2')

%%

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
c1 = 1; d1 =-2;
c2 = 1; d2 = 2;

A  = 1;
BCL= 0;
BCR= 0;

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);

% Initialization
L  = 1;
yy = linspace(-L,L,1001);
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
%         F1 = D21*u1 - u1 .* (D11 * u1);
%         J1 = D21 - spdiags(u1,0,l1,l1)*D11 - spdiags(D11*u1,0,l1,l1);
%         F1 = D21*u1 - u1.^A;
%         J1 = D21 - spdiags(A*u1.^(A-1),0,l1,l1);
        F1 = D21*u1 - A*sin(u1);
        J1 = D21 - spdiags(A*cos(u1),0,l1,l1);
        J1 = [J1BC(1,:);J1(2:end-1,:);J1BC(2,:)];
        F1(1) = J1BC(1,:)*u1 - BCL; F1(end) = J1BC(2,:)*u1 - yy(k);
        u1 = u1 - J1 \ F1;
    end
    u1a= c2*D11(x1==a,:)*u1 + d2*u1(x1==a);
    g1 = J1 \ [ zeros(l1-1,1); 1];
    g1a= c2*D11(x1==a,:)*g1 + d2*g1(x1==a);
    
    % 2nd domain
    for i=1:nonlinsolves
%         F2 = D22*u2 - u2 .* (D12*u2);
%         J2 = D22 - spdiags(u2,0,l2,l2)*D12 - spdiags(D12*u2,0,l2,l2);
%         F2 = D22*u2 - u2.^A;
%         J2 = D22 - spdiags(A*u2.^(A-1),0,l2,l2);
        F2 = D22*u2 - A*sin(u2);
        J2 = D22 - spdiags(A*cos(u2),0,l2,l2);
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
    
%     plot(x1,u1,x2,u2)
%     axis([-1,1,-1,1])
%     pause(0.01)
    
end

figure(1)
plot(yy,G,'b',yy,N,'k',yy,yy,'--',yy,-yy,'--','Linewidth',2)
axis([-L,L,-L,L])
axis square
xlabel('\gamma')
ylabel('G(\gamma)')
legend('G(\gamma)','NR')
set(gca,'linewidth',2,'fontsize',26)

C  = -10:0.5:10;
yC = bsxfun(@times,C',sqrt(abs(yy))); yC = bsxfun(@plus,yy,yC);

figure(2)
plot(yy,G,'b',yy,yC,'k')%,'linewidth',2)
axis([-L,L,-L,L])
axis square
xlabel('\gamma')
ylabel('G(\gamma)')