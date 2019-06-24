%% Experiment 19 - modifying the nonlinear portion of the equation
% Rerun experiment 13 with a modulated sinusoid as the nonlinear function
% F. The hope is that this change will allow us to see how G will be
% affected.

%% Modulated function
a = 3.7;
% F = @(x) sin(x);
% Fp= @(x) cos(x);
% for k = 1:5
%     F = @(x) F(F(x));
%     Fp= @(x) Fp(x).*Fp(F(x));
% end
% F = @(x) F(a*x);
% Fp= @(x) a*Fp(a*x);

n = 9;
F = @(x) sin(a*x).^n;
Fp= @(x) n*a*cos(a*x) .* sin(a*x).^(n-1);

%% Implementing test

% Problem parameters
P = 1;

% Grid
nx = 1001;
a  =-0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
h  = x(2) - x(1);
dx = 1./h;
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

% BC parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
c1 = 0; d1 = 1; BCb = [-c1*dx, d1, c1*dx];
c2 = 0; d2 = 1; BCa = [-c2*dx, d2, c2*dx];
BCL= 0;
BCR= 0;

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);

% Initialization
N = 101;
k = 0;
tol = 1e-8;
itermax = P+1;
G = zeros(1,N);
Gp= G;
itersaveNewton = G;
itersaveReg = G;
L = 4;
testu2 = linspace(-L,L,N);
% testu2 = linspace(-1.8,-1.4,1001);
% testu2=0.03;
nonlinsolves = 10;
for u2b0 = testu2
% u2b0 = 0;
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = 0*ones(size(x1));
    u2 = 0*ones(size(x2));
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1*u1 - F(u1);
            J1 = [J1BC(1,:);J1(2:end-1,:);J1BC(2,:)];
            F1(1) = J1BC(1,:)*u1 - BCL; F1(end) = J1BC(2,:)*u1 - u2b;
            u1 = u1 - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=h);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2*u2 - F(u2);
            J2 = [J2BC(1,:);J2(2:end-1,:);J2BC(2,:)];
            F2(1) = J2BC(1,:)*u2 - u1a; F2(end) = J2BC(2,:)*u2 - BCR;
            u2 = u2 - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=h);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(Fp(u1),0,l1,l1);
        dF1= [J1BC(1,:);dF1(2:end-1,:);J1BC(2,:)];
        g1 = dF1 \ [ zeros(l1-1,1); 1];
        g1a= BCa * g1(abs(x1-a)<=h);

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(Fp(u2),0,l2,l2);
        dF2= [J2BC(1,:);dF2(2:end-1,:);J2BC(2,:)];
        g2 = dF2 \ [ g1a; zeros(l2-1,1)];
        g2b= BCb * g2(abs(x2-b)<=h);

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
    u1 = 0*ones(size(x1));
    u2 = 0*ones(size(x2));
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1*u1 - F(u1);
            J1 = [J1BC(1,:);J1(2:end-1,:);J1BC(2,:)];
            F1(1) = J1BC(1,:)*u1 - BCL; F1(end) = J1BC(2,:)*u1 - u2b;
            u1 = u1 - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=h);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2*u2 - F(u2);
            J2 = [J2BC(1,:);J2(2:end-1,:);J2BC(2,:)];
            F2(1) = J2BC(1,:)*u2 - u1a; F2(end) = J2BC(2,:)*u2 - BCR;
            u2 = u2 - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=h);
        
        if iter==1
            G(k) = u2b;
        end

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
        plot(x1,u1,'r',x2,u2,'b')
        axis([-1,1,-4,4])
        pause(0.1)
        
    end
    
    if iter < itermax && error < tol
        itersaveReg(k) = iter;
    else
        itersaveReg(k) = NaN;
    end
        
end

yy = linspace(-L,L,101);
% C  = -10:1:10;
% yC = bsxfun(@times,C',sqrt(abs(yy))); yC = bsxfun(@plus,yy,yC);
% xC = bsxfun(@plus,C',yy);
% 
% figure(1)
% plot(testu2,G,'r',testu2,Gp,'b',yy,yC,'k',yy,xC,'k','linewidth',2)
% xlabel('\gamma')
% ylabel('G(\gamma)')
% legend('FP','NR')
% % axis([-L,L,-L,L])
% set(gca,'fontsize',26,'linewidth',2)
%%
figure(2)
plot(testu2,G,'r.',testu2,Gp,'b.',yy,yy,'k--',yy,-yy,'k--','linewidth',2)
xlabel('\gamma')
ylabel('G(\gamma)')
legend('FP','NP')
axis equal
axis([-L,L,-L,L])
% axis([-1.8,-1.4,-1.8,-1.4])
set(gca,'fontsize',20,'linewidth',2)