%% Experiment 16 - 2D Lu-sin(au)=0 fixed point function
% Repeating Exp13 in 2D

%% Implementing test

% Problem parameters
% 2-cycles: 3.6
% 4-cycles: 3.72
% 8-cycles: 3.731
%    Chaos: 3.735
C = 1;  a = 3.6;
F = @(x) C * sin(a*x);
Fp= @(x) C*a*cos(a*x);
P = 1;

% Grid
nx = 21;
a  =-0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
h  = x(2) - x(1);
dx = 1./h;
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);
l1 = length(x1); l2 = length(x2);
ind1 = 1:l1; ind1 = ind1(abs(x1-a)<=(h+2*eps));
ind2 = 1:l2; ind2 = ind2(abs(x2-b)<=(h+2*eps));

% Diff. mat.
X1 = [ [x1(2:end) ; 0], x1, [0 ; x1(1:end-1)] ];
X2 = [ [x2(2:end) ; 0], x2, [0 ; x2(1:end-1)] ];

I1 = speye(l1); I1int = I1; I1int(1) = 0; I1int(end) = 0;
D1 = ones(l1,1);
D1 = dx^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[-1,0,1],l1,l1);
I2 = speye(l2); I2int = I2; I2int(1) = 0; I2int(end) = 0;
D2 = ones(l2,1);
D2 = dx^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[-1,0,1],l2,l2);

pert    = 0;
epsilon = 1; %1e-3 gives period doubling
II = speye(nx); IIint = II; IIint(1) = 0; IIint(end) = 0;
DD = epsilon * ones(nx,1); % perturbed y derivatives
DD = dx^2 * [ DD, -2*DD, DD ];
DD = spdiags(DD,[-1,0,1],nx,nx);

% BC parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
c1 = 0; d1 = 1; BC1 = sparse([1,1,1],ind1,[-c1*dx, d1, c1*dx],1,l1);
c2 = 0; d2 = 1; BC2 = sparse([1,1,1],ind2,[-c2*dx, d2, c2*dx],1,l2);
BCL= zeros(1,nx-2);
BCR= zeros(1,nx-2);

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);
% DDBC = sparse([1,1,1,2,2,2],[1,2,3,nx-2,nx-1,nx],...
%     0.5 * dx * [-3, 4, -1, 1, -4, 3],2,nx);
DDBC = [II(1,:) ; II(end,:)];
D1 = [J1BC(1,:) ; D1(2:end-1,:) ; J1BC(2,:)];
D2 = [J2BC(1,:) ; D2(2:end-1,:) ; J2BC(2,:)];
DD = [DDBC(1,:) ; DD(2:end-1,:) ; DDBC(2,:)];

% Laplacian
Lap1 = kron(IIint,D1) + kron(DD,I1int) + kron((II-IIint)*DD,I1-I1int);
Lap2 = kron(IIint,D2) + kron(DD,I2int) + kron((II-IIint)*DD,I2-I2int);
ind1x= kron(IIint,I1-I1int)* ones(nx*l1,1) == 1;
ind1y= kron(II-IIint,   I1)* ones(nx*l1,1) == 1;
ind2x= kron(IIint,I2-I2int)* ones(nx*l2,1) == 1;
ind2y= kron(II-IIint,   I2)* ones(nx*l2,1) == 1;

% Initialization
tol = 1e-8;
itermax = P+1;
% fx= 1-x(2:end-1)'.^2;
fx= ones(1,nx-2);
L = 10/norm(fx); N = 101;
testu2 = linspace(-L,L,N);
G = zeros(N,1);
Gp= G;
itersaveNewton= G;
itersaveReg   = G;
nonlinsolves  =10;
for k = 1:N
    u2b    = testu2(k)*fx;
    u2bold = u2b;
    error  = 1;
    iter   = 1;
    u1 = zeros(l1,nx); g1 = u1;
    u2 = zeros(l2,nx); g2 = u2;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1        = Fp(u1(:));
            J1(ind1x) = 0; J1(ind1y) = 0;
            J1        = Lap1 - spdiags(J1,0,l1*nx,l1*nx);
            F1        = Lap1*u1(:) - F(u1(:));
            BCy1      =(DDBC * u1'- pert)';
            BCx1      = J1BC * u1(:,2:end-1) - [ BCL ; u2b ];
            F1(ind1x) = BCx1(:); F1(ind1y) = BCy1(:);
            u1(:)     = u1(:) - J1 \ F1;
        end
        u1a= BC1 * u1; u1a = u1a(2:end-1);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2        = Fp(u2(:));
            J2(ind2x) = 0; J2(ind2y) = 0;
            J2        = Lap2 - spdiags(J2,0,l2*nx,l2*nx);
            F2        = Lap2*u2(:) - F(u2(:));
            BCy2      =(DDBC * u2' - pert)';
            BCx2      = J2BC * u2(:,2:end-1) - [ u1a ; BCR ];
            F2(ind2x) = BCx2(:); F2(ind2y) = BCy2(:);
            u2(:)     = u2(:) - J2 \ F2;
        end
        u2b = BC2 * u2; u2b = u2b(2:end-1);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= Fp(u1(:)); dF1(ind1x) = 0; dF1(ind1y) = 0;
        dF1= Lap1 - spdiags(dF1,0,l1*nx,l1*nx);
        g1 = dF1 \ kron(IIint,[ zeros(l1-1,1); 1]);
        g1a= kron(IIint,BC1) * g1;

        % Step 4: solve g in second domain
        dF2= Fp(u2(:)); dF2(ind2x) = 0; dF2(ind2y) = 0;
        dF2= Lap2 - spdiags(dF2,0,l2*nx,l2*nx);
        g2 = dF2 \ kron(IIint,[ 1; zeros(l2-1,1)]);
        g2b= kron(IIint,BC2) * g2; g2b = g2b * g1a - II; g2b = g2b(2:end-1,2:end-1);

        % Step 5: update u2b
        u2b = u2bold - ( g2b \ (u2b' - u2bold') )';
        u2bold= u2b;
        
        if iter==P
            Gp(k) = norm(u2b);
        end

        error = norm(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
    end
    
    if iter < itermax && error < tol
        itersaveNewton(k) = iter;
    else
        itersaveNewton(k) = NaN;
    end
        
end

% k=0;
for k = 1:N
    u2b    = testu2(k)*fx;
    u2bold = u2b;
    error  = 1;
    iter   = 1;
    u1 = zeros(l1,nx); g1 = u1;
    u2 = zeros(l2,nx); g2 = u2;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1        = Fp(u1(:));
            J1(ind1x) = 0; J1(ind1y) = 0;
            J1        = Lap1 - spdiags(J1,0,l1*nx,l1*nx);
            F1        = Lap1*u1(:) - F(u1(:));
            BCy1      =(DDBC * u1'- pert)';
            BCx1      = J1BC * u1(:,2:end-1) - [ BCL ; u2b ];
            F1(ind1x) = BCx1(:); F1(ind1y) = BCy1(:);
            u1(:)     = u1(:) - J1 \ F1;
        end
        u1a= BC1 * u1; u1a = u1a(2:end-1);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2        = Fp(u2(:));
            J2(ind2x) = 0; J2(ind2y) = 0;
            J2        = Lap2 - spdiags(J2,0,l2*nx,l2*nx);
            F2        = Lap2*u2(:) - F(u2(:));
            BCy2      =(DDBC * u2' - pert)';
            BCx2      = J2BC * u2(:,2:end-1) - [ u1a ; BCR ];
            F2(ind2x) = BCx2(:); F2(ind2y) = BCy2(:);
            u2(:)     = u2(:) - J2 \ F2;
        end
        u2b = BC2 * u2; u2b = u2b(2:end-1);
        
        if iter==1
            G(k) = norm(u2b);
        end

        error = norm(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
    end
    
    if iter < itermax && error < tol
        itersaveReg(k) = iter;
    else
        itersaveReg(k) = NaN;
    end
        
end

L  = norm(fx) * L;
testu2 = norm(fx) * testu2;
root = mean(testu2(abs(testu2-G')<0.1));

%%
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

figure(2)
% plot(testu2,G(:,round(nx/2)-1),'r.',testu2,Gp(:,round(nx/2)-1),'k.',yy,yy,yy,-yy,'linewidth',2)
plot(testu2,G,'r.',testu2,Gp,'k.',yy,yy,yy,2*root-yy,'linewidth',2)
xlabel('$$\Vert \gamma \Vert$$','interpreter','latex')
ylabel('$$\Vert G(\gamma) \Vert$$','interpreter','latex')
legend('FP','NR')
axis([-L,L,-L,L])
% axis([-1.8,-1.4,-1.8,-1.4])
set(gca,'fontsize',26,'linewidth',2)