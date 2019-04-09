%% Experiment 16 - 2D Lu-sin(au)=0 fixed point function
% Repeating Exp13 in 2D

%% Implementing test

% Problem parameters
% 2-cycles: 3.6
% 4-cycles: 3.72
% 8-cycles: 3.731
%    Chaos: 3.735
C = 1;  a = 4;
F = @(x) C * sin(a*x);
Fp= @(x) C*a*cos(a*x);
P = 1;

% Grid
nx = 51;
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

I1 = speye(l1); I1int = I1; I1int(1) = 0; I1int(end) = 0;
D1 = ones(l1,1);
D1 = dx^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[-1,0,1],l1,l1);
I2 = speye(l2); I2int = I2; I2int(1) = 0; I2int(end) = 0;
D2 = ones(l2,1);
D2 = dx^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[-1,0,1],l2,l2);

II = speye(nx); IIint = II; IIint(1) = 0; IIint(end) = 0;
DD = zeros(nx,1); % removed y derivatives
DD = dx^2 * [ DD, -2*DD, DD ];
DD = spdiags(DD,[-1,0,1],nx,nx);

% BC parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
c1 = 0; d1 = 1; BCb = [-c1*dx, d1, c1*dx];
c2 = 0; d2 = 1; BCa = [-c2*dx, d2, c2*dx];
BCL= zeros(1,nx);
BCR= zeros(1,nx);

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);
D1 = [J1BC(1,:) ; D1(2:end-1,:) ; J1BC(2,:)];
D2 = [J2BC(1,:) ; D2(2:end-1,:) ; J2BC(2,:)];
DD = [II(1,:)   ; DD(2:end-1,:) ; II(end,:)];

% Laplacian
Lap1 = kron(II,D1) + kron(DD,I1int);
Lap2 = kron(II,D2) + kron(DD,I2int);
ind1x= kron(II,I1-I1int)   * ones(nx*l1,1) == 1;
ind1y= kron(II-IIint,I1int)* ones(nx*l1,1) == 1;
ind2x= kron(II,I2-I2int)   * ones(nx*l2,1) == 1;
ind2y= kron(II-IIint,I2int)* ones(nx*l2,1) == 1;

% Initialization
tol = 1e-8;
itermax = P+1;
fx= 1-x'.^2;
L = 10/norm(fx); N = 1001;
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
            J1 = Fp(u1(:));
            J1(ind1x) = 0; J1(ind1y) = 0;
            J1 = Lap1 - spdiags(J1,0,l1*nx,l1*nx);
            F1 = Lap1*u1(:) - F(u1(:));
            F1(ind1y) = 0;
            BCx1 = J1BC * u1 - [ BCL ; u2b ];
            F1(ind1x) = BCx1(:);
            u1(:) = u1(:) - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=(h+2*eps),:);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = Fp(u2(:));
            J2(ind2x) = 0; J2(ind2y) = 0;
            J2 = Lap2 - spdiags(J2,0,l2*nx,l2*nx);
            F2 = Lap2*u2(:) - F(u2(:));
            F2(ind2y) = 0;
            BCx2 = J2BC * u2 - [ u1a ; BCR ];
            F2(ind2x) = BCx2(:);
            u2(:) = u2(:) - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=(h+2*eps),:);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= Fp(u1(:)); dF1(ind1x) = 0; dF1(ind1y) = 0;
        dF1= Lap1 - spdiags(dF1,0,l1*nx,l1*nx);
        g1(:) = dF1 \ kron([0 ; ones(nx-2,1) ; 0],[ zeros(l1-1,1); 1]);
        g1a= BCa * g1(abs(x1-a)<=(h+2*eps),:);

        % Step 4: solve g in second domain
        dF2= Fp(u2(:)); dF2(ind2x) = 0; dF2(ind2y) = 0;
        dF2= Lap2 - spdiags(dF2,0,l2*nx,l2*nx);
        g2(:) = dF2 \ kron([0 ; ones(nx-2,1) ; 0],[ 1; zeros(l2-1,1)]);
        g2b= BCb * g2(abs(x2-b)<=(h+2*eps),:); g2b = g2b .* g1a;

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)./(g2b - 1);
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
            J1 = Fp(u1(:));
            J1(ind1x) = 0; J1(ind1y) = 0;
            J1 = Lap1 - spdiags(J1,0,l1*nx,l1*nx);
            F1 = Lap1*u1(:) - F(u1(:));
            F1(ind1y) = 0;
            BCx1 = J1BC * u1 - [ BCL ; u2b ];
            F1(ind1x) = BCx1(:);
            u1(:) = u1(:) - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=(h+2*eps),:);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = Fp(u2(:));
            J2(ind2x) = 0; J2(ind2y) = 0;
            J2 = Lap2 - spdiags(J2,0,l2*nx,l2*nx);
            F2 = Lap2*u2(:) - F(u2(:));
            F2(ind2y) = 0;
            BCx2 = J2BC * u2 - [ u1a ; BCR ];
            F2(ind2x) = BCx2(:);
            u2(:) = u2(:) - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=(h+2*eps),:);
        
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

%%
L  = norm(fx) * L;
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
plot(testu2*norm(fx),G,'r.',testu2*norm(fx),Gp,'k.',yy,yy,yy,-yy,'linewidth',2)
xlabel('$$\Vert \gamma \Vert$$','interpreter','latex')
ylabel('$$\Vert G(\gamma) \Vert$$','interpreter','latex')
legend('FP','NR')
axis([-L,L,-L,L])
% axis([-1.8,-1.4,-1.8,-1.4])
set(gca,'fontsize',26,'linewidth',2)