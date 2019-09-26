%% Experiment 18 - 2D Lu-sin(au)=0 fixed point function
% Repeating Exp16 for a course mesh: 2 points in the y-direction allows us
% to visualize the regions of convergence and divergence, should they exist

%% Implementing test

% Problem parameters
% 2 pts: 7.703
% 3 pts: ?
C = 1;  a = 8.45;
F = @(x) C * sin(a*x);
Fp= @(x) C*a*cos(a*x);
P = 1;

% Grid
nx = 21; ny = 6;
a  =-0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
y  = linspace(-1,1,ny)';
h  = x(2) - x(1);
hy = y(2) - y(1);
dx = 1./h;
dy = 1./hy;
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
II = speye(ny); IIint = II; IIint(1) = 0; IIint(end) = 0;
DD = epsilon * ones(ny,1); % perturbed y derivatives
DD = dy^2 * [ DD, -2*DD, DD ];
DD = spdiags(DD,[-1,0,1],ny,ny);

% BC parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
c1 = 0; d1 = 1; BC1 = sparse([1,1,1],ind1,[-c1*dx, d1, c1*dx],1,l1);
c2 = 0; d2 = 1; BC2 = sparse([1,1,1],ind2,[-c2*dx, d2, c2*dx],1,l2);
BCL= zeros(1,ny-2);
BCR= zeros(1,ny-2);

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);
% DDBC = sparse([1,1,1,2,2,2],[1,2,3,ny-2,ny-1,ny],...
%     0.5 * dx * [-3, 4, -1, 1, -4, 3],2,ny);
DDBC = [II(1,:) ; II(end,:)];
D1 = [J1BC(1,:) ; D1(2:end-1,:) ; J1BC(2,:)];
D2 = [J2BC(1,:) ; D2(2:end-1,:) ; J2BC(2,:)];
DD = [DDBC(1,:) ; DD(2:end-1,:) ; DDBC(2,:)];

% Laplacian
Lap1 = kron(IIint,D1) + kron(DD,I1int) + kron((II-IIint)*DD,I1-I1int);
Lap2 = kron(IIint,D2) + kron(DD,I2int) + kron((II-IIint)*DD,I2-I2int);
ind1x= kron(IIint,I1-I1int)* ones(ny*l1,1) == 1;
ind1y= kron(II-IIint,   I1)* ones(ny*l1,1) == 1;
ind2x= kron(IIint,I2-I2int)* ones(ny*l2,1) == 1;
ind2y= kron(II-IIint,   I2)* ones(ny*l2,1) == 1;

% Initialization
itermax = P+1;
L = 20; N = 101;
% testu2 = linspace(10,L,N);
% testu2 = linspace(-16.645,-16.63,N);
testu2x= linspace(10,13,N);
testu2y= linspace(37,40,N);
G = zeros(N,N);
Gp= zeros(N,N,P);
UNR = G;
VNR = G;
UFP = G;
VFP = G;
nonlinsolves  =10;
% G1 = repmat(G(:),1,P);
% G2 = G1;

for j = 1:N
    for k = 1:N
        fx     = [testu2x(k), testu2y(j), testu2y(j) testu2x(k)]; %try some kind of cosine function
        u2b    = fx;
        u2bold = u2b;
        iter   = 1;
        u1 = zeros(l1,ny); g1 = u1;
        u2 = zeros(l2,ny); g2 = u2;
        while iter < itermax

            % Step 1: solve u in first domain
            for i = 1:nonlinsolves
                J1        = Fp(u1(:));
                J1(ind1x) = 0; J1(ind1y) = 0;
                J1        = Lap1 - spdiags(J1,0,l1*ny,l1*ny);
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
                J2        = Lap2 - spdiags(J2,0,l2*ny,l2*ny);
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
            dF1= Lap1 - spdiags(dF1,0,l1*ny,l1*ny);
            g1 = dF1 \ kron(IIint,[ zeros(l1-1,1); 1]);
            g1a= kron(IIint,BC1) * g1;

            % Step 4: solve g in second domain
            dF2= Fp(u2(:)); dF2(ind2x) = 0; dF2(ind2y) = 0;
            dF2= Lap2 - spdiags(dF2,0,l2*ny,l2*ny);
            g2 = dF2 \ kron(IIint,[ 1; zeros(l2-1,1)]);
            g2b= kron(IIint,BC2) * g2; g2b = g2b * g1a - II; g2b = g2b(2:end-1,2:end-1);

            % Step 5: update u2b
            u2b = u2bold - ( g2b \ (u2b' - u2bold') )';
            
            Gp(k,j,iter) = norm(u2b)/norm(fx);
            if iter==P
                UNR(j,k) = u2b(1)-fx(1);
                VNR(j,k) = u2b(2)-fx(2);
            end
%             if Gp(k,j,iter)>=0.5
%                 plot(y,[0, u2b, 0],y,[0,-fx,0])
%                 pause(5)
%             end

            u2bold= u2b;
            iter  = iter+1;
            
%             if mod(iter,2)==1 && norm(u2b)/norm(fx)>=1
%                 G1(k+N*(j-1),floor(iter/2)) = fx(1);
%                 G2(k+N*(j-1),floor(iter/2)) = fx(2);
%             end

        end

    end
end

% testu2x = repmat(testu2,N,1);
% testu2y = repmat(testu2',1,N);
% figure(1)
% plot(testu2x(:),testu2y(:),'.',G1(:,1),G2(:,1),'.',G1(:,2),G2(:,2),'.',...
%     G1(:,3),G2(:,3),'.',G1(:,4),G2(:,4),'.',G1(:,5),G2(:,5),'.',G1(:,6),G2(:,6),'.')
% axis([-16.645,-16.63,-16.645,-16.63])
% set(gca,'fontsize',26,'linewidth',2)

for j = 1:N
    for k = 1:N
        fx     = [testu2x(k), testu2y(j), testu2y(j), testu2x(k)];
        u2b    = fx;
        u2bold = u2b;
        iter   = 1;
        u1 = zeros(l1,ny); g1 = u1;
        u2 = zeros(l2,ny); g2 = u2;
        while iter < itermax

            % Step 1: solve u in first domain
            for i = 1:nonlinsolves
                J1        = Fp(u1(:));
                J1(ind1x) = 0; J1(ind1y) = 0;
                J1        = Lap1 - spdiags(J1,0,l1*ny,l1*ny);
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
                J2        = Lap2 - spdiags(J2,0,l2*ny,l2*ny);
                F2        = Lap2*u2(:) - F(u2(:));
                BCy2      =(DDBC * u2' - pert)';
                BCx2      = J2BC * u2(:,2:end-1) - [ u1a ; BCR ];
                F2(ind2x) = BCx2(:); F2(ind2y) = BCy2(:);
                u2(:)     = u2(:) - J2 \ F2;
            end
            u2b = BC2 * u2; u2b = u2b(2:end-1);

            if iter==P
                G(k,j) = norm(u2b)/norm(fx);
                UFP(j,k) = u2b(1)-fx(1);
                VFP(j,k) = u2b(2)-fx(2);
            end

            u2bold= u2b;
            iter  = iter+1;

        end

    end
end

%%
figure(1)
Mat1 = repmat(testu2x,N,1);
Mat2 = repmat(testu2y',1,N);
contourf(testu2x,testu2y,G',0:0.5:1)
hold on
quiver(Mat1,Mat2,UFP,VFP,'linewidth',2)
hold off
xlabel('\gamma_1')
ylabel('\gamma_2')
title('FP')
axis square
set(gca,'fontsize',26,'linewidth',2)

figure(2)
contourf(testu2x,testu2y,Gp(:,:,end)',0:0.1:2)
hold on
quiver(Mat1,Mat2,UNR,VNR,'linewidth',2)
hold off
xlabel('\gamma_1')
ylabel('\gamma_2')
title('NR')
axis square
set(gca,'fontsize',26,'linewidth',2)

figure(3)
contourf(testu2x,testu2y,(2*UNR.^2 + VNR.^2)./(2*(Mat1.^2) + Mat2.^2),0:0.1:2,'linewidth',2)
xlabel('\gamma_1')
ylabel('\gamma_2')
axis square
set(gca,'fontsize',26,'linewidth',2)