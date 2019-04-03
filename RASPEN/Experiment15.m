%% Experiment 15 - plotting the bifurcation diagram for Lu-sin(au)=0
% L is the Laplacian, so we are checking the problem in the second
% dimension. We use homogeneous BCs on all sides of the square [-1,1] x
% [-1,1].

%% Bifurcation diagram

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
DD = ones(nx,1);
DD = dx^2 * [ DD, -2*DD, DD ];
DD = spdiags(DD,[-1,0,1],nx,nx);

% BC parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
c1 = 0; d1 = 1; BCb = [-c1*dx, d1, c1*dx];
c2 = 0; d2 = 1; BCa = [-c2*dx, d2, c2*dx];
BCL= zeros(1,nx-2);
BCR= zeros(1,nx-2);

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);
D1 = [J1BC(1,:) ; D1(2:end-1,:) ; J1BC(2,:)];
D2 = [J2BC(1,:) ; D2(2:end-1,:) ; J2BC(2,:)];
DD = [II(1,:)   ; DD(2:end-1,:) ; II(end,:)];

% Laplacian
Lap1 = kron(IIint,D1) + kron(DD,I1int) + kron(II-IIint,I1-I1int);
Lap2 = kron(IIint,D2) + kron(DD,I2int) + kron(II-IIint,I2-I2int);
ind1x= kron(IIint,I1-I1int) * ones(nx*l1,1) == 1;
ind1y= (kron(II-IIint,I1int) + kron(II-IIint,I1-I1int)) * ones(nx*l1,1) == 1;
ind2x= kron(IIint,I2-I2int) * ones(nx*l2,1) == 1;
ind2y= (kron(II-IIint,I2int) + kron(II-IIint,I1-I2int)) * ones(nx*l2,1) == 1;

% Plotting bifurcation region
strpt = 20.5;
endpt = 21.5;
res   = 0.1;
% strpt = 9.10;
% endpt = 9.155;
% res   = 0.001;

ind = round((endpt - strpt) / res);
gam = zeros(ind,64,nx-2);
nonlinsolves = 10;

for k = 1:ind
    C = strpt + k*res;
    F = @(x)   sin(C*x);
    Fp= @(x) C*cos(C*x);
    
    u1 = zeros(l1,nx); g1 = u1;
    u2 = zeros(l2,nx); g2 = u2;
    u2b= 1.6 * ones(1,nx-2);
%     u2b= 1.45;
    u2bold = u2b;
    for iter = 1:(50+64)
        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1 = Fp(u1(:));
            J1(ind1x) = 0; J1(ind1y) = 0;
            J1 = Lap1 - spdiags(J1,0,l1*nx,l1*nx);
            F1 = Lap1*u1(:) - F(u1(:));
            F1(ind1y) = 0;
            BCx1 = J1BC * u1(:,2:end-1) - [ BCL ; u2b ];
            F1(ind1x) = BCx1(:);
            u1(:) = u1(:) - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=(h+2*eps),2:end-1);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = Fp(u2(:));
            J2(ind2x) = 0; J2(ind2y) = 0;
            J2 = Lap2 - spdiags(J2,0,l2*nx,l2*nx);
            F2 = Lap2*u2(:) - F(u2(:));
            F2(ind2y) = 0;
            BCx2 = J2BC * u2(:,2:end-1) - [ u1a ; BCR ];
            F2(ind2x) = BCx2(:);
            u2(:) = u2(:) - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=(h+2*eps),2:end-1);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= Fp(u1(:)); dF1(ind1x) = 0; dF1(ind1y) = 0;
        dF1= Lap1 - spdiags(dF1,0,l1*nx,l1*nx);
        g1(:) = dF1 \ kron([0 ; ones(nx-2,1) ; 0],[ zeros(l1-1,1); 1]);
        g1a= BCa * g1(abs(x1-a)<=(h+2*eps),2:end-1);

        % Step 4: solve g in second domain
        dF2= Fp(u2(:)); dF2(ind2x) = 0; dF2(ind2y) = 0;
        dF2= Lap2 - spdiags(dF2,0,l2*nx,l2*nx);
        g2(:) = dF2 \ kron([0 ; ones(nx-2,1) ; 0],[ 1; zeros(l2-1,1)]);
        g2b= BCb * g2(abs(x2-b)<=(h+2*eps),2:end-1); g2b = g2b .* g1a;

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)./(g2b - 1);
        u2bold= u2b;
        
        if iter > 50
            gam(k,iter-50,:) = u2b;
        end
        
    end
end

gam_neg = zeros(size(gam));

for k = 1:ind
    C = strpt + k*res;
    F = @(x)   sin(C*x);
    Fp= @(x) C*cos(C*x);
    
    u1 = zeros(l1,nx); g1 = u1;
    u2 = zeros(l2,nx); g2 = u2;
    u2b=-1.6 * ones(1,nx-2);
    u2bold = u2b;
    for iter = 1:(50+64)
        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1 = Fp(u1(:));
            J1(ind1x) = 0; J1(ind1y) = 0;
            J1 = Lap1 - spdiags(J1,0,l1*nx,l1*nx);
            F1 = Lap1*u1(:) - F(u1(:));
            F1(ind1y) = 0;
            BCx1 = J1BC * u1(:,2:end-1) - [ BCL ; u2b ];
            F1(ind1x) = BCx1(:);
            u1(:) = u1(:) - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=(h+2*eps),2:end-1);
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = Fp(u2(:));
            J2(ind2x) = 0; J2(ind2y) = 0;
            J2 = Lap2 - spdiags(J2,0,l2*nx,l2*nx);
            F2 = Lap2*u2(:) - F(u2(:));
            F2(ind2y) = 0;
            BCx2 = J2BC * u2(:,2:end-1) - [ u1a ; BCR ];
            F2(ind2x) = BCx2(:);
            u2(:) = u2(:) - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=(h+2*eps),2:end-1);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= Fp(u1(:)); dF1(ind1x) = 0; dF1(ind1y) = 0;
        dF1= Lap1 - spdiags(dF1,0,l1*nx,l1*nx);
        g1(:) = dF1 \ kron([0 ; ones(nx-2,1) ; 0],[ zeros(l1-1,1); 1]);
        g1a= BCa * g1(abs(x1-a)<=(h+2*eps),2:end-1);

        % Step 4: solve g in second domain
        dF2= Fp(u2(:)); dF2(ind2x) = 0; dF2(ind2y) = 0;
        dF2= Lap2 - spdiags(dF2,0,l2*nx,l2*nx);
        g2(:) = dF2 \ kron([0 ; ones(nx-2,1) ; 0],[ 1; zeros(l2-1,1)]);
        g2b= BCb * g2(abs(x2-b)<=(h+2*eps),2:end-1); g2b = g2b .* g1a;

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)./(g2b - 1);
        u2bold= u2b;
        
        if iter > 50
            gam(k,iter-50,:) = u2b;
        end
        
    end
end

%%
figure(1)
aa = strpt + (1:ind)*res;
% subplot(2,1,1)
plot(aa,gam(:,:,round(nx/2)-1),'k.',aa,gam_neg(:,:,round(nx/2)-1),'r.')
% axis([strpt,endpt,1.4,1.5])
xlabel('a')
ylabel('\gamma')
set(gca,'fontsize',26,'linewidth',2)
% 
% subplot(2,1,2)
% plot(aa,gam,'k.',aa,gam_neg,'r.')
% % axis([strpt,endpt,-1.75,-1.35])
% axis([strpt,endpt,-1.5,-1.4])
% xlabel('a')
% ylabel('\gamma')
% set(gca,'fontsize',26,'linewidth',2)