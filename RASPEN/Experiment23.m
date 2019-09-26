%% Experiment 23 - measuring basin of attraction as function of matrix properties
% Repeat of exp15 and comparing basin of attraction size to certain matrix
% properties.

% It is clear from experiments that the appropriate measure of cycling
% behaviour is the lowest magnitude of norm(u2b) within a given number of
% iterations for a single initial guess. There is a steep drop-off in this
% measure corresponding to loss of cycling behaviour (expected). Also of
% note is the largest magnitude of norm(u2b) for the same initial guess
% within a given number of iterations after some amount of stabilization.
% In some sense this is the measure of how many iterations is required to
% achieve the correct result (0) which cannot be used in this instance as
% the answer is infinite in the presence of cycling.

% The appropriate measure of the problem that provides insight into why
% this occurs remains elusive.

%% Bifurcation diagram

% Grid
nx = 21; ny = 21;
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
II = speye(ny); IIint = II; IIint(1) = 0; IIint(end) = 0;


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
% DDBC = sparse([1,1,1,2,2,2],[1,2,3,nx-2,nx-1,nx],...
%     0.5 * dx * [-3, 4, -1, 1, -4, 3],2,nx); % Neumann BCs in y
DDBC = [II(1,:) ; II(end,:)]; % Dirichlet BCs in y
D1 = [J1BC(1,:) ; D1(2:end-1,:) ; J1BC(2,:)];
D2 = [J2BC(1,:) ; D2(2:end-1,:) ; J2BC(2,:)];

ind1x= kron(IIint,I1-I1int)* ones(ny*l1,1) == 1;
ind1y= kron(II-IIint,   I1)* ones(ny*l1,1) == 1;
ind2x= kron(IIint,I2-I2int)* ones(ny*l2,1) == 1;
ind2y= kron(II-IIint,   I2)* ones(ny*l2,1) == 1;

% Plotting bifurcation region
strpt = 1.6;
endpt = 1.61;
res   = 0.01;

fx  = ones(1,ny-2);
ind = round((endpt - strpt) / res);
aa  = strpt + (1:ind)*res;
nonlinsolves = 10;
stab= 50;

N = 100;
basin_a = zeros(1,N);
basin_b = basin_a;
basin_c = basin_a;
basin_d = basin_a;
Cond1   = basin_a;
Cond2   = basin_a;
% E = linspace(0.0004,0.0007,N);
E = logspace(-5,-2,N);

C = 3.65;
F = @(x)   sin(C*x);
Fp= @(x) C*cos(C*x);

itermax = 64;
gam = zeros(ind,itermax,N);

for j = 1:N
    k = 0; konst = strpt;
    epsilon = E(j); %1e-3 gives period doubling
    DD = epsilon * ones(ny,1); % perturbed y derivatives
    DD = dy^2 * [ DD, -2*DD, DD ];
    DD = spdiags(DD,[-1,0,1],ny,ny);
    DD = [DDBC(1,:) ; DD(2:end-1,:) ; DDBC(2,:)];

    % Laplacian
    Lap1 = kron(IIint,D1) + kron(DD,I1int) + kron((II-IIint)*DD,I1-I1int);
    Lap2 = kron(IIint,D2) + kron(DD,I2int) + kron((II-IIint)*DD,I2-I2int);
    Cond1(j)= (norm(Lap1-diag(diag(Lap1)),1))./norm(Lap1,1); % change this to reflect a different aspect of the problem
    
    while basin_b(j)==0 && k<ind
        k = k+1;
        konst = konst+res;
        iter  = 0;

        u1 = zeros(l1,ny);
        u2 = zeros(l2,ny);
        u2b= konst*fx; %multiply by appropriate number depending on k
        u2bold = u2b;
                
        error  = norm(u2b);
        while iter<(stab+itermax) && error>eps
            iter = iter+1;
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
            u2bold= u2b;
            error = norm(u2b,1);

            if iter > stab
                gam(k,iter-stab,j) = error;
            end

        end
        if k>1
            if gam(k-1,end,j)<eps && gam(k,end,j)>eps
                basin_a(j) = konst-res;
                k_keep  = k;
            elseif gam(k-1,end,j)>eps && gam(k,end,j)<eps
                gam_keep = gam(k_keep:k-1,:,j);
                basin_b(j) = konst;
                basin_c(j) = min(gam_keep(:));
                basin_d(j) = max(gam_keep(:));
                Cond2(j)= condest(dF1);
                
%                 figure(1)
%                 plot(aa,gam(:,:,j)','k.')
%                 axis([1,3,0,10])
%                 xlabel('a')
%                 ylabel('$$\Vert \gamma \Vert$$','interpreter','latex')
%                 set(gca,'fontsize',26,'linewidth',2)
%                 pause(eps)
            end
        end
    end
end

gam = permute(gam,[2,3,1]);

%%
figure(1)
plot(E,basin_d)

figure(2)
semilogx(E,Cond1-0.5)

%%
figure(3)
plot(E,gam,'.')

figure(4)
plot(E,min(gam,[],1))