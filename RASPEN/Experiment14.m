%% Experiment 14 - plotting the bifurcation diagram for u''-sin(au)=0
% Having found a period doubling example, we plot out its bifurcation
% diagram. We may change the overlap and transmission condition to see how
% this bifurcation changes.

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
ind1 = 1:l1; ind1 = ind1(abs(x1-a)<=(h+2*eps));
ind2 = 1:l2; ind2 = ind2(abs(x2-b)<=(h+2*eps));

% Diff. mat.
X1 = [ [x1(2:end) ; 0], x1, [0 ; x1(1:end-1)] ];
X2 = [ [x2(2:end) ; 0], x2, [0 ; x2(1:end-1)] ];

I1 = eye(l1);
D1 = ones(l1,1);
D1 = dx^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[-1,0,1],l1,l1);
I2 = eye(l2);
D2 = ones(l2,1);
D2 = dx^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[-1,0,1],l2,l2);

% BC parameters
a1 = 0; b1 = 1;
a2 = 0; b2 = 1;
c1 = 0; d1 = 1; BCb = sparse([1,1,1],ind2,[-c1*dx, d1, c1*dx],1,l1);
c2 = 0; d2 = 1; BCa = sparse([1,1,1],ind1,[-c2*dx, d2, c2*dx],1,l2);
BCL= 0;
BCR= 0;

% Jacobian BCs
J1BC = sparse([1,1,1,2,2,2],[1,2,3,l1-2,l1-1,l1],...
    [-3*a1/(2*h) + b1,4*a1/(2*h),-a1/(2*h),c1/(2*h),-4*c1/(2*h),3*c1/(2*h)+d1],2,l1);
J2BC = sparse([1,1,1,2,2,2],[1,2,3,l2-2,l2-1,l2],...
    [-3*c2/(2*h) + d2,4*c2/(2*h),-c2/(2*h),a2/(2*h),-4*a2/(2*h),3*a2/(2*h)+b2],2,l2);

% Plotting bifurcation region
strpt = 3.52;
endpt = 3.82;
res   = 0.01;
% strpt = 9.10;
% endpt = 9.155;
% res   = 0.001;

ind = round((endpt - strpt) / res);
gam = zeros(ind,64);
nonlinsolves = 10;

for k = 1:ind
    C = strpt + k*res;
    F = @(x)   sin(C*x);
    Fp= @(x) C*cos(C*x);
    
    u1 = 0*ones(size(x1));
    u2 = 0*ones(size(x2));
%     u2b= 1.6;
    u2b= 1.68;
    u2bold = u2b;
    for iter = 1:(50+64)
        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1*u1 - F(u1);
            J1 = [J1BC(1,:);J1(2:end-1,:);J1BC(2,:)];
            F1(1) = J1BC(1,:)*u1 - BCL; F1(end) = J1BC(2,:)*u1 - u2b;
            u1 = u1 - J1 \ F1;
        end
        u1a= BCa * u1;
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2*u2 - F(u2);
            J2 = [J2BC(1,:);J2(2:end-1,:);J2BC(2,:)];
            F2(1) = J2BC(1,:)*u2 - u1a; F2(end) = J2BC(2,:)*u2 - BCR;
            u2 = u2 - J2 \ F2;
        end
        u2b = BCb * u2;

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(Fp(u1),0,l1,l1);
        dF1= [J1BC(1,:);dF1(2:end-1,:);J1BC(2,:)];
        g1 = dF1 \ [ zeros(l1-1,1); 1];
        g1a= BCa * g1;

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(Fp(u2),0,l2,l2);
        dF2= [J2BC(1,:);dF2(2:end-1,:);J2BC(2,:)];
        g2 = dF2 \ [ g1a; zeros(l2-1,1)];
        g2b= BCb * g2;

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)/(g2b - 1);
        u2bold= u2b;
        
        plot(x1,u1)
        pause
        
        if iter > 50
            gam(k,iter-50) = u2b;
        end
        
    end
end

gam_neg = zeros(size(gam));

for k = 1:ind
    C = strpt + k*res;
    F = @(x)   sin(C*x);
    Fp= @(x) C*cos(C*x);
    
    u1 = 0*ones(size(x1)); u1(1)  = 0;
    u2 = 0*ones(size(x2)); u2(end)= 0;
%     u2b= -1.6;
    u2b= -1.45;
    u2bold = u2b;
    for iter = 1:(50+64)
        % Step 1: solve u in first domain
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1*u1 - F(u1);
            J1 = [J1BC(1,:);J1(2:end-1,:);J1BC(2,:)];
            F1(1) = J1BC(1,:)*u1 - BCL; F1(end) = J1BC(2,:)*u1 - u2b;
            u1 = u1 - J1 \ F1;
        end
        u1a= BCa * u1(abs(x1-a)<=(h+2*eps));
        
        % Step 2: solve u in second domain
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2*u2 - F(u2);
            J2 = [J2BC(1,:);J2(2:end-1,:);J2BC(2,:)];
            F2(1) = J2BC(1,:)*u2 - u1a; F2(end) = J2BC(2,:)*u2 - BCR;
            u2 = u2 - J2 \ F2;
        end
        u2b = BCb * u2(abs(x2-b)<=(h+2*eps));

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(Fp(u1),0,l1,l1);
        dF1= [J1BC(1,:);dF1(2:end-1,:);J1BC(2,:)];
        g1 = dF1 \ [ zeros(l1-1,1); 1];
        g1a= BCa * g1(abs(x1-a)<=(h+2*eps));

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(Fp(u2),0,l2,l2);
        dF2= [J2BC(1,:);dF2(2:end-1,:);J2BC(2,:)];
        g2 = dF2 \ [ g1a; zeros(l2-1,1)];
        g2b= BCb * g2(abs(x2-b)<=(h+2*eps));

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)/(g2b - 1);
        u2bold= u2b;
        
        if iter > 50
            gam_neg(k,iter-50) = u2b;
        end
        
    end
end

%%
figure(1)
aa = strpt + (1:ind)*res;
subplot(2,1,1)
plot(aa,gam,'k.',aa,gam_neg,'r.')
% axis([strpt,endpt,1.35,1.75])
axis([strpt,endpt,1.4,1.5])
xlabel('a')
ylabel('\gamma')
set(gca,'fontsize',26,'linewidth',2)

subplot(2,1,2)
plot(aa,gam,'k.',aa,gam_neg,'r.')
% axis([strpt,endpt,-1.75,-1.35])
axis([strpt,endpt,-1.5,-1.4])
xlabel('a')
ylabel('\gamma')
set(gca,'fontsize',26,'linewidth',2)