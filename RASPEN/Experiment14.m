%% Experiment 14 - plotting the bifurcation diagram for u''-sin(au)=0
% Having found a period doubling example, we plot out its bifurcation
% diagram.

%% Bifurcation diagram

% Grid
nx = 1001;
a  = -0.1;
b  = 0.1;
x  = linspace(-1,1,nx)';
dx = 1./abs(x(2) - x(1));
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);
l1 = length(x1); l2 = length(x2);

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

% strpt = 3.52;
% endpt = 3.82;
% res   = 0.0001;
strpt = 2.1;
endpt = 2.4;
res   = 0.001;

ind = round((endpt - strpt) / res);
gam = zeros(ind,64);
nonlinsolves = 10;

for k = 1:ind
    C = strpt + k*res;
    F = @(x)   sin(C*x);
    Fp= @(x) C*cos(C*x);
    
    u1 = 0*ones(size(x1)); u1(1)  = 0;
    u2 = 0*ones(size(x2)); u2(end)= 0;
%     u2b= 1.6;
    u2b= 2;
    u2bold = u2b;
    for iter = 1:(50+64)
        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - F(u1(2:end-1));
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        u1a = u1(x1==a);
        
        % Step 2: solve u in second domain
        u2(1) = u1a;
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2(2:end-1,:)*u2 - F(u2(2:end-1));
            u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
        end
        u2b = u2(x2==b);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(Fp(u1),0,l1,l1);
        g1 = dF1(2:end-1,2:end-1) \ ( -dF1(2:end-1,end) );
        g1 = [0 ; g1 ; 1 ];
        g1a= g1(x1==a);

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(Fp(u2),0,l2,l2);
        g2 = dF2(2:end-1,2:end-1) \ ( -g1a*dF2(2:end-1,1) );
        g2 = [ g1a ; g2 ; 0];
        g2b= g2(x2==b);

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)/(g2b - 1);
        u2bold= u2b;
        
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
    u2b= -2;
    u2bold = u2b;
    for iter = 1:(50+64)
        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(Fp(u1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - F(u1(2:end-1));
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        u1a = u1(x1==a);
        
        % Step 2: solve u in second domain
        u2(1) = u1a;
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(Fp(u2),0,l2,l2);
            F2 = D2(2:end-1,:)*u2 - F(u2(2:end-1));
            u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
        end
        u2b = u2(x2==b);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(Fp(u1),0,l1,l1);
        g1 = dF1(2:end-1,2:end-1) \ ( -dF1(2:end-1,end) );
        g1 = [0 ; g1 ; 1 ];
        g1a= g1(x1==a);

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(Fp(u2),0,l2,l2);
        g2 = dF2(2:end-1,2:end-1) \ ( -g1a*dF2(2:end-1,1) );
        g2 = [ g1a ; g2 ; 0];
        g2b= g2(x2==b);

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
axis([strpt,endpt,1.5,2.5])
xlabel('a')
ylabel('\gamma')
set(gca,'fontsize',26,'linewidth',2)
subplot(2,1,2)
plot(aa,gam,'k.',aa,gam_neg,'r.')
% axis([strpt,endpt,-1.75,-1.35])
axis([strpt,endpt,-2.5,-1.5])
xlabel('a')
ylabel('\gamma')
set(gca,'fontsize',26,'linewidth',2)