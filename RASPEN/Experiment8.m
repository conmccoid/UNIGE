%% Experiment 8 - seeking to break the method
% Having shown examples where preconditioning by Newton severly improves
% the convergence, we'd now like to find an example where the precondition
% fails. We do so by first finding a function G(y) for which Newton does
% poorly, or has a very small basin of attraction. We then find a PDE that
% would produce such a G(y) as part of the preconditioning.

%% Finding a function G(y)

% Grid
y = linspace(-1,1,1001);
a = 0.25;
b = 0.5;
G = @(y,c) -a*y - (b/c) * sin(c*y);
Gp= @(y,c) -a   - b*cos(c*y);

c = 10;
tol = 1e-6;
itermax = 100;
save = zeros(size(y));

for k = 1:length(y)
    y0 = y(k);
    error= 1;
    iter = 1;
    while error > tol && iter < itermax
%         y1 = y0 - (G(y0,c) - y0)./(Gp(y0,c) - 1);
%         y1 = y0 - G(y0,c)./Gp(y0,c);
        y1 = G(y0,c)+y0;
        error = abs(y1 - y0);
        y0 = y1;
        iter = iter+1;
        
%         plot(y,G(y,c),'k',y0,G(y0,c),'ro')
%         axis([-1,1,-1,1])
%         pause(0.1)
    end
    
    if error < tol && iter < itermax
        save(k) = iter;
    else
        save(k) = NaN;
    end

end

plot(y,save,'b*--')

%% Implementing test

% Problem parameters
C = 3.6;

% Grid
nx = 1001;
a  = -0.2;
b  = 0.2;
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
D1 = dx^2 * [ D1, -2*D1, D1 ];% - dx/2 * [ -D1, zeros(l1,1), D1 ].*X1;
D1 = spdiags(D1,[-1,0,1],l1,l1);% + spdiags(ones(l1,1),0,l1,l1);
I2 = eye(l2);
D2 = ones(l2,1);
D2 = dx^2 * [ D2, -2*D2, D2 ];% - dx/2 * [ -D2, zeros(l2,1), D2 ].*X2;
D2 = spdiags(D2,[-1,0,1],l2,l2);% + spdiags(ones(l2,1),0,l2,l2);

% Initialization
k = 0;
tol = 1e-8;
itermax = 100;
G = zeros(1,101);
Gp= G;
itersaveNewton = G;
itersaveReg = G;
testu2 = linspace(-2,2,101);
% testu2=0.03;
nonlinsolves = 50;
for u2b0 = testu2
% u2b0 = 0;
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = 0*ones(size(x1)); u1(1)  = 0;
    u2 = 0*ones(size(x2)); u2(end)= 0;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(C*cos(C*u1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - sin(C*u1(2:end-1));
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        u1a = u1(x1==a);
        
        % Step 2: solve u in second domain
        u2(1) = u1a;
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(C*cos(C*u2),0,l2,l2);
            F2 = D2(2:end-1,:)*u2 - sin(C*u2(2:end-1));
            u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
        end
        u2b = u2(x2==b);

        % Preconditioning with Newton
        % Step 3: solve g in first domain
        dF1= D1 - spdiags(C*cos(C*u1),0,l1,l1);
        g1 = dF1(2:end-1,2:end-1) \ ( -dF1(2:end-1,end) );
        g1 = [0 ; g1 ; 1 ];
        g1a= g1(x1==a);

        % Step 4: solve g in second domain
        dF2= D2 - spdiags(C*cos(C*u2),0,l2,l2);
        g2 = dF2(2:end-1,2:end-1) \ ( -g1a*dF2(2:end-1,1) );
        g2 = [ g1a ; g2 ; 0];
        g2b= g2(x2==b);
        
%         if iter==2
%             Gp(k) = g2b;
%         end

        % Step 5: update u2b
        u2b = u2bold - (u2b - u2bold)/(g2b - 1);
        
        if iter==2
            Gp(k) = u2b;
        end

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
%         if iter==2
% %         if error < tol
%             subplot(1,2,1)
%             plot(x1,g1,x2,g2)
%             axis([-1,1,-3,3])
%             subplot(1,2,2)
%             plot(x1,u1,x2,u2)
% %             axis([-1,1,-3,3])
%             pause(0.1)
%         end
        
    end
    
    if iter < itermax && error < tol
        itersaveNewton(k) = iter;
    else
        itersaveNewton(k) = NaN;
    end
    
%     plot(x1,u1,x2,u2)
%     axis([-1,1,-10,10])
%     title(num2str(G(k)))
%     pause(0.01)
        
end

k=0;
for u2b0 = testu2
% u2b0 = 0;
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = 0*ones(size(x1)); u1(1)  = 0;
    u2 = 0*ones(size(x2)); u2(end)= 0;
    while error > tol && iter < itermax

        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
            J1 = D1 - spdiags(C*cos(C*u1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - sin(C*u1(2:end-1));
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        u1a = u1(x1==a);
        
        % Step 2: solve u in second domain
        u2(1) = u1a;
        for i = 1:nonlinsolves
            J2 = D2 - spdiags(C*cos(C*u2),0,l2,l2);
            F2 = D2(2:end-1,:)*u2 - sin(C*u2(2:end-1));
            u2(2:end-1) = u2(2:end-1) - J2(2:end-1,2:end-1) \ F2;
        end
        u2b = u2(x2==b);
        
        if iter==2
            G(k) = u2b;
        end

        error = abs(u2b - u2bold);
        u2bold= u2b;
        iter  = iter+1;
        
    end
    
    if iter < itermax && error < tol
        itersaveReg(k) = iter;
    else
        itersaveReg(k) = NaN;
    end
        
end

figure(1)
plot(testu2,itersaveReg,'b',testu2,itersaveNewton,'r')
xlabel('\gamma')
ylabel('Number of iterations to convergence')
title('Newton precond. on transmission condition')

figure(2)
plot(testu2,G,'r',testu2,Gp,'k',testu2,testu2,'g','linewidth',2)
xlabel('\gamma')
ylabel('G(\gamma)')
legend('FP','NR','\gamma')
axis([-2,2,-2,2])
set(gca,'fontsize',26,'linewidth',2)

%% Mapping out the function G(x,y)

% Problem parameters
e = 1;

% Grid
nx = 1001;
a  = -0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
dx = 1./abs(x(2) - x(1));
x1 = x(x<=b); b = x1(end);
l1 = length(x1);

% Diff. mat.
X1 = [ [x1(2:end) ; 0], x1, [0 ; x1(1:end-1)] ];

I1 = eye(l1);
D1 = ones(l1,1);
D1 = e * dx^2 * [ D1, -2*D1, D1 ];% - dx/2 * [ -D1, zeros(l1,1), D1 ].*X1;
D1 = spdiags(D1,[-1,0,1],l1,l1);% + spdiags(ones(l1,1),0,l1,l1);
d1 = dx/2 * [-ones(l1,1), zeros(l1,1), ones(l1,1)];
d1 = spdiags(d1,[-1,0,1],l1,l1);

% Initialization
k = 0;
k0= 1.5;
tol = 1e-8;
itermax = 100;
G = zeros(l1,101);
Gp= G;
itersaveNewton = G;
testu2 = linspace(-5,5,101);
nonlinsolves = 50;
for u2b0 = testu2
% u2b0 = 0;
    u2b    = u2b0;
    u2bold = u2b0;
    k      = k+1;
    error  = 1;
    iter   = 1;
    u1 = 1*ones(size(x1)); u1(1)  = 1;

        % Step 1: solve u in first domain
        u1(end) = u2b;
        for i = 1:nonlinsolves
%             J1 = D1 - spdiags(u1,0,l1,l1) * d1 - spdiags(d1*u1,0,l1,l1);
%             F1 = D1(2:end-1,:)*u1 - u1(2:end-1).* (d1(2:end-1,:)*u1);
            J1 = D1 - spdiags(u1.^(k0-1),0,l1,l1);
            F1 = D1(2:end-1,:)*u1 - u1(2:end-1).^k0/k0;
            u1(2:end-1) = u1(2:end-1) - J1(2:end-1,2:end-1) \ F1;
        end
        
        G(:,k) = u1;
        
        % Preconditioning with Newton
        % Step 3: solve g in first domain
%         dF1= D1 - spdiags(u1,0,l1,l1) * d1 - spdiags(d1*u1,0,l1,l1);
        dF1= D1 - spdiags(u1.^(k0-1),0,l1,l1);
        g1 = dF1(2:end-1,2:end-1) \ ( -dF1(2:end-1,end) );
        g1 = [0 ; g1 ; 1];
        g1a= g1(x1==a);
        
        Gp(:,k) = g1;
        
end

figure(1)
contour(testu2,x1,G-testu2,[0 0])
ylabel('x')
xlabel('\gamma')

figure(2)
plot(testu2,G(1:100:end,:))
axis([-7,7,-7,7])
xlabel('\gamma')
ylabel('G(\gamma)')

figure(3)
surf(testu2,x1,G)

figure(4)
contour(testu2,x1,Gp,[1,1])