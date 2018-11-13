function [u2b, g, u1, u2] = G(gamma,u1,u2)

% Problem parameters
e = 0.1;

% Grid
nx = 101;
a  = -0.2;
b  = 0.2;
x  = linspace(-1,1,nx)';
x1 = x(x<=b); b = x1(end);
x2 = x(x>=a); a = x2(1);
l1 = length(x1); l2 = length(x2);

% Diff. mat.
I1 = eye(length(x1));
D1 = ones(length(x1)-2,1);
D1 = (nx-1)^2 * [ D1, -2*D1, D1 ];
D1 = spdiags(D1,[0,1,2],l1-2,l1);
D1 = [I1(1,:) ; e*D1 ; I1(end,:) ];
I2 = eye(length(x2));
D2 = ones(length(x2)-2,1);
D2 = (nx-1)^2 * [ D2, -2*D2, D2 ];
D2 = spdiags(D2,[0,1,2],l2-2,l2);
D2 = [I2(1,:) ; e*D2 ; I2(end,:) ];

d1 = ones(l1-2,1);
d1 = (nx-1)/2 * [ d1, zeros(l1-2,1), -d1 ];
d1 = spdiags(d1,[0,1,2],l1-2,l1);
d1 = [ zeros(1,l1) ; d1 ; zeros(1,l1) ];
d2 = ones(l2-2,1);
d2 = (nx-1)/2 * [ d2, zeros(l2-2,1), -d2 ];
d2 = spdiags(d2,[0,1,2],l2-2,l2);
d2 = [ zeros(1,l2) ; d2 ; zeros(1,l2) ];

% Jacobian
J1 = ones(l1-2,1);
J1 = [J1, -2*J1, J1];
J1 = e * J1 * (nx-1)^2;
J1 = spdiags(J1,[-1,0,1],l1-2,l1-2);
J2 = ones(l2-2,1);
J2 = [J2, -2*J2, J2];
J2 = e * J2 * (nx-1)^2;
J2 = spdiags(J2,[-1,0,1],l2-2,l2-2);

% Initialization
nonlinsolves = 1;

% Standard additive Schwarz
u2b = gamma;

if nargin==2
    u1 = u1(x<=b); u1(1) = 1;
    u2 = u1(x>=a); u2(end) = -1;
end

    % Step 1: solve u in first domain
    u1(end) = u2b;
    for j = 1:nonlinsolves
        f1 = -D1 * u1 + d1 * (u1.^2 / 2);
        JJ = spdiags([u1(2:end-1) , -u1(3:end)+u1(1:end-2) , -u1(2:end-1)],[-1,0,1],l1-2,l1);
        JJ = JJ(:,2:end-1);
        JJ = J1 + JJ*(nx-1)/2;
%             u1(2:end-1) = u1(2:end-1) + (JJ \ f1(2:end-1));
        JJ = JJ \ eye(l1-2); S1 = JJ(x1(2:end-1)==a,end);
        u1(2:end-1) = u1(2:end-1) + JJ * f1(2:end-1);
    end
    u1a = u1(x1==a);

    % Step 2: solve u in second domain
    u2(1) = u1a;
    for j = 1:nonlinsolves
        f2 = -D2 * u2 + d2 * (u2.^2/2);
        JJ = spdiags([u2(2:end-1) , -u2(3:end)+u2(1:end-2) , -u2(2:end-1)],[-1,0,1],l2-2,l2);
        JJ = JJ(:,2:end-1);
        JJ = J2 + JJ*(nx-1)/2;
%             u2(2:end-1) = u2(2:end-1) + (JJ \ f2(2:end-1));
        JJ = JJ \ eye(l2-2); S2 = JJ(x2(2:end-1)==b,1);
        u2(2:end-1) = u2(2:end-1) + JJ * f2(2:end-1);
    end
    u2b = u2(x2==b);
    
g = -( -e*(nx-1)^2 + u1(end-1)*(nx-1)/2 ) * ( e*(nx-1)^2 + u2(2)*(nx-1)/2 ) * S1 * S2;