function Inverse = IOMcc(P,U1,U2,v,N)
%IOMcc Inverse operator matrix for constant coefficients.
%   R = IOMcc(P,U1,U2,v,N) is the inverse integration matrix for the linear
%   operator with constant coefficients defined in P, with rows indexed by
%   v removed and boundary conditions prescribed by U1 (at x = 1) and U2
%   (at x = -1). The length of v defines the order of the linear operator.
%   The input N defines the number of Chebyshev points being used.
%
%   R = IOMcc(P,U1,U2,[a,b,c,...],N) replaces the rows indexed by a, b, c, ...
%   with the boundary conditions specified in U1 and U2. The rows are
%   replaced in order, ie. if v = [a b c] and U1 has 2 rows, rows a and b
%   will be replaced with boundary conditions at x = 1, while row c will be
%   replaced with boundary conditions at x = -1.

%---Homog. Solns & their derivatives---%
% Size of inputs P and v
% for this version of IOM P should provide the constant coefficients of the
% operator

[T,dT,cmat] = formCheb(N);      % form necessary pieces of Chebyshev spectral methods
N = N+1;                        % avoid an off-by-one error

m = length(P); m = m-1;         % check length of P and v are related
if m~=length(v)
    disp('Error: dimension mismatch between P and v')
end

%---Roots and their multiplicities---%
% need to find the roots with their multiplicities and order them so mk(1)
% is the largest multiplicity
Y = roots(P);
y = zeros(1,m);
mk= y;
ind  = 1;
y(1) = Y(1);
mk(1)= sum(Y==y(1));
Y = Y(Y~=y(1));
while isempty(Y)~=1
    ind    = ind + 1;
    y(ind) = Y(1);
    mk(ind)= sum(Y==Y(1));
    Y = Y(Y~=Y(1));
end
M = ind;
y = y(1:M);
mk=mk(1:M);

%---Pascal's Triangle and Polynomials---%
PT = zeros(mk(1),mk(1)); PT(:,1) = ones(mk(1),1);
xx = ones(N,mk(1));
for i = 2:mk(1)
    PT(i,2:mk(1)) = PT(i-1,2:mk(1))+PT(i-1,1:mk(1)-1);
    xx(:,i) = xx(:,i-1) .* T(:,2) / (i-1);
end

%---Omega matrices---%
Om = zeros(m,m); ind = 0;
for i = 1:M
    indk = 1:mk(i);
    Om(:,ind+indk) = toeplitz(y(i).^(0:m-1),eye(1,mk(i))) .* PT(:,indk);
    ind = ind + mk(i);
end
z = Om \ [ zeros(m-1,1); 1];

%---Fundamental matrix, Wronskian and Beta coefficients---%
W = zeros(m,N);
Im= toeplitz(eye(1,mk(1)),(-1).^(0:mk(1)-1));
for i = 1:N
    F = toeplitz(eye(mk(1),1),xx(i,:));
    ind = 0;
    for k = 1:M
        indk = 1:mk(k);
        W(ind+indk,i) = exp(-y(k)*T(i,2)) * (Im(indk,indk).*F(indk,indk)) * z(ind+indk');
        ind = ind + mk(k);
    end
    if i==1
        F1 = F;
    elseif i==N
        F2 = F;
    end
end
Beta = W(:,v) \ W;

%---Fundamental solution set---%
% Arbitrary
P = ones(N,m);
ind = 0;
for i = 1:M
    indk =     1:mk(i);
    P(:,ind + indk) = exp(T(:,2)*y(i)) .* xx(:,indk);
    ind  = ind + mk(i);
end
% Particular
P = P * W(:,v);
% Derivatives at x=1 and -1
Pd = zeros(m,m,2);
ind = 0;
for k = 1:M
    indk =     1:mk(k);
    Pd(ind+indk,:,1  ) = exp( y(k)) * F1(indk,indk) * W(ind+indk,v);
    Pd(ind+indk,:,end) = exp(-y(k)) * F2(indk,indk) * W(ind+indk,v);
    ind  = ind + mk(k);
end

%---Interpolants---%
% G coefficient functions, homog. solns
B = zeros(N,N,m);
Inverse = zeros(N,N,m);
for k = 1:m
    B(:,:,k) = T' - T(v(k),:)'*T(:,N)'/T(v(k),N);     % (eq. 3.13)
    B(:,:,k) = cmat.*B(:,:,k);
    B(:,:,k) =   dT *B(:,:,k);
    Inverse(:,:,k) = ( P(:,k)*Beta(k,:) ).*B(:,:,k);  % (eq. 3.18)
end
Inverse = sum(Inverse,3);

%---Boundary conditions---% (eq. 3.19)
B  = permute( B,[3,2,1] );
if isempty(U2)
    Lhs = U1*Om*Pd(:,:,1);
    Rhs = Lhs*( Beta.*B(:,:,1) );
elseif isempty(U1)
    Lhs = U2*Om*Pd(:,:,end);
    Rhs = Lhs*( Beta.*B(:,:,end) );
else
    BCplus  = U1*Om*Pd(:,:,1); 
    BCminus = U2*Om*Pd(:,:,end);
    Lhs = [ BCplus ; BCminus ];
    Rhs = [ BCplus*( Beta.*B(:,:,1) ) ; BCminus*( Beta.*B(:,:,end) ) ];
end

Inverse = Inverse - P*( Lhs\Rhs );
Inverse(:,v) = P*( Lhs\eye(m) );

end

function [T,dT,cmat,D] = formCheb(N)
% formCheb - forms the Chebyshev polynomials and integrals as well as the
%       differentiation matrix and the Chebyshev points

% Differentiation matrix (eq. 1.2)
if nargout==4
    D = chebdifmat(N,1,1);
end

% Polynomials (eq. 1.5)
T = cos( (0:N)'*(0:N+2)*pi/N );

% Integrals (eq. 2.11)
dT = zeros(N+1,N+1);
dT(:,1)=T(:,2);
dT(:,2)=T(:,3)/4;
for k=3:N+1
    dT(:,k)=(1/2)*((T(:,k+1)/k) - (T(:,k-1)/(k-2)));
end

T = T(1:N+1,1:N+1);

% Weights (eq. 1.3 with modification)
c = ones(1,N+1);
c(1) = 2*c(1); c(end) = 2*c(end);
cmat = c'*c; cmat = (2/N)./cmat;

end