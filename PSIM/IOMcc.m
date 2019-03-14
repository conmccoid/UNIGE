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
y = sort(roots(P));                   % calculates the characteristics(?) of the operator (column vector)
Y = y(1);
M = length(y==Y);
y = y(y~=Y);
while isempty(y)~=1
    newy = y(1);
    Y = [ Y ; newy ];
    M = [ M ; length(y==newy) ];
    y = y(y~=newy);
end

Y = bsxfun(@power,y',(0:m-1)');

w = zeros(m,1);
ind = 1:m;
for j = 1:m
    w(j) = (-1)^(j+m) * det(Y(1:m-1,ind~=j));
end
w = w./det(Y);

P = exp(T(:,2)*y');             % this needs to change to account for roots with multiplicity
Pm= P * diag(w) * (1./P(v,:)');

%---Pascal's Triangle---%
PT = spalloc(m,m,m*(m+1)/2); PT(:,1) = ones(m,1);
for i = 2:m
    PT(i,2:m) = PT(i-1,2:m)+PT(i-1,1:m-1);
end

%---Omega matrices---%
Om = sparse(m,m); ind = 0;
for i = 1:M
    Omtemp = spdiags(0:-1:1-m,y(i).^(0:m-1),m,mk(i));
    Om(:,ind+(1:mk(i))) = Omtemp.*PT(:,1:mk(i));
    ind = ind + mk(i);
end

detOm = det(Om);
detOmmat = zeros(m,m);
xx = ones(N,m); ee = zeros(N,m);
ind = 0;
for i = 1:M
    for j = 1:mk(i)
        detOmmat(ind + (j:-1:1),ind + (1:j)) = det(Om(1:m-1,(1:m)~=(ind+j))); %block diagonal and symmetric
        if ind==0 || j>mk(1)
            xx(:,ind+j+1) = xx(:,ind+j)*T(:,2)/j;
        else
            xx(:,ind+j+1) = xx(:,j+1);
        end
    end
    ee(ind + 1:mk(i)) = exp(-T(:,2)*y(i));
    ind = ind + mk(i);
end

%---Wronskians---%
W = ee.*(xx * detOmmat)./ detOm;
Gam = W(v,:);

%---Fundamental solution set(s)---%
Phat = ee .* xx; %hopefully unnecessary in final version? although this probably doesn't cost much
P = Phat * Gam;

% %---Beta coefficients---%
% Coeffs = zeros(m,m);
% for j = 1:m
%     for k = 1:m
%         Coeffs(j,k) = (-1)^(j+k) * det( 1./P(v(ind~=k),ind(ind~=j))' );
%     end
% end
% W = Coeffs' * (1./P') ./ det(1./P(v,:)');

%---Interpolants---%
% G coefficient functions, homog. solns
B = zeros(N,N,m);
Inverse = zeros(N,N,m);
for k = 1:m
    B(:,:,k) = T' - T(v(k),:)'*T(:,N)'/T(v(k),N);     % (eq. 3.13)
    B(:,:,k) = cmat.*B(:,:,k);
    B(:,:,k) = dT*B(:,:,k);
    Inverse(:,:,k) = ( Pm(:,k)*W(k,:) ).*B(:,:,k);  % (eq. 3.18)
end
Inverse = sum(Inverse,3);

%---Boundary conditions---% (eq. 3.19)
% Pm = permute( Pm,[3,2,1] );
Pd = zeros(m,m,2);
Pd(1,:,1) = Pm(1,:);
Pd(1,:,end) = Pm(end,:);
for l=2:m
    Pd(l,:,:) = P([1,end],:) * diag(w .* y.^(l-1)) * (1./P(v,:)');
end
B  = permute( B,[3,2,1] );

if isempty(U2)
    Lhs = U1*Pd(:,:,1);
    Rhs = Lhs*( W.*B(:,:,1) );
elseif isempty(U1)
    Lhs = U2*Pd(:,:,end);
    Rhs = Lhs*( W.*B(:,:,end) );
else
    BCplus  = U1*Pd(:,:,1); 
    BCminus = U2*Pd(:,:,end);
    Lhs = [ BCplus ; BCminus ];
    Rhs = [ BCplus*( W.*B(:,:,1) ) ; BCminus*( W.*B(:,:,end) ) ];
end

Inverse = Inverse - Pm*( Lhs\Rhs );
Inverse(:,v) = Pm*( Lhs\eye(m) );

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