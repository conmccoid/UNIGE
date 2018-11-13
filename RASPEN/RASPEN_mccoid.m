%% RASPEN pseudocode

% define the problem: F, J, Ri, Pi, PPi

% u0 = initial guess;

% while abs(un - un-1) < tol

%   for i = run over all subdomains
%       solve RiF( v ) = 0 using Newton (to presumably some tolerance?):
%       v0 = 0;
%       while abs(vm - vm-1) < tol
%           Ji = inv( Ri J(vm) Pi ); [solve using GMRES]
%           w = -Ji * RiF(vm-1);
%           vm = vm-1 + w;
%       end
%       PiCi(un) = vm - un;
%       JJ = JJ + ( Pi Ji Ri J(vm) ); [probably a better way to do this]
%   end

%   FF(un) = sum( PPi Ri PiCi(un) );
%   use Newton's method again on FF() [need the Jacobian for FF]:
%   w = ( JJ ) \ FF(un); [do this by GMRES?]
%   un+1 = un + w;

% end