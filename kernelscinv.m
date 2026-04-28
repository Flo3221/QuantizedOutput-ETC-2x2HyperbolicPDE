function [err,Kuus,Kuvs,Kvus,Kvvs] = kernelscinv(lam,mu,c1,c2,p,N)
%KERNELSK Solves control kernels for 2x2 system with constant parameters
%using power series approximation of order N. Returns the backstepping 
%control gains k1 and k2.

syms x z % symbolic variables for power series

% parameters
L = lam;
M = mu;
W = c1;
TH = c2;
Q = p;

% initalizations
NK = (N+1)*(N+2)/2; % number of terms in power series
KC = sym('k', [1 4*NK],'real'); % symbols for coefficients
rind = 0;
Kuu = 0; Kuv = 0;
Kvu = 0; Kvv = 0;
for tm=0:N
  for tn=0:tm
    rind = rind + 1;
    Kuu = Kuu + KC(rind)*x^tn*z^(tm-tn);
    Kuv = Kuv + KC(NK+rind)*x^tn*z^(tm-tn);
    Kvu = Kvu + KC(2*NK+rind)*x^tn*z^(tm-tn);
    Kvv = Kvv + KC(3*NK+rind)*x^tn*z^(tm-tn);
  end
end

% construct equations and boundary conditions with power series
KE1 = L*diff(Kuu,x) + diff(Kuu,z)*subs(L,x,z) ... 
   + Kuu*subs(diff(L,x),x,z) - Kvu*W;
KE2 = L*diff(Kuv,x) - diff(Kuv,z)*subs(M,x,z) ... 
   - Kuv*subs(diff(M,x),x,z) - Kvv*W;
KE3 = M*diff(Kvu,x) - diff(Kvu,z)*subs(L,x,z) ...
  - Kvu*subs(diff(L,x),x,z) + Kuu*TH;
KE4 = M*diff(Kvv,x) + diff(Kvv,z)*subs(M,x,z) ...
  + Kvv*subs(diff(M,x),x,z) + Kuv*TH;
BC1 = Q*L*subs(Kuu,z,0) - M*subs(Kuv,z,0);
BC2 = (L+M)*subs(Kuv,z,x) - W;
BC3 = (L+M)*subs(Kvu,z,x) + TH;
BC4 = M*subs(Kvv,z,0) - Q*L*subs(Kvu,z,0);
% extract coefficients for different powers of spatial variables
KE1C = coeffs(KE1,[x z]);
KE2C = coeffs(KE2,[x z]);
KE3C = coeffs(KE3,[x z]);
KE4C = coeffs(KE4,[x z]);
BC1C = coeffs(BC1,x);
BC2C = coeffs(BC2,x);
BC3C = coeffs(BC3,x);
BC4C = coeffs(BC4,x);
% construct set linear equations for the coefficients based on the above%
neqs = [numel(KE1C), numel(KE2C), numel(KE3C), numel(KE4C), ...
  numel(BC1C), numel(BC2C), numel(BC3C), numel(BC4C)];
A = zeros(sum(neqs), 4*NK); % initialize A...
B = zeros(sum(neqs),1); % ...and b
EQS = {KE1C,KE2C,KE3C,KE4C,BC1C,BC2C,BC3C,BC4C}; % cell array for equations
% disp(['Parsing data: ',num2str(sum(neqs)),' equations for ',...
%   num2str(4*NK),' unknonws']);
rind = 0; % row index
for ii = 1:8
  DAT = EQS{ii};
  for m=1:neqs(ii)
    rind = rind + 1;
    [cc, kc] = coeffs(DAT(m),KC);
    for mk=1:numel(kc) % check which K's are present in the coefficients
      str = char(kc(mk)); 
      sl = numel(str);
      if sl > 1 % look up the index and insert to A
        A(rind,str2double(str(2:sl))) = cc(mk);
      else % if no index, it's a contant; insert to B
        B(rind) = -cc(mk);
      end
    end
  end
  % disp(num2str(ii))
end
% solve equations and compute residual
Kc = A\B;
err = norm(A*Kc-B);
% insert the obtained values into the power series and compute gains
Kuus = subs(Kuu,KC(1:NK),Kc(1:NK)');
Kuvs = subs(Kuv,KC(NK+1:2*NK),Kc(NK+1:2*NK)');
Kvus = subs(Kvu,KC(2*NK+1:3*NK),Kc(2*NK+1:3*NK)');
Kvvs = subs(Kvv,KC(3*NK+1:4*NK),Kc(3*NK+1:4*NK)');
end