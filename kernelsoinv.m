function [err,Puus,Puvs,Pvus,Pvvs] = kernelsoinv(lam,mu,c1,c2,p,N)
%KERNELSO Solves observer kernels for 2x2 system with constant parameters
%using power series approximation of order N. Returns the backstepping 
%output injection gains, i.e., the kernels along z=1 at points x=xu or x=xv.

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
Puu = 0; Puv = 0;
Pvu = 0; Pvv = 0;
for tm=0:N
  for tn=0:tm
    rind = rind + 1;
    Puu = Puu + KC(rind)*x^tn*z^(tm-tn);
    Puv = Puv + KC(NK+rind)*x^tn*z^(tm-tn);
    Pvu = Pvu + KC(2*NK+rind)*x^tn*z^(tm-tn);
    Pvv = Pvv + KC(3*NK+rind)*x^tn*z^(tm-tn);
  end
end

% construct equations and boundary conditions with power series
KE1 = L*diff(Puu,x) + diff(Puu,z)*subs(L,x,z) ... 
   + Puu*subs(diff(L,x),x,z) + Puv*subs(TH,x,z);
KE2 = L*diff(Puv,x) - diff(Puv,z)*subs(M,x,z) ... 
   - Puv*subs(diff(M,x),x,z) + Puu*subs(W,x,z);
KE3 = M*diff(Pvu,x) - diff(Pvu,z)*subs(L,x,z) ...
  - Pvu*subs(diff(L,x),x,z) - Pvv*subs(TH,x,z);
KE4 = M*diff(Pvv,x) + diff(Pvv,z)*subs(M,x,z) ...
  + Pvv*subs(diff(M,x),x,z) - Pvu*subs(W,x,z);
BC1 = subs(Puu,z,0) - Q*subs(Pvu,z,0);
BC2 = (L+M)*subs(Puv,z,x) - W;
BC3 = (L+M)*subs(Pvu,z,x) + TH;
BC4 = subs(Puv,z,0) - Q*subs(Pvv,z,0);
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
Puus = subs(Puu,KC(1:NK),Kc(1:NK)');
Puvs = subs(Puv,KC(NK+1:2*NK),Kc(NK+1:2*NK)');
Pvus = subs(Pvu,KC(2*NK+1:3*NK),Kc(2*NK+1:3*NK)');
Pvvs = subs(Pvv,KC(3*NK+1:4*NK),Kc(3*NK+1:4*NK)');
end