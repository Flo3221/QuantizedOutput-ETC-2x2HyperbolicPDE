% ================================================================
% Computes all theoretical bounds, verifies feasibility conditions,
% selects parameters (sigma, epsilon, delta_sg) for the small-gain
% condition, and saves constants for Simulink.
%% System parameters
a  = @(x) ones(size(x));
b  = @(x) 2*ones(size(x));
c1 = @(x) 2*ones(size(x));
c2 = @(x) 2*ones(size(x));
p  = 0.6;
r  = 0.9;

xg  = 128;
tmp = linspace(0,1,xg+1);
x0  = tmp(1:xg);  x1 = tmp(2:xg+1);

%% Kernels
N = 13;
[k1,k2,ec,Kuu,Kuv,Kvu,Kvv]      = kernelsc(a(1),b(1),c1(1),c2(1),p,r,x1,x0,N);
[eci,Laa,Lab,Lba,Lbb]            = kernelscinv(a(1),b(1),c1(1),c2(1),p,N);
[p1,p2,eo,Puu,Puv,Pvu,Pvv]      = kernelso(a(1),b(1),c1(1),c2(1),p,x1,x0,N);
[eoi,Raa,Rab,Rba,Rbb]            = kernelsoinv(a(1),b(1),c1(1),c2(1),p,N);

%% Symbolic bounds
syms x z
phi1 = @(x) x/a(1);  phi2 = @(x) x/b(1);                   % (1)

p1s   = -a(1)*subs(Puu,z,1);                                % (2)
p2s   = -a(1)*subs(Pvu,z,1);                                % (3)
P1bar = p1s - int(Kuu*subs(p1s,x,z)+Kuv*subs(p2s,x,z),z,0,x);
P2bar = p2s - int(Kvu*subs(p1s,x,z)+Kvv*subs(p2s,x,z),z,0,x);
[~,tmp] = fminbnd(matlabFunction(-P1bar),0,1);  P1 = -tmp;  % (5)
[~,tmp] = fminbnd(matlabFunction(-P2bar),0,1);  P2 = -tmp;

Na = subs(Lba,x,1) - r*subs(Laa,x,1);                      % (7)
Nb = subs(Lbb,x,1) - r*subs(Lab,x,1);

a0  = double(abs(subs(Nb,z,1)*b(1)));
C2  = double(int(abs(diff(Na*a(z),z)),z,0,1) + abs(subs(Na,z,1)*a(1) - subs(Nb,z,1)*b(1)*r));
C3  = double(int(abs(diff(Nb*b(z),z)),z,0,1) + abs(subs(Nb,z,0)*b(0) - subs(Na,z,0)*a(0)*p));
C1  = double(max([C2,C3]));
Cth = 1/(1-abs(p*r)) * max([1+abs(r), 1+abs(p), (1+abs(r))*P1*phi1(1)+(1+abs(p))*P2*phi2(1)+abs(r)*(1+abs(p))]);
Kq  = P1/a(1) + (1+abs(p))*(abs(r)+P2/b(1));
Cp  = double(abs(int(Na*subs(P1bar,x,z)+Nb*subs(P2bar,x,z),z,0,1)) + abs(subs(Nb,z,1)*b(1)*r));

Pi = int(abs(Puu)+abs(Puv)+abs(Pvu)+abs(Pvv),z,x,1);       % (10)
Ri = int(abs(Raa)+abs(Rab)+abs(Rba)+abs(Rbb),z,x,1);       % (11)
Ki = int(abs(Kuu)+abs(Kuv)+abs(Kvu)+abs(Kvv),z,0,x);       % (12)
Li = int(abs(Laa)+abs(Lab)+abs(Lba)+abs(Lbb),z,0,x);       % (13)
[~,tmp] = fminbnd(matlabFunction(-Pi),0,1);  Pbar = -tmp;
[~,tmp] = fminbnd(matlabFunction(-Ri),0,1);  Rbar = -tmp;
[~,tmp] = fminbnd(matlabFunction(-Ki),0,1);  Kbar = -tmp;
[~,tmp] = fminbnd(matlabFunction(-Li),0,1);  Lbar = -tmp;
M1 = max([1+Rbar, 1+Kbar]);                                  % (8)
M2 = min([1/(1+Pbar), 1/(1+Lbar)]);                          % (9)

%% Quantizer / STC design parameters
M     = 3;
Delta = M/10;
mu0   = 8;
Omega = 0.6;
T     = 2;
Cetm  = 0.14;
delta = -log(Omega)/T;   % positive decay rate  (Omega < 1 => delta > 0)

%% Condition (14): upper bound on Cetm
% Cetm < Cetmax = (1 - |p*r|) / (2*(1+|p|))
Cetmax = (1 - abs(p*r)) / (2*(1+abs(p)));
fprintf('\n--- Condition (14) ---\n');
fprintf('  Cetmax = %.6f\n', Cetmax);
fprintf('  Cetm   = %.6f  -->  satisfied: %d\n', Cetm, Cetm < Cetmax);
assert(Cetm < Cetmax, 'Cetm must be < Cetmax (condition 14).');

%% Condition (15): phi(0,0,0)
% phi(0,0,0) = term1 + (numer_frac / denom15) * factor2
% denom15 > 0 is guaranteed by condition (14)
amin    = a(1);
bmin    = b(1);
denom15 = 1 - abs(p*r) - 2*Cetm*(1+abs(p));
term1   = P1/amin + (2+abs(p))*(abs(r)+P2/bmin);
numer_frac = (1+abs(r))*P1/amin ...
           + (1+abs(p))*P2/bmin ...
           + (abs(p)+1)*(abs(r)+P2/bmin);
factor2 = 1 + P1/amin + (abs(p)+1)*(P2/bmin+abs(r));
phi000  = term1 + (numer_frac/denom15)*factor2;

fprintf('\n--- Condition (15) ---\n');
fprintf('  denom15    = %.6f  (> 0 by (14))\n', denom15);
fprintf('  phi(0,0,0) = %.6f\n', phi000);

%% Lambda feasibility window
% lambda must simultaneously satisfy:
%   (15): phi000 < 1+lambda  =>  lambda > phi000/target - 1
%   (16): Delta/M < M2/(1+lambda)*exp(-2*delta*T)/Omega^2
%         => lambda < M2*exp(-2*delta*T)/(Omega^2*(Delta/M)) - 1
%
% With delta = -log(Omega)/T:
%   exp(-2*delta*T) = exp(2*log(Omega)) = Omega^2
%   => lambda_max16 = M2/(Delta/M) - 1 = M2*M/Delta - 1
%
% To open the window we therefore need:  M2*M/Delta > phi000/target
% i.e.  Delta < M2*M*target/phi000

target       = 0.95;
lambda_max16 = M2*exp(-2*delta*T) / (Omega^2*(Delta/M)) - 1;
lambda_min15 = phi000/target - 1;

fprintf('\n--- Lambda feasibility window ---\n');
fprintf('  lambda_min needed  (from (15), target=%.2f) = %.6f\n', target, lambda_min15);
fprintf('  lambda_max allowed (from (16))              = %.6f\n', lambda_max16);

if lambda_min15 >= lambda_max16
    % Compute the Delta required to open the window
    Delta_needed = M2 * M * target / phi000;
    fprintf('\n  *** Window closed with current Delta=%.4f ***\n', Delta);
    fprintf('  To open it, need Delta < %.2e\n', Delta_needed);
    fprintf('  Setting Delta = %.2e  (factor 10 below threshold)\n\n', Delta_needed/10);
    Delta        = Delta_needed / 10;
    lambda_max16 = M2*exp(-2*delta*T) / (Omega^2*(Delta/M)) - 1;
    lambda_min15 = phi000/target - 1;
    fprintf('  Revised lambda_min = %.6f,  lambda_max = %.6f\n', ...
            lambda_min15, lambda_max16);
    assert(lambda_min15 < lambda_max16, ...
        ['Window still closed after Delta adjustment. ' ...
         'Reduce Cetm closer to 0 or reduce p*r.']);
end

lambda = (lambda_min15 + lambda_max16) / 2;   % midpoint for maximum margin
fprintf('  lambda chosen      = %.6f  (midpoint)\n', lambda);
fprintf('  phi000/(1+lambda)  = %.6f  (must be < 1)\n', phi000/(1+lambda));

%% Condition (16): Ccond
Ccond = M2 / (1+lambda) * exp(-2*delta*T) / Omega^2;
fprintf('\n--- Condition (16) ---\n');
fprintf('  Ccond   = %.6f\n', Ccond);
fprintf('  Delta/M = %.6f  -->  satisfied: %d\n', Delta/M, Delta/M < Ccond);

%% Self-triggering interval
Ak  = @(wk) Cth*C1*wk + (Cp+Cth*C1)*(M2*M/(Omega*exp(delta*T))+(1+Kq)*Delta)*mu0;
tdk = @(wk) min(T, 1/max([a0,C1*Cth]) .* ...
                log((1 + sqrt(1 + 4*a0*Cetm.*wk./Ak(wk)))/2));
fprintf('\n--- Self-triggering ---\n');
fprintf('  tdk(0) = %.6f  (lower bound on inter-event interval)\n', tdk(0));

%% Parameter selection: (sigma, epsilon, delta_sg)
% Finds the triple satisfying the small-gain condition (iv)
% phi(sigma,eps,delta_sg)/(1+lambda) < 1
% using a two-pass grid search.
fprintf('\n=== Parameter selection (sigma, epsilon, delta_sg) ===\n');
assert(phi000/(1+lambda) < 1, ...
    'phi(0,0,0)/(1+lambda) >= 1 even after lambda adjustment.');

phi1val = phi1(1);
phi2val = phi2(1);
M1bar   = M1;
eta     = p;
rho     = r;

% --- Coarse search ---
sig_vec  = linspace(0,    0.5, 30);
eps_vec  = linspace(0,    0.3, 30);
dsg_vec  = linspace(1e-4, 0.4, 30);

best_ratio = Inf;
sig_opt = 0;  eps_opt = 0;  dsg_opt = 1e-4;

for sig  = sig_vec
for eps  = eps_vec
for dsg  = dsg_vec
    if sig > 0 && dsg >= sig, continue; end
    v = compute_phi(sig,eps,dsg,P1,P2,amin,bmin, ...
                    phi1val,phi2val,eta,rho,Cetm,lambda);
    if v/(1+lambda) < best_ratio
        best_ratio = v/(1+lambda);
        sig_opt = sig;  eps_opt = eps;  dsg_opt = dsg;
    end
end
end
end
fprintf('Coarse: sigma=%.4f  eps=%.4f  delta_sg=%.4f  ratio=%.6f\n', ...
        sig_opt,eps_opt,dsg_opt,best_ratio);

% --- Refined search ---
ds = sig_vec(2)-sig_vec(1);
de = eps_vec(2)-eps_vec(1);
dd = dsg_vec(2)-dsg_vec(1);

for sig  = linspace(max(0,sig_opt-ds), sig_opt+ds, 40)
for eps  = linspace(max(0,eps_opt-de), eps_opt+de, 40)
for dsg  = linspace(max(1e-6,dsg_opt-dd), dsg_opt+dd, 40)
    if sig > 0 && dsg >= sig, continue; end
    v = compute_phi(sig,eps,dsg,P1,P2,amin,bmin, ...
                    phi1val,phi2val,eta,rho,Cetm,lambda);
    if v/(1+lambda) < best_ratio
        best_ratio = v/(1+lambda);
        sig_opt = sig;  eps_opt = eps;  dsg_opt = dsg;
    end
end
end
end

% --- Final report ---
phi_opt   = compute_phi(sig_opt,eps_opt,dsg_opt,P1,P2,amin,bmin, ...
                        phi1val,phi2val,eta,rho,Cetm,lambda);
gam_opt   = (1+eps_opt)*abs(eta*rho);
Ttrans    = phi1val + phi2val;
k_eta_opt = abs(eta)*Cetm*(1+exp(dsg_opt*Ttrans))*(1+eps_opt)^2/(1-gam_opt);
k_c_opt   =         Cetm*(1+exp(dsg_opt*Ttrans))*(1+eps_opt)  /(1-gam_opt);
k_bb_opt  = k_eta_opt*k_c_opt/((1-k_eta_opt)*(1-k_c_opt));

fprintf('\n=== Optimal parameters ===\n');
fprintf('  sigma    = %.6f\n', sig_opt);
fprintf('  epsilon  = %.6f\n', eps_opt);
fprintf('  delta_sg = %.6f\n', dsg_opt);
fprintf('\n=== Verification ===\n');
fprintf('  (i)   gamma=(1+eps)|pr| = %.6f < 1       [%s]\n', gam_opt,          tf(gam_opt<1));
fprintf('  (ii)  Cetm=%.4f < Cetmax=%.4f            [%s]\n', Cetm,Cetmax,       tf(Cetm<Cetmax));
fprintf('  (iii) k_eta+k_c = %.4f < 1               [%s]\n', k_eta_opt+k_c_opt, tf(k_eta_opt+k_c_opt<1));
fprintf('        k_bb = %.6f < 1                     [%s]\n', k_bb_opt,          tf(k_bb_opt<1));
fprintf('  (iv)  phi(sig,eps,dsg)/(1+lambda) = %.6f  [%s]\n', phi_opt/(1+lambda),tf(phi_opt/(1+lambda)<1));

%% Parameter vectors for Simulink
param  = [a0, C1, Cth, Cp, Kq, M2, Cetm];
qparam = [M, Omega, T, mu0, Delta];

save('stc_constants.mat', ...
     'a0','C1','Cth','Cp','Kq','M2','Cetm', ...
     'M','Delta','Omega','T','mu0','delta',  ...
     'Cetmax','phi000','lambda','Ccond',      ...
     'sig_opt','eps_opt','dsg_opt',           ...
     'param','qparam');
fprintf('\nAll constants saved to stc_constants.mat\n');

% %% Sensitivity plot: phi/(1+lambda) over (sigma, eps) at delta_sg*
% SIG = linspace(0,0.5,60);
% EPS = linspace(0,0.3,60);
% PHI = zeros(numel(EPS),numel(SIG));
% for is = 1:numel(SIG)
%     for ie = 1:numel(EPS)
%         PHI(ie,is) = compute_phi(SIG(is),EPS(ie),dsg_opt, ...
%                          P1,P2,amin,bmin,phi1val,phi2val, ...
%                          eta,rho,Cetm,lambda) / (1+lambda);
%     end
% end
% figure;
% contourf(SIG,EPS,min(PHI,2),20,'LineColor','none'); colorbar; hold on;
% contour(SIG,EPS,PHI,[1 1],'r-','LineWidth',2);
% plot(sig_opt,eps_opt,'w*','MarkerSize',12,'LineWidth',2);
% xlabel('\sigma','Interpreter','tex','FontSize',12);
% ylabel('\epsilon','Interpreter','tex','FontSize',12);
% title(sprintf('\\phi/(1+\\lambda) at \\delta_{sg}^*=%.4f  (red=boundary)',dsg_opt), ...
%       'Interpreter','tex','FontSize',12);
% legend('\phi/(1+\lambda)','= 1 boundary','Optimum', ...
%        'Location','northeast','Interpreter','tex');
% grid on;

% ================================================================
%  Local functions
% ================================================================
function val = compute_phi(sig,eps,dsg,P1,P2,amin,bmin, ...
                            phi1v,phi2v,eta,rho,Cetm,lambda)
% Evaluates phi(sigma,epsilon,delta_sg)
    gam = (1+eps)*abs(eta*rho);
    if gam >= 1, val = Inf; return; end

    e1 = exp(sig*phi1v);  e2 = exp(sig*phi2v);

    g1 = e1/amin*P1 + abs(eta)*(abs(rho) + e2/bmin*P2);
    g3 = abs(rho) + e2/bmin*P2;
    g2 = abs(rho)*(1+eps)*e1/amin*P1 + abs(rho) + e2/bmin*P2;
    g4 = abs(eta)*(1+eps)*(abs(rho)+e2/bmin*P2) + e1/amin*P1;

    Ttrans = phi1v + phi2v;
    k_eta  = abs(eta)*Cetm*(1+exp(dsg*Ttrans))*(1+eps)^2/(1-gam);
    k_c    =         Cetm*(1+exp(dsg*Ttrans))*(1+eps)  /(1-gam);

    if k_eta>=1 || k_c>=1 || k_eta+k_c>=1, val=Inf; return; end

    denom = (1-gam)*(1-k_eta-k_c);
    if denom <= 0, val = Inf; return; end

    val = (1+eps)*(g1+2*g3) ...
        + (1+eps)*(1+(1+eps)*(g1+g3))*(g2+g4)/denom;
end

function s = tf(cond)
    if cond, s='OK'; else, s='FAIL'; end
end

%% Compute C0 = max{c_alpha0tilde, c_beta0tilde, c_alpha0hat, c_beta0hat}
% Uses sig_opt, eps_opt, dsg_opt from the parameter selection above.
% All k-constants follow from the formulas in the paper.

fprintf('\n=== Computing C0 ===\n');

% Shorthand evaluated at optimal parameters
sig  = sig_opt;
eps  = eps_opt;
dsg  = dsg_opt;

e1   = exp(sig*phi1val);   % exp(sigma*phi1(1))
e2   = exp(sig*phi2val);   % exp(sigma*phi2(1))
gam  = (1+eps)*abs(eta*rho);            % gamma(eps)
igam = 1/(1-gam);

% gamma_1, gamma_2, gamma_3, gamma_4
g1 = e1/amin*P1 + abs(eta)*(abs(rho) + e2/bmin*P2);
g3 = abs(rho) + e2/bmin*P2;
g2 = abs(rho)*(1+eps)*e1/amin*P1 + abs(rho) + e2/bmin*P2;
g4 = abs(eta)*(1+eps)*(abs(rho)+e2/bmin*P2) + e1/amin*P1;

% k_eta, k_c  (eqs. kp, kc)
Ttrans = phi1val + phi2val;
k_eta  = abs(eta)*Cetm*(1+exp(dsg*Ttrans))*(1+eps)^2/(1-gam);
k_c    =          Cetm*(1+exp(dsg*Ttrans))*(1+eps)  /(1-gam);
ike    = 1/(1-k_eta);
ikc    = 1/(1-k_c);

% k_bb = k_{betahat betahat}  (eq. khatbeta)
k_bb   = ike * ikc * k_c * k_eta;
ikbb   = 1/(1-k_bb);

% ---- k_{betahat0 betahat}  (eq. khatbetabeta0)
k_b0b  = igam * ikc * e2 ...
        * (1 + k_c*ike*(1+eps)*abs(eta));

% ---- k_{alphahat0 betahat}  (eq. kalphahat0betahat)
k_a0b  = ikc * igam * e1 ...
        * (abs(rho)*(1+eps) + ike*k_c);

% ---- k_{alphatilde0 betahat}  (eq. ktildealpha0betahat)
k_ta0b = igam * ikc * (1+eps) * M1bar ...
        * (g2 + ike*g4*k_c);

% ---- k_{betatilde0 betahat}  (eq. ktildebeta0betahat)
k_tb0b = igam * ikc * (1+eps) * M1bar ...
        * (g2 + ike*g4*k_c);

% ---- k_{qbetahat}  (eq. kqbetahat)
k_qb   = ikc * igam * (1+eps) ...
        * (1+(1+eps)*(g1+g3)) ...
        * (ike*k_c*g4 + g2);

% ---- k_{alphahat0 alphahat}  (eq. kalphahat0alphahat)
k_a0a  = ike * igam * e1 ...
        * (1 + ikc*k_eta*abs(rho)*(1+eps));

% ---- k_{betahat0 alphahat}  (eq. kbetahat0alphahat)
k_b0a  = ike * igam * e2 ...
        * (abs(eta)*(1+eps) + ikc*k_eta);

% ---- k_{alphatilde0 alphahat}  (eq. ktildealpha0alphahat)
k_ta0a = igam * ike * (1+eps) * M1bar ...
        * (g4 + ikc*g2*k_eta);

% ---- k_{betatilde0 alphahat}  (eq. ktildebeta0alphahat)
k_tb0a = ike * igam * (1+eps) * M1bar ...
        * (g4 + ikc*g2*k_eta);

% ---- k_{qalphahat}  (eq. kqalphahat)
k_qa   = ike * igam ...
        * (1+(1+eps)*(g1+g3)) ...
        * (1+eps) * (ikc*k_eta*g2 + g4);

% ================================================================
% c constants
% ================================================================
c_ta0 = M1bar + ikbb*(k_ta0b + k_ta0a);          % c_{alphatilde0}
c_tb0 = M1bar + e2 + ikbb*(k_tb0b + k_tb0a);     % c_{betatilde0}
c_a0  = ikbb*(k_a0b  + k_a0a);                   % c_{alphahat0}
c_b0  = ikbb*(k_b0b  + k_b0a);                   % c_{betahat0}

C0 = max([c_ta0, c_tb0, c_a0, c_b0]);

fprintf('  k_bb       = %.6f  (< 1: %s)\n',  k_bb,  tf(k_bb<1));
fprintf('  c_alpha0~  = %.6f\n', c_ta0);
fprintf('  c_beta0~   = %.6f\n', c_tb0);
fprintf('  c_alpha0^  = %.6f\n', c_a0);
fprintf('  c_beta0^   = %.6f\n', c_b0);
fprintf('  C0         = %.6f\n', C0);

% Add C0 to saved constants
save('stc_constants.mat', 'C0', 'c_ta0','c_tb0','c_a0','c_b0', ...
     'k_bb','k_eta','k_c','k_qb','k_qa', '-append');
fprintf('  C0 appended to stc_constants.mat\n');
%% G0
G0 = C0 * (1 - phi000/(1+lambda))^(-1);
fprintf('G0 = %.6f\n', G0);