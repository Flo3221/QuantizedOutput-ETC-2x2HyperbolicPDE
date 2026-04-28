%% =========================
%  States and observer
%% =========================
TT = out.tout;                 % general simulation time (Nt x 1)
UV = out.uv.Data;              % full state data (Nt x 2*Nx)
%tauks=out.tauk.Data;

% --- define Nx safely ---
if exist('xg','var') && ~isempty(xg)
    Nx = xg;
else
    Nx = size(UV,2)/2;
end

UU = UV(:,1:Nx);               % u data
VV = UV(:,Nx+1:2*Nx);          % v data

UVhat = out.uvhat.Data;        % observer data (Nt x 2*Nx)
UUhat = UVhat(:,1:Nx);         % u hat
VVhat = UVhat(:,Nx+1:2*Nx);    % v hat

% --- Ensure your plotting vectors exist and are numeric vectors ---
t = TT(:);                     % use your TT for plotting

% If you have x1 in workspace (as in your commented surf lines), use it:
if exist('x1','var') && ~isempty(x1)
    x = x1(:);
elseif exist('x','var') && ~isempty(x)
    x = x(:);
else
    % fallback if x not available (adapt if your domain is different)
    x = linspace(0,1,Nx).';
end
%% =========================
%  Plots 
%% =========================
pstep = 10;   % spacing between spatial grid points in plot
tstep = 10;   % spacing between time steps in plot
idxX = 1:pstep:length(x);
idxT = 1:tstep:length(t);
% Create meshgrid (NO 'end' keyword)
[Y, Z] = meshgrid(x(idxX), t(idxT));
% -------------------------
% Figure 1: v(t,x)
% -------------------------
figure(1); clf
V_plot = VV(idxT, idxX);
mesh(Y, Z, V_plot, 'LineWidth', 1, 'EdgeColor', 'black');
view(83,10);
ax = gca;
ax.FontSize = 40;
xlabel('$x$','Interpreter','latex','LineWidth',7,'FontSize',50);
ylabel('$t$','Interpreter','latex','LineWidth',7,'FontSize',50);
zlabel('$v(t,x)$','Interpreter','latex','LineWidth',7,'FontSize',50);
set(gca,'FontSize',40);
% grid on
%% =========================
%  Phi(t)
%% =========================
Nt = length(t);
Phi = zeros(Nt,1);

for k = 1:Nt
    Phi(k) = norm(UU(k,:),inf) + norm(VV(k,:),inf) + norm(UUhat(k,:),inf) + norm(VVhat(k,:),inf);
end
tmu=out.mu.Time;
mu=out.mu.Data;

figure(2); clf
plot(t, Phi,'k', 'LineWidth', 5);
% grid on
xlabel('$t$','Interpreter','latex','FontSize',50);
ylabel('$\Psi(t)$','Interpreter','latex','FontSize',50);
set(gca,'FontSize',40);
% --- Get signals (timeseries) ---
tnom = out.Unom.Time;
Unom = out.Unom.Data;

tetc = out.Uetc.Time;
Uetc = out.Uetc.Data;

% --- Plot ---
figure(3); clf
plot(tnom, Unom,'r', 'LineWidth', 7); hold on
stairs(tetc, Uetc,'k', 'LineWidth', 7);  
 grid on
xlabel('$t$','Interpreter','latex','FontSize',50)
ylabel('','Interpreter','latex','FontSize',50)
 legend({'$U_{\mathrm{nom}}(t)$','$U_d(t)$'}, 'Interpreter','latex','Location','best')
set(gca,'FontSize',40)
%title('Nominal vs event-triggered control','Interpreter','latex')

% --- Get signals (timeseries) ---
ty  = out.y.Time;
y   = out.y.Data;

tyq = out.yq.Time;
yq  = out.yq.Data;

tmmu = out.mu.Time;     % out.mu = M*mu(t) (your block output "Mmu")
Mmu  = out.mu.Data;

% --- Plot ---
figure(4); clf
plot(ty,  y,  'r', 'LineWidth', 5); hold on
stairs(tyq, yq, 'k', 'LineWidth', 5);
 plot(tmmu,  Mmu,'k:', tmmu,-Mmu, 'k:', 'LineWidth', 5);
 % plot(tmmu, -Mmu, 'k', 'LineWidth', 5);   % optional but usually meaningful

% grid on
xlabel('$t$','Interpreter','latex','FontSize',50)
ylabel('','Interpreter','latex','FontSize',50)
legend({'$y(t)$','$y_q(t)$','$\pm M\nu(t)$'},'Interpreter','latex','Location','best')
set(gca,'FontSize',40)
% title('Output, quantized output, and range','Interpreter','latex')

% tdk=zeros(1,Nt);
% for k = 1:Nt
%     wk= norm(UU(k,:),inf) + norm(VV(k,:),inf);
%     tdk(k)=  1/max([a0,C1*Cth])*log((1 + sqrt(1 + 4*a0*Cetm*wk./Ak(wk)))/2);
% 
% end
% figure(5); clf
% plot(t, tdk,'k', 'LineWidth', 5);
% % grid on
% xlabel('$t$','Interpreter','latex','FontSize',50);
% ylabel(' ','Interpreter','latex','FontSize',50);
% legend({'$\frac{1}{\max\{a_0,C_1C\}} \ln\left(\frac{1+\sqrt{1+4\frac{a_0 C_{\rm etm} w(\tau_k)}{A_k}}}{2}\right) $'},'Interpreter','latex','Location','best')
% set(gca,'FontSize',40);

% % Nombre total de triggers ETC
% N_triggers_etc = out.count_etc.Data(end);
% fprintf('ETC - Nombre de triggers : %d\n', N_triggers_etc);
% 
% % Instants de triggering
% count_vec     = out.count_etc.Data;
% time_vec      = out.count_etc.Time;
% trigger_times = time_vec(diff([0; count_vec]) > 0);
% fprintf('ETC - Inter-event moyen : %.4f s\n', mean(diff(trigger_times)));
% fprintf('ETC - Inter-event min   : %.4f s\n', min(diff(trigger_times)));