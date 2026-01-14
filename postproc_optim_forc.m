%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-processing of forced system                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all;
clc
addpath(genpath('aux_functions'));

set(groot,'defaulttextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultColorbarTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultfigureunits','inches');
set(0,'defaultfigureposition',[7 7 6 4.5])
set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultFigureColor','remove')
set(0,'DefaultFigureColor',[1 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading model');
tic;
load('../galerkin_model_cntrl_2_36_2_Lx_2_Lz_1.mat');
toc
Re                   = 1000;
L0                   = L;
L                    = L0 + Re*L2;           % Complete linear term (diffusion + energy extraction from base flow)
y                    = squeeze(Y(:,1,1));
nmodes               = size(phi,5);
nalphas              = 2; ngammas = 2;
Lx                   = (max(max(max(X)))+X(1,2)-X(1,1))/pi;
Lz                   = (max(max(max(Z)))+Z(1,1,2)-Z(1,1,1))/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load and reconstruct optimal forcing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['control_rom_',num2str(nmodes),'_Re',num2str(Re),'.mat']);
n_stokes                = length(id_x);

fx                      = zeros(Ny,Nx,Nz); 
for ns=2:2:n_stokes
    fx                  = fx + phi_u_pos(:,:,:,id_x(ns))*b_optim(pos(id_x(ns)));
end

fx                      = mean(mean(fx,2),3);

figure
plot(y,fx)
xlabel('$y$','Interpreter','latex')
ylabel('$F(y)$','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Energy analysis (turning forcing on and off)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tspan                  = 0:0.5:4000;
idt1                   = 3001;      % Time index 1: turn forcing on
idt2                   = 6001;      % Time index 2: turn forcing off
t_targ1 = tspan(idt1);
t_targ2 = tspan(idt2);
options                = odeset('RelTol',1e-6);% Options for ODE solver

[t,q_fnf]              = ode45(@(t,X) galerkinsys_forc_tramp(t,X,L,QQ,FF*b_optim,Re,t_targ1,t_targ2),tspan,q0_turb,options);

% Coefficients
figure
plot(t,q_fnf)
hold on
plot(t_targ1*ones(1,100), linspace(-0.5,0.5,100),'--k')
plot(t_targ2*ones(1,100),linspace(-0.5,0.5,100),'--k')
xlabel('$t$','Interpreter','latex')
ylabel('$a_i$','Interpreter','latex')

% Input and dissipation
I                   = sum(q_fnf'.*((L2*Re)*q_fnf')); % Energy input: extraction of energy from base flow
I2                  = sum((FF*b_optim*(0.5+0.5*tanh(0.1*(tspan-t_targ1)))-FF*b_optim*(0.5+0.5*tanh(0.1*(tspan-t_targ2))))*Re.*q_fnf.'); % Energy input: body forcing
D                   = -sum(q_fnf'.*(L0*q_fnf')); % Energy dissipation
E = sum(q_fnf.^2,2)/2;
I = I.'+1;
D = D.'+1;
I2 = I2.'+1;
figure;
plot(t,I,t,I2,t,D);
xlabel('$t$','Interpreter','latex')
legend('$I_{1}$','$I_{2}$','$D$')

figure
plot(linspace(0,10,100),linspace(0,10,100),'--r')
hold on
h=plot(I,D,'k','linewidth',1)
[X,Y]=deal(h.XData, h.YData);
qq=(1:0.5:numel(X))';
X=interp1(X(:),qq); Y=interp1(Y(:),qq);
hold on 
quiver(X(1:end-1),Y(1:end-1), diff(X), diff(Y), 0,'k')
hold off
xlabel('$I$','Interpreter','latex')
ylabel('$D$','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forced laminar solution (da_i/dt=0 and throw away quadratic terms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_lam           = -(L/Re)\FF*b_optim;
U0_lam          = phi0(:,1,1,1);
ns_ind = [2:2:n_stokes];
for i=pos(id_x(ns_ind))
    U0_lam      = U0_lam+squeeze(phi(:,1,1,1,i))*a_lam(i);
end

figure
plot(y,U0_lam)
xlabel('$y$','Interpreter','latex')
ylabel('$U_{lam}$','Interpreter','latex')
set(gca,'FontSize',20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stability of the forced flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A           = L/Re;
for i=1:nmodes
    A       = A + squeeze(Q(:,i,:))*a_lam(i) + squeeze(Q(:,:,i))*a_lam(i);
end
[V,lambda] = eig(A);
lambda     = diag(lambda);
[~,posl]=sort(real(lambda),'descend'); lambda = lambda(posl);
V          = V(:,posl);
lcount     = find(real(lambda)>1e-4);
ninst      = length(lcount);

A_couette  = L/Re;
[V_couette,lambda_couette] = eig(A_couette);
lambda_couette = diag(lambda_couette);
[~,posl]   = sort(real(lambda_couette),'descend'); lambda_couette = lambda_couette(posl);
V_couette  = V_couette(:,posl);

figure
plot(real(lambda_couette),imag(lambda_couette),'o')
hold on
plot(real(lambda),imag(lambda),'+')
xlabel('$\lambda_r$','Interpreter','latex')
ylabel('$\lambda_i$','Interpreter','latex')
legend('Laminar Couette','Forced flow','box','off','location','northeast')
xlim([-0.25 0.15])
set(gca,'FontSize',20)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transient growth based on forced flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n                   = ceil(length(tspan)/2);
cont                = 1;
for i=1:10:n
    [U,S,V]         = svd(expm(A*tspan(i)));
    [U_c,S_c,V_c]   =   svd(expm(A_couette*tspan(i)));
    [~,posl]=sort(real(S),'descend'); S     = S(posl);
    [~,posl]=sort(real(S_c),'descend'); S_c = S_c(posl);
    G(cont)         = max(max(S));
    G_c(cont)       = max(max(S_c));
    cont            = cont+1;
end

figure
plot(tspan(1:10:n),G_c.^2,'o')
hold on
plot(tspan(1:10:n),G.^2,'+')
xlabel('$t$','Interpreter','latex')
ylabel('$G(t)$','Interpreter','latex')
legend('Laminar Couette','Forced flow','box','off','location','northeast')
set(gca,'FontSize',20)
