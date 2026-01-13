%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Galerkin system of ODE's and reconstruct velocity field
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
load('galerkin_model_cntrl_2_36_2_Lx_2_Lz_1.mat');
toc

nmodes               = size(phi,5);
nalphas              = 2; ngammas = 2;
Lx                   = (max(max(max(X)))+X(1,2)-X(1,1))/pi;
Lz                   = (max(max(max(Z)))+Z(1,1,2)-Z(1,1,1))/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F                   = zeros(nmodes,1);      % Body forcing term; zero for standandard Couette flow (w/o pressure gradient)
Re                  = 1000;                 % Reynolds number for ROM simulation
L0                  = L;
L                   = L0 + Re*L2;           % Complete linear term (diffusion + energy extraction from base flow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ROM simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan               = 0:0.5:1000;           % Time vector      
q0                  = 0.1*(rand(nmodes,1)); % Initial condition
options             = odeset('RelTol',1e-6);% Options for ODE solver
disp(['Nmodes = ' int2str(nmodes)])
disp('Simulation');
% Time iontegration
tic;
[t,q]               = ode45(@(t,X) galerkinsys(t,X,L,QQ,F/Re,Re),tspan,q0,options);
toc

% Plot mode coefficients
figure
plot(t,q)
xlabel('$t$')
ylabel('$a_i$')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct velocity field
% u(x,y,z,t) =\sum q(t)\phi(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt                        = length(t); 
posend                    = floor(nt/2):1:nt;
posstats                  = floor(nt/2):1:nt;
qfield                    = zeros(size(phi,1),size(phi,2),size(phi,3),size(phi,4),length(posend));
for j=1:length(posstats)
    disp(['Statistics, ' int2str(j) '/' int2str(length(posstats))]);
    qfield(:,:,:,:,j)     = phi0;
    for i=1:nmodes
        qfield(:,:,:,:,j) = qfield(:,:,:,:,j)+squeeze(phi(:,:,:,:,i))*q(posend(j),i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute flow statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos                     = find(kvec(:,1)==0&kvec(:,2)==0); % Find mean flow modes (with zero streamwise and spanwise wavenumbers)
U0                      = phi0(:,1,1,1);                   % Base flow
x                       = squeeze(X(1,:,1));
y                       = squeeze(Y(:,1,1));
z                       = squeeze(Z(1,1,:));

% Compute mean flow
for i=pos.'
    U0                  = U0+squeeze(phi(:,1,1,1,i))*mean(q(posend,i));
end

qmean                   = mean(qfield,5);
qprime                  = qfield - qmean;

uv                      = qprime(:,:,:,1,:).*qprime(:,:,:,2,:);

q_rms                   = squeeze(sqrt(mean(mean(mean(qprime.^2,5),3),2)));

uv_rms                  = mean(uv,5);
uv_rms                  = mean(mean(uv_rms,3),2);

figure
plot(y,U0)
hold on
xlabel('$y$','Interpreter','latex')
ylabel('$\overline{U}$','Interpreter','latex')

figure
plot(y,q_rms)
xlabel('$y$','Interpreter','latex')
ylabel('$\left<u,v,w\right>_{xzt}$','Interpreter','latex')
