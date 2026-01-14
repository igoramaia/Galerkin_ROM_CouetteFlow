%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing optimisation using ROMs                   
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

nmodes               = size(phi,5);
nalphas              = 2; ngammas = 2;
Lx                   = (max(max(max(X)))+X(1,2)-X(1,1))/pi;
Lz                   = (max(max(max(Z)))+Z(1,1,2)-Z(1,1,1))/pi;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Stokes modes for u and w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos                 = find(kvec(:,1)==0&kvec(:,2)==0);
phi_u_pos           = phi(:,:,:,1,pos);
phi_v_pos           = phi(:,:,:,2,pos);
phi_w_pos           = phi(:,:,:,3,pos);
y                   = squeeze(Y(:,1,1));

target_u            = squeeze(phi_u_pos(:,:,:,:));
target_w            = squeeze(phi_w_pos(:,:,:,:));

[~, id_x]           = find(max(max(max(target_u)))>1e-2);
[~, id_w]           = find(max(max(max(target_w)))>1e-2);

% Making sure forcing modes only have u component of velocity
if ~isempty(find(id_x==id_w))
    error('Mixed Stokes modes for u and w')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run ROM simulation (first iteration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F                   = zeros(nmodes,1);      % Body forcing term; zero for standandard Couette flow (w/o pressure gradient)
Re                  = 1000;                 % Reynolds number for ROM simulation
L0                  = L;
L                   = L0 + Re*L2;           % Complete linear term (diffusion + energy extraction from base flow)

% Simulation
tspan               = 0:0.5:2500;
t_targ1             = tspan(1);
t_targ2             = tspan(3001);
q0                  = 0.1*(rand(nmodes,1)); % Initial condition
options             = odeset('RelTol',1e-6);% Options for ODE solver
disp(['Nmodes = ' int2str(nmodes)])
disp('Simulation');
% Time iontegration
tic;
[t,q_ref]           = ode45(@(t,X) galerkinsys(t,X,L,QQ,F/Re,Re),tspan,q0,options);
toc

% Compute functional
J_ref               = compute_J(t,q_ref,pos);

% Select turbulent trajectory for control (after initial transient)
std_q               = std((q_ref-mean(q_ref,2)),1,2);
id                  = find(std_q<1e-2);
if isempty(id)
    nt              = length(t);
else
    nt              = id(1);
end
pos_init            = floor(nt/2);
q0_turb             = q_ref(pos_init,:); % Initial condition for control (turbulent trajectory)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining forcing basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_cntrl_x                      = zeros(Ny,Nx,Nz,nmodes);
phi_cntrl                        = zeros(Ny,Nx,Nz,3,nmodes);

phi_cntrl_x(:,:,:,pos(id_x))     = phi(:,:,:,1,pos(id_x)); % Streamwise forcing in the u component alpha=0 and beta=0
phi_cntrl(:,:,:,1,:)             = phi_cntrl_x; % Assigning the forcing modes

% Forcing operator
FF                               = zeros(nmodes,nmodes);
for i=1:nmodes
    for j=1:nmodes
      FF(i,j)                    = inprod_s(phi(:,:,:,:,i),phi_cntrl(:,:,:,:,j),Wf);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_stokes                = length(id_x);
nit                     = 10;
epsilon                 = 3e-4;
stokes_ind              = length(id_x);
dx                      = 1e-2;

% Initial Guess
b_i                     = zeros(nmodes,nit);  
b_i(pos(id_x))          = 0;

J                       = zeros(1,nit);

x2                      = [0 1];
y2                      = [1.5 J_ref];
p                       = plot(x2,y2,'o');
xlim([0 nit])
ylim([0 4])
xlabel('Iteration')
ylabel('$J$')
p(1).XDataSource        = 'x2';
p(1).YDataSource        = 'y2';
id_laminar              = zeros(1,nit);
parpool(6);
ns_index                = 2:2:n_stokes;
ind_F                   = pos(id_x(ns_index));
for i=1:nit
    % Compute f(x)
    tic
    [t,q]              = ode45(@(t,X) galerkinsys_forc_tramp(t,X,L,QQ,FF*b_i(:,i),Re,t_targ1,t_targ2),tspan,q0_turb,options);
    toc
    J(i)               = compute_J(t,q,pos);
    std_q              = std((q-mean(q,2)),1,2);
    id                 = find(std_q<1e-3);
    if isempty(id)
        id_laminar(i)  = length(t);
    else
        id_laminar(i)  = id(1);  
    end
    if i==1
        baux                                    = b_i(:,i);
        ref                                     = J(i);
        parfor ns=1:length(ns_index)
            x_dx                                = baux;
            x_dx(ind_F(ns))                     = x_dx(ind_F(ns))+dx; % x+dx
            % Compute f(x+dx)
            [t,q]                               = ode45(@(t,X) galerkinsys_forc_tramp(t,X,L,QQ,FF*x_dx,Re,t_targ1,t_targ2),tspan,q0_turb,options);
            J_x_dx                              = compute_J(t,q,pos);
            %Compute grad
            dtau_dF                             = (J_x_dx - ref)/dx
            auxvar(ns)                          = epsilon*dtau_dF
            ns
        end
        b_i(ind_F,i+1)                          = b_i(ind_F,i)- auxvar.';
    else
        if id_laminar(i)<length(t)
            epsilon                             = 1e-5;
            dx                                  = 1e-4;
        else 
            epsilon                             = 3e-4;
            dx                                  = 1e-2;
        end
        baux                                    = b_i(:,i);
        ref                                     = J(i);
        parfor ns=1:length(ns_index)
            x_dx                                = baux;
            x_dx(ind_F(ns))                     = x_dx(ind_F(ns))+dx; % x+dx
            % Compute f(x+dx)
            [t,q]                               = ode45(@(t,X) galerkinsys_forc_tramp(t,X,L,QQ,FF*x_dx,Re,t_targ1,t_targ2),tspan,q0_turb,options);
            J_x_dx                              = compute_J(t,q,pos);
            %Compute grad
            dtau_dF                             = (J_x_dx - ref)/dx;
            auxvar(ns)                          = epsilon*dtau_dF;
            ns
        end
        auxvar
        b_i(ind_F,i+1)                          = b_i(ind_F,i) -auxvar.'; % Go in the direction opposite to the gradient
    end
    x2(i)                                       = i;
    y2(i)                                       = J(i);
    ylim([0 max(J*1.1)])
    refreshdata
    drawnow
    if id_laminar(i) <nt && abs(J(i)-J(i-1))<1e-3
        break
    end
end

b_optim = b_i(:,i);
save(['control_rom_',num2str(nmodes),'_Re',num2str(Re),'.mat'],'b_optim','q0_turb','Re','pos','id_x','J', 'phi_u_pos', 'FF');