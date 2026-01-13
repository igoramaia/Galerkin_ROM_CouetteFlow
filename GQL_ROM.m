%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform generalised-quasilinear approximations (GQL) within ROM framework                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all;
clc
addpath(genpath('aux_functions'));
addpath(genpath('aux_files'));

set(groot,'defaulttextInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaultColorbarTickLabelInterpreter','latex'); 
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultfigureunits','inches');
set(0,'defaultfigureposition',[7 7 6 4.5])
set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',18);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading model');
tic;
load('galerkin_model_cntrl_2_24_2_Lx_2_Lz_1.mat');
toc

nmodes               = size(phi,5);
nalphas              = 2; ngammas = 2;
Lx                   = (max(max(max(X)))+X(1,2)-X(1,1))/pi;
Lz                   = (max(max(max(Z)))+Z(1,1,2)-Z(1,1,1))/pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Choose GQL formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driven_gql           = false; % Use 1 for Driven-GQL formulation and 0 for standard GQL
                          % See the details in Maia & Cavalieri (Theor. Comp. Fluid Dyn. 38, 313-330 2024)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F                   = zeros(nmodes,1);      % Body forcing term; zero for standandard Couette flow (w/o pressure gradient)
Re                  = 1000;                 % Reynolds number for ROM simulation
L0                  = L;
L                   = L0 + Re*L2;           % Complete linear term (diffusion + energy extraction from base flow)

tspan               = 0:0.5:4000;           % Time vector      
q0                  = 0.1*(rand(nmodes,1)); % Initial condition
options             = odeset('RelTol',1e-6);% Options for ODE solver
disp(['Nmodes = ' int2str(nmodes)])
disp('Simulation');
% Time iontegration
tic;
[t,q_ref]           = ode45(@(t,X) galerkinsys(t,X,L,QQ,F/Re,Re),tspan,q0,options);
toc

% Compute Retau
I_ref               = sum(q_ref'.*((L2*Re)*q_ref'))+1;
Retau_ref           = mean(Re.*(sqrt(I_ref(end/2:end)/Re)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct velocity field and compute statistics (reference simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt                          = length(t); 
posend                      = floor(nt/2):1:nt;
posstats                    = floor(nt/2):1:nt;
qfield_ref                  = zeros(size(phi,1),size(phi,2),size(phi,3),size(phi,4),length(posend));
for j=1:length(posstats)
    disp(['Statistics, ' int2str(j) '/' int2str(length(posstats))]);
    qfield_ref(:,:,:,:,j)     = phi0;
    for i=1:nmodes
        qfield_ref(:,:,:,:,j) = qfield_ref(:,:,:,:,j)+squeeze(phi(:,:,:,:,i))*q_ref(posend(j),i);
    end
end

pos                     = find(kvec(:,1)==0&kvec(:,2)==0); % Find mean flow modes (with zero streamwise and spanwise wavenumbers)
U0_ref                  = phi0(:,1,1,1);                   % Base flow
x                       = squeeze(X(1,:,1));
y                       = squeeze(Y(:,1,1));
z                       = squeeze(Z(1,1,:));

% Compute mean flow
for i=pos.'
    U0_ref              = U0_ref+squeeze(phi(:,1,1,1,i))*mean(q_ref(posend,i));
end

qmean_ref               = mean(qfield_ref,5);
qprime_ref              = qfield_ref - qmean_ref;

uv_ref                  = qprime_ref(:,:,:,1,:).*qprime_ref(:,:,:,2,:);

q_rms_ref               = squeeze(sqrt(mean(mean(mean(qprime_ref.^2,5),3),2)));

uv_rms_ref              = mean(uv_ref,5);
uv_rms_ref              = mean(mean(uv_rms_ref,3),2);

ytau_ref                = (y+1)*Retau_ref;

figure
plot(ytau_ref,(q_rms_ref(:,1)*(Re/Retau_ref)).^2,'k-')
xlim([0 Retau_ref])
xlabel('$y^+$','Interpreter','latex')
ylabel('$\left<u^+u^+\right>_{xz}$','Interpreter','latex')

figure
plot(ytau_ref,(q_rms_ref(:,2)*(Re/Retau_ref)).^2,'k-')
xlim([0 Retau_ref])
xlabel('$y^+$','Interpreter','latex')
ylabel('$\left<u^+u^+\right>_{xz}$','Interpreter','latex')

figure
plot(ytau_ref,uv_rms_ref(:,1)*(Re/Retau_ref).^2,'k-')
xlim([0 Retau_ref])
xlabel('$y^+$','Interpreter','latex')
ylabel('$\left<u^+u^+\right>_{xz}$','Interpreter','latex')

figure
plot(ytau_ref,(q_rms_ref(:,3)*(Re/Retau_ref)).^2,'k-')
xlim([0 Retau_ref])
xlabel('$y^+$','Interpreter','latex')
ylabel('$\left<u^+u^+\right>_{xz}$','Interpreter','latex') 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GQL with controllability criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('cntrl_eigenvalues_2_24_2.mat')       % Load controllability eigenvalues

Ny             = size(phi,1);
[y,A]          = chebdif(Ny,2);
[~,W]          = clencurt(Ny-1);
Wcol           = repmat(W.',Nx*Nz,1);
W              = diag(W(1:end));

thresh         = 0.5; % Controllability threshold typically in the range [1e-1 1e3];
for n = 1:length(thresh)

    % Non-zero wavenumber modes
    posh_aux1   = find(~(kvec(:,1)==0*2*pi/Lx & abs(kvec(:,2))==0*2*pi/Lz));
    posh1       = posh_aux1(lvec(posh_aux1)<=thresh(n)); % Small controllability mode set
    posl1       = posh_aux1(lvec(posh_aux1)>thresh(n));  % Large controllability mode set

    % Proportionality factor
    fac         = max(lvec(posh_aux1))/thresh(n);

    % Mean-flow modes (use the same factor for the mean-flow modes)
    posh_aux2   = find((kvec(:,1)==0*2*pi/Lx & abs(kvec(:,2))==0*2*pi/Lz));
    posh2       = posh_aux2(lvec(posh_aux2)<=max(lvec(posh_aux2))/fac); % Small controllability mode set
    posl2       = posh_aux2(lvec(posh_aux2)>max(lvec(posh_aux2))/fac);  % Large controllability mode set

    % Large- and small-controllability mode sets
    posh        = [posh1;posh2];
    posl        = [posl1;posl2];

    % Corresponding linearisation parameter
    T(n)                       = size(posh,1)/nmodes;

    % Discard triadic interactions
    Q2 = Q;

    if driven_gql
        Q2(posh, posh, posh)    = 0;
        disp('Driven GQL')
    else
        Q2(posh, posh, posh)    = 0;
        Q2(posh, posl, posl)    = 0;
        Q2(posl, posl, posh)    = 0;
        Q2(posl, posh, posl)    = 0;
        disp('Standard GQL')
    end


    % Flatten modified quadratic term
    QQ                          = zeros(nmodes,nmodes^2);
    for i=1:nmodes
        Qaux                    = squeeze(Q2(i,:,:)); Qaux = Qaux(:).';
        QQ(i,:)                 = (Qaux);
    end

    QQ                          = sparse((abs(QQ)>1e-8).*QQ);

    % Run simulation of GQL system
    q0                          = 0.1*(rand(nmodes,1));
    options                     = odeset('RelTol',1e-6);
    disp(['Nmodes = ' int2str(nmodes)])
    disp('Simulation');
    tic;
    [t,q]                       = ode45(@(t,X) galerkinsys(t,X,L,QQ,F/Re,Re),tspan,q0,options);
    toc

    % Compute statistics (GQL system)
    nt                          = length(t);
    posend                      = floor(nt/2):1:nt;
    posstats                    = floor(nt/2):1:nt;
    qfield                      = zeros(size(phi,1),size(phi,2),size(phi,3),size(phi,4),length(posend));
    for j=1:length(posstats)
        qfield(:,:,:,:,j)       = phi0;
        for i=1:nmodes
            qfield(:,:,:,:,j)   = qfield(:,:,:,:,j)+squeeze(phi(:,:,:,:,i))*q(posend(j),i);
        end
    end

    % Mean flow
    pos                         = find(kvec(:,1)==0&kvec(:,2)==0);
    U0                          = phi0(:,1,1,1);
    y                           = squeeze(Y(:,1,1));

    for i=pos.'
        U0                      = U0+squeeze(phi(:,1,1,1,i))*mean(q(posend,i));
    end

    qmean                       = mean(qfield,5);
    qprime                      = qfield - qmean;

    uv                          = qprime(:,:,:,1,:).*qprime(:,:,:,2,:);

    q_rms                       = squeeze(sqrt(mean(mean(mean(qprime.^2,5),3),2)));

    uv_rms                      = mean(uv,5);
    uv_rms                      = mean(mean(uv_rms,3),2);

    % Compare original and modified systems
    ci = 0.15;
    ls = '--';
    figure
    plot(y,q_rms_ref,'k-')
    hold on
    plot(y,q_rms,ls, 'Color',[ci ci ci])
    xlabel('$y$','Interpreter','latex')
    ylabel('$u'', v'', w''$','Interpreter','latex')

    figure
    plot(y,U0_ref, 'k')
    hold on
    plot(y,U0,ls, 'Color',[ci ci ci])
    xlabel('$y$','Interpreter','latex')
    ylabel('$\bar{U}$','Interpreter','latex')

    % Compute Retau
    I                       = sum(q'.*((L2*Re)*q'))+1;
    figure; plot(t,Re.*(sqrt(I/Re)))
    Retau(n)                = mean(Re.*(sqrt(I(end/2:end)/Re)));

    aux                     = Re.*(sqrt(I(end/2:end)/Re))-mean(Re.*(sqrt(I(end/2:end)/Re)));
    it                      = sqrt(mean(aux.^2));

    ytau_ref= (y+1)*Retau_ref;
    ytau= (y+1)*Retau(n);

    ci = 0.25;
    ls = '--';
    figure
    plot(ytau_ref,(q_rms_ref(:,1)*(Re/Retau_ref)).^2,'k-')
    hold on
    plot(ytau,(q_rms(:,1)*(Re/Retau(n))).^2,ls, 'Color',[ci(n) ci(n) ci(n)])
    xlim([0 Retau_ref])
    xlabel('$y^+$','Interpreter','latex')
    ylabel('$\left<u^+u^+\right>_{xz}$','Interpreter','latex')

    figure
    plot(ytau_ref,uv_rms_ref(:,1)*(Re/Retau_ref).^2,'k-')
    hold on
    plot(ytau,uv_rms(:,1)*(Re/Retau(n)).^2,ls, 'Color',[ci(1) ci(1) ci(1)])
    xlim([0 Retau_ref])
    xlabel('$y^+$','Interpreter','latex')
    ylabel('$-\left<u^+v^+\right>_{xz}$','Interpreter','latex')

    figure
    plot(ytau_ref,(q_rms_ref(:,2)*(Re/Retau_ref)).^2,'k-')
    hold on
    plot(ytau,(q_rms(:,2)*(Re/Retau(n))).^2,ls, 'Color',[ci(n) ci(n) ci(n)])
    xlim([0 Retau_ref])
    xlabel('$y^+$','Interpreter','latex')
    ylabel('$\left<v^+v^+\right>_{xz}$','Interpreter','latex')

    figure
    plot(ytau_ref,(q_rms_ref(:,3)*(Re/Retau_ref)).^2,'k-')
    hold on
    plot(ytau,(q_rms(:,3)*(Re/Retau(n))).^2,ls, 'Color',[ci(n) ci(n) ci(n)])
    xlim([0 Retau_ref])
    xlabel('$y^+$','Interpreter','latex')
    ylabel('$\left<w^+w^+\right>_{xz}$','Interpreter','latex')

    figure
    semilogx(ytau,(U0+1)*(Re/Retau(n)),ls, 'Color',[ci(n) ci(n) ci(n)])
    hold on
    semilogx(ytau_ref,(U0_ref+1)*(Re/Retau_ref), 'k')
    xlabel('$y^+$','Interpreter','latex')
    xlim([1 Retau_ref])
    ylabel('$U^+$','Interpreter','latex')
    xticks([1e-1 1e0 1e1 1e2])

    % Compute errors with respect to reference simulation

    aux1                       = ((q_rms_ref(:,1)*(Re/Retau_ref)).^2);
    aux2                       = ((q_rms(:,1)*(Re/Retau(n))).^2);
    amp1                       = aux1'*W*aux1;
    amp2                       = aux2'*W*aux2;
    f1                         = (aux1-aux2)'*W*(aux1-aux2);
    f2                         = (aux1)'*W*(aux1);
    err_u_rms(n)               = f1/f2;


    aux1                       = ((q_rms_ref(:,2)*(Re/Retau_ref)).^2);
    aux2                       = ((q_rms(:,2)*(Re/Retau(n))).^2);
    amp1                       = aux1'*W*aux1;
    amp2                       = aux2'*W*aux2;
    f1                         = (aux1-aux2)'*W*(aux1-aux2);
    f2                         = (aux1)'*W*(aux1);
    err_v_rms(n)               = f1/f2;

    aux1                       = (uv_rms_ref(:,1)*(Re/Retau_ref).^2);
    aux2                       = (uv_rms(:,1)*(Re/Retau(n)).^2);
    amp1                       = aux1'*W*aux1;
    amp2                       = aux2'*W*aux2;
    f1                         = (aux1-aux2)'*W*(aux1-aux2);
    f2                         = (aux1)'*W*(aux1);
    err_uv_rms(n)              = f1/f2;


    aux1                       = ((q_rms_ref(:,3)*(Re/Retau_ref)).^2);
    aux2                       = ((q_rms(:,3)*(Re/Retau(n))).^2);
    amp1                       = aux1'*W*aux1;
    amp2                       = aux2'*W*aux2;
    f1                         = (aux1-aux2)'*W*(aux1-aux2);
    f2                         = (aux1)'*W*(aux1);
    err_w_rms(n)               = f1/f2;
    
end
