%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate ROM for Couette flow using controllability and Stokes modes                   
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


%% Domain parameters (adapt to avoid aliasing)

Nx      = 10;                   % Number of grid points in the streamwise direction
Nz      = 10;                   % Number of grid points in the spanwise direction
Ny      = 65;                   % Number of grid points in the wall-normal direction
Lx      = 2.0*pi;               % Streamwise domain size
Lz      = 1.0*pi;               % Spanwise domain size


%% Grid generation

% Wall-normal grid (Chebyshev discretisation)
[y,A]   = chebdif(Ny,2);        % Chebyshev grid in the wall-normal direction
[~,W]   = clencurt(Ny-1);       % Clenshaw-Curtis quadrature weights
Wcol    = repmat(W.',Nx*Nz,1);
W       = diag(W);
Dy      = squeeze(A(:,:,1));    % Differentiation matrix
Dy2     = squeeze(A(:,:,2));    % Second-order differentiation


% x-z grid (homegeneous directions, Fourier discretisation)
[x,Dx]  = fourdif(Nx,1);
[~,Dx2] = fourdif(Nx,2);
[z,Dz]  = fourdif(Nz,1);
[~,Dz2] = fourdif(Nz,2);

% Scaling x-z grid
facx    = Lx/(2*pi); x = x*facx; Dx = Dx/facx; Dx2 = Dx2/facx^2;
facz    = Lz/(2*pi); z = z*facz; Dz = Dz/facz; Dz2 = Dz2/facz^2;
dx      = x(2)-x(1);
dz      = z(2)-z(1);

%% ROM parameters

Ulam    = y;                % Laminar flow solution (used to compute basis)
Re_bas  = 100;              % Reynolds number for basis computation
nalphas = 2;                % Number of non-zero streamwise wavenumbers
ngammas = 2;                % Number of non-zero spanwise wavenumbers
nmodesy = 36;               % Number of controllability modes or Stokes modes per wavenumber pair
alphas  = (0:nalphas)*2*pi/Lx;
gammas  = (-ngammas:1:ngammas)*2*pi/Lz;
gref    = 2*pi/Lz;
[a,g]   = meshgrid(alphas,gammas);
[X,Y,Z] = meshgrid(x,y,z);
Wf      = 0.5/(Nx*Nz)*spdiags(Wcol,0,length(Wcol),length(Wcol)); % Global quadrature weight matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create modal basis                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise basis
phiaux  = zeros(size(X,1),size(X,2),size(X,3),3,size(a,1),size(a,2),nmodesy,2);
avec    = zeros(size(a,1),size(a,2),nmodesy,2);
gvec    = avec; bvec = avec;
Phases  = zeros(size(a,1),size(a,2),nmodesy,2);

for i=1:size(a,1)
    for j=1:size(a,2)
    	if((a(i,j)==0)&&(g(i,j)==0))                                 % if streamwise and spanwise wavenumbers are null, use Stokes modes
            q = modalbasis_s(a(i,j),g(i,j),Ny,nmodesy);              % Stokes modes
    	else
            q = modalbasis_cntrl(a(i,j),g(i,j),Ulam,Re_bas,nmodesy); % Controllability modes
    	end
        avec(i,j,:,:) = a(i,j);
        gvec(i,j,:,:) = g(i,j);
        for k=1:nmodesy
            bvec(i,j,k,:) = k;
            fu = q(1:Ny,k);
            fv = q(Ny+1:2*Ny,k);
            fw = q(2*Ny+1:3*Ny,k);
            poscentre = round((Ny+1)/2);
            if(abs(fu(poscentre))>1e-5)
                phasemode = angle(fu(poscentre));
            elseif(abs(fw(poscentre))>1e-5)
                phasemode = angle(fw(poscentre));
            else
                phasemode = angle(fv(poscentre));
            end

            fu = fu*exp(-1i*phasemode+1i*pi/2);
            fv = fv*exp(-1i*phasemode+1i*pi/2);
            fw = fw*exp(-1i*phasemode+1i*pi/2);

            phiaux(:,:,:,1,i,j,k,1) = real(exp(1i*a(i,j)*X).*fu.*exp(1i*g(i,j)*Z));
            phiaux(:,:,:,2,i,j,k,1) = real(exp(1i*a(i,j)*X).*fv.*exp(1i*g(i,j)*Z));
            phiaux(:,:,:,3,i,j,k,1) = real(exp(1i*a(i,j)*X).*fw.*exp(1i*g(i,j)*Z));
            Phases(i,j,k,1)         = 0.0;

            phiaux(:,:,:,1,i,j,k,2) = imag(exp(1i*a(i,j)*X).*fu.*exp(1i*g(i,j)*Z));
            phiaux(:,:,:,2,i,j,k,2) = imag(exp(1i*a(i,j)*X).*fv.*exp(1i*g(i,j)*Z));
            phiaux(:,:,:,3,i,j,k,2) = imag(exp(1i*a(i,j)*X).*fw.*exp(1i*g(i,j)*Z));
            Phases(i,j,k,2)         = pi/2;
        end
    end
end

% Organise basis
phi         = phiaux(:,:,:,:,:);
nmodes      = size(phi,5);
avec        = avec(:);
gvec        = gvec(:);
kvec        = [avec gvec];      % Vector with wavenumber pairs contained in the basis
bvec        = bvec(:);          % Vector with controllability order of each mode
Phases      = Phases(:);

% throw away (0,-gamma) modes (Cavalieri & Nogueira, PRF 2022);
pos         = find(~((kvec(:,1)==0)&kvec(:,2)<0));
phi         = phi(:,:,:,:,pos);
kvec        = kvec(pos,:);
bvec        = bvec(pos,:);
Phases      = Phases(pos);
nmodes      = length(pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check orthogonality                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = zeros(nmodes,nmodes);
for i=1:nmodes
    disp(['Inner products in basis, ' num2str(i) '/' num2str(nmodes)])
    for j=1:nmodes
        M(i,j) = inprod_s(squeeze(phi(:,:,:,:,i)),squeeze(phi(:,:,:,:,j)),Wf);
    end
end
disp(['Maximum non-orthogonality: ' num2str(max(max(M-diag(diag(M)))))]);

% Remove modes with zero norm
normmode            = diag(M);
pos                 = find((normmode>1e-4));
phi                 = squeeze(phi(:,:,:,:,pos));
M                   = M(pos,pos);
nmodes              = length(pos);
kvec                = kvec(pos,:);
bvec                = bvec(pos,:);
Phases              = Phases(pos);

% Normalise modes
for i=1:nmodes
   phi(:,:,:,:,i)   = phi(:,:,:,:,i)/sqrt(M(i,i));
   M(i,i)           = M(i,i)/M(i,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute ROM operators (see mathematical details in papers)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = (1/Re)L_{ij}a_{j} + L2_{ij}a_{j} +Q_{ikj}a_{j}a_{k}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phix                = phi;
phixx               = phi;
phiy                = phi;
phiyy               = phi;
phiz                = phi;
phizz               = phi;

% Initialise base flow 
phi0                = zeros(size(phi,1),size(phi,2),size(phi,3),3);
phi0(:,:,:,1)       = Y; % Laminar solution

% Differential operators
for i=1:nmodes
    disp(['Differentiation, ' num2str(i) '/' num2str(nmodes)])

    for j=1:Ny
        for k=1:Nz
            phix(j,:,k,:,i) = Dx*squeeze(phi(j,:,k,:,i));
            phixx(j,:,k,:,i) = Dx2*squeeze(phi(j,:,k,:,i));
        end
    end
    for j=1:Nx
        for k=1:Nz
            phiy(:,j,k,:,i) = Dy*squeeze(phi(:,j,k,:,i));
            phiyy(:,j,k,:,i) = Dy2*squeeze(phi(:,j,k,:,i));
        end
    end
    for j=1:Ny
        for k=1:Nx
            phiz(j,k,:,:,i) = Dz*squeeze(phi(j,k,:,:,i));
            phizz(j,k,:,:,i) = Dz2*squeeze(phi(j,k,:,:,i));
        end
    end
end

% Base flow differential operators
phi0x       = zeros(Ny, Nx, Nz,3);
phi0y       = zeros(Ny, Nx, Nz,3);
phi0z       = zeros(Ny, Nx, Nz,3);

for j=1:Ny
    for k=1:Nz
        phi0x(j,:,k,:) = Dx*squeeze(phi0(j,:,k,:));
    end
end
for j=1:Nx
    for k=1:Nz
        phi0y(:,j,k,:) = Dy*squeeze(phi0(:,j,k,:));
    end
end
for j=1:Ny
    for k=1:Nx
        phi0z(j,k,:,:) = Dz*squeeze(phi0(j,k,:,:));
    end
end

% Laplacian
Lapl                   = phixx+phiyy+phizz;

% Check if modes are divergence-free
divU                    = zeros(nmodes,1);
for i=1:nmodes
divU(i)                 = max(max(max(abs(squeeze(phix(:,:,:,1,i)) + squeeze(phiy(:,:,:,2,i)) + + squeeze(phiz(:,:,:,3,i))))));
end
disp(['Max. divergence = ' num2str(max(max(max(max(abs(divU)))))) ]);
pos                     = find(divU>0.01);
disp('Wavenumbers with non-zero divergence: ');
kvec(pos,:) 

% Operators in Galerkin projection
L                       = zeros(nmodes,nmodes);
L2                      = zeros(nmodes,nmodes);

% Linear diffusion term
for i=1:nmodes
    disp(['Linear term, ' num2str(i) '/' num2str(nmodes)])
    for j=1:nmodes
        L(i,j)          = inprod_s(phi(:,:,:,:,i),Lapl(:,:,:,:,j),Wf);
    end
end

% Linear convective term (energy extraction from the base flow)
for i=1:nmodes
    disp(['Linear term 2, ' num2str(i) '/' num2str(nmodes)])
    f1              = squeeze(phi(:,:,:,:,i));

    for j=1:nmodes
        f2          = zeros(size(f1));
        f2(:,:,:,1) = phi(:,:,:,1,j).*phi0x(:,:,:,1) + phi(:,:,:,2,j).*phi0y(:,:,:,1) + phi(:,:,:,3,j).*phi0z(:,:,:,1);
        f2(:,:,:,2) = phi(:,:,:,1,j).*phi0x(:,:,:,2) + phi(:,:,:,2,j).*phi0y(:,:,:,2) + phi(:,:,:,3,j).*phi0z(:,:,:,2);
        f2(:,:,:,3) = phi(:,:,:,1,j).*phi0x(:,:,:,3) + phi(:,:,:,2,j).*phi0y(:,:,:,3) + phi(:,:,:,3,j).*phi0z(:,:,:,3);

        f2(:,:,:,1) = f2(:,:,:,1) + phix(:,:,:,1,j).*phi0(:,:,:,1) + phiy(:,:,:,1,j).*phi0(:,:,:,2) + phiz(:,:,:,1,j).*phi0(:,:,:,3);
        f2(:,:,:,2) = f2(:,:,:,2) + phix(:,:,:,2,j).*phi0(:,:,:,1) + phiy(:,:,:,2,j).*phi0(:,:,:,2) + phiz(:,:,:,2,j).*phi0(:,:,:,3);
        f2(:,:,:,3) = f2(:,:,:,3) + phix(:,:,:,3,j).*phi0(:,:,:,1) + phiy(:,:,:,3,j).*phi0(:,:,:,2) + phiz(:,:,:,3,j).*phi0(:,:,:,3);
        L2(i,j) = -inprod_s(f1,f2,Wf);
    end
end

% Quadratic convective term term
Q      = zeros(nmodes,nmodes,nmodes);
parpool(24); % Parallel computing
parfor i=1:nmodes
    tic;
    disp(['Quadratic term, ' num2str(i) '/' num2str(nmodes)])
    f1                  = squeeze(phi(:,:,:,:,i));
    for j=1:nmodes
        for k=1:nmodes
            f2          = zeros(size(f1));
            f2(:,:,:,1) = phi(:,:,:,1,j).*phix(:,:,:,1,k) + phi(:,:,:,2,j).*phiy(:,:,:,1,k) + phi(:,:,:,3,j).*phiz(:,:,:,1,k);
            f2(:,:,:,2) = phi(:,:,:,1,j).*phix(:,:,:,2,k) + phi(:,:,:,2,j).*phiy(:,:,:,2,k) + phi(:,:,:,3,j).*phiz(:,:,:,2,k);
            f2(:,:,:,3) = phi(:,:,:,1,j).*phix(:,:,:,3,k) + phi(:,:,:,2,j).*phiy(:,:,:,3,k) + phi(:,:,:,3,j).*phiz(:,:,:,3,k);
            Q(i,j,k) = -inprod_s(f1,f2,Wf);
        end
    end
    toc
end

for i=1:nmodes
    for j=1:nmodes
        for k=1:(j-1)
            Q(i,j,k) = Q(i,j,k)+Q(i,k,j);
            Q(i,k,j) = 0;
        end
    end
end

% Flatten quadratic term for better computation later (see the runROM.m and galerkinsys.m scripts)
QQ          = zeros(nmodes,nmodes^2);
for i=1:nmodes
    Qaux    = squeeze(Q(i,:,:)); Qaux = Qaux(:).';
    QQ(i,:) = (Qaux);
end
QQ          = sparse((abs(QQ)>1e-8).*QQ);

% Write output file
outfile = ['galerkin_model_cntrl_' int2str(nalphas) '_' int2str(nmodesy) '_' int2str(ngammas) '_Lx_',num2str(Lx/pi),'_Lz_',num2str(Lz/pi),'.mat'];
save(outfile,'L','L2','Q','QQ','phi','phi0','phix','phiy','phiz','X','Y','Z','Wf','Phases','avec','bvec','gvec','kvec','alphas','gammas','Nx','Ny','Nz','-v7.3');

