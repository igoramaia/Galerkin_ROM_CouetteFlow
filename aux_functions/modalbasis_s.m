function q = modalbasis_s(alpha,beta,Ny,nmodes)
% This function computes Stokes modes for the case of zero streamwise and
% spanwise wavenumbers


[y,A]           = chebdif(Ny,2);
[~,W]           = clencurt(Ny-1);

Dy              = squeeze(A(:,:,1));
Dy2             = squeeze(A(:,:,2));


% Stokes operator laplacian(u) - grad(p) = lambda*p
II              = eye(size(Dy)); ZZ = zeros(size(II));
k               = sqrt(alpha^2+beta^2);
Lapl            = Dy2-k^2*II;

L               = [Lapl ZZ ZZ -1i*alpha*II;
                   ZZ Lapl ZZ -Dy;
                   ZZ ZZ Lapl -1i*beta*II;
                   1i*alpha*II Dy 1i*beta*II ZZ];
F               = [II ZZ ZZ ZZ;
                   ZZ II ZZ ZZ;
                   ZZ ZZ II ZZ;
                   ZZ ZZ ZZ ZZ];

% BCs
L(1,:)          = 0; L(1,1) = 1.; F(1,:) = 0; %u(pm1)=0
L(Ny,:)         = 0; L(Ny,Ny) = 1.; F(Ny,:) = 0;
L(Ny+1,:)       = 0; L(Ny+1,Ny+1) = 1.; F(Ny+1,:) = 0; %v(pm1)=0
L(2*Ny,:)       = 0; L(2*Ny,2*Ny) = 1.; F(2*Ny,:) = 0;
L(2*Ny+1,:)     = 0; L(2*Ny+1,2*Ny+1) = 1.; F(2*Ny+1,:) = 0; %v(pm1)=0
L(3*Ny,:)       = 0; L(3*Ny,3*Ny) = 1.; F(3*Ny,:) = 0;

[V,lambda]      = eig(L,F); %Stokes modes
lambda          = diag(lambda);
pos             = find(isfinite(lambda));
lambda          = lambda(pos); V=V(:,pos);
pos             = find(abs(lambda)<1e6);
lambda          = lambda(pos); V=V(:,pos);

[~,pos]         = sort(abs(lambda),'ascend');
lambda          = lambda(pos);
V               = V(:,pos);

q               = V(:,1:nmodes);
if(max(lambda(1:nmodes)<0)>0) disp(['Error in Stokes modes, alpha,beta = ' num2str(alpha) num2str(beta)]); end
