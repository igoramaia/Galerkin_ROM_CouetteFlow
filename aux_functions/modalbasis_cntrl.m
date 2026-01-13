function q = modalbasis_cntrl(alpha,beta,U,Re,nmodes)
% This function computes controllability modes from the linearised
% Navier-Stokes equations, doubly homegeneous in the streamwise and spanwise directions,
% forced with with noise, as described by Jovanovic & Bamieh (J. Fluid Mech., 2005)

Ny              = length(U);
[y,A]           = chebdif(Ny,2);
[~,W]           = clencurt(Ny-1);
Dy              = squeeze(A(:,:,1));
Dy2             = squeeze(A(:,:,2));

% Base flow
dU              = Dy*U; ddU = Dy2*U;

% Dirichlet BCs for everything
Dy              = Dy(2:Ny-1,2:Ny-1);
Dy2             = Dy2(2:Ny-1,2:Ny-1);

y               = y(2:Ny-1);
U               = U(2:Ny-1); 
dU              = dU(2:Ny-1);
ddU             = ddU(2:Ny-1);

W               = W(2:Ny-1); Wh = sqrt(W);
W               = diag(W); 
ZZ              = zeros(size(W));

Wchol           = [diag(Wh) ZZ ZZ;
                ZZ diag(Wh) ZZ;
                ZZ ZZ diag(Wh)];
iWchol          =[diag(1./Wh) ZZ ZZ;
                ZZ diag(1./Wh) ZZ;
                ZZ ZZ diag(1./Wh)];

% Clamped BCs for D4
[~,Dy4]         = cheb4c(Ny);

% Jovanovic & Bamieh 2005
II              = eye(size(Dy)); ZZ = zeros(size(II));
k               = sqrt(alpha^2+beta^2);
Lapl            = Dy2-k^2*II;
Lapl2           = Dy4-2*k^2*Dy2+k^4*II;

iLapl           = inv(Lapl);

A11             = -1i*alpha*iLapl*diag(U)*Lapl + 1i*alpha*iLapl*diag(ddU) + 1/Re*iLapl*Lapl2;
A12             = ZZ;
A21             = -1i*beta*diag(dU);
A22             = -1i*alpha*diag(U) + 1/Re*Lapl;

A               = [A11 A12;
                   A21 A22];

Ba              = [iLapl ZZ;
                   ZZ II];
Bb              = [-1i*alpha*Dy -k^2*II -1i*beta*Dy;
                   1i*beta*II ZZ -1i*alpha*II];

B               = Ba*Bb;
Bxtil           = [-1i*alpha*Dy;
                  1i*beta*II];
Bytil           = [-k^2*II;
                   ZZ];
Bztil           = [-1i*beta*Dy;
                   -1i*alpha*II];
Bx              = Ba*Bxtil;
By              = Ba*Bytil;
Bz              = Ba*Bztil;


A11d            = 1i*alpha*diag(U) -1i*alpha*iLapl*diag(ddU) + 1/Re*iLapl*Lapl2;
A12d            = ZZ;
A21d            = -1i*beta*iLapl*diag(dU);
A22d            = 1i*alpha*diag(U) + 1/Re*Lapl;

Ad              = [A11d A21d;
                   A12d A22d];

Cu              = 1/k^2*[1i*alpha*Dy -1i*beta*II];
Cv              = 1/k^2*[k^2*II ZZ;];
Cw              = 1/k^2*[1i*beta*Dy 1i*alpha*II];

C               = [Cu; Cv; Cw];

BxA             = Cu;
ByA             = Cv;
BzA             = Cw;

CuA             = Bx;
CvA             = By;
CwA             = Bz;

BBA             = [II ZZ;...
                   ZZ II];
CAC             = BBA;

Xr              = lyap(A,Ad,BBA); % Lyapunov equation, controllability Gramian

% Eigenvalue decomposition
[Vx,lx]         = eig(Xr);
lx              = diag(lx);
[~,order]       = sort(-real(lx));
lx              = lx(order);
Vx              = Vx(:,order);
np              = Ny-2;
q               = zeros(3*Ny,nmodes);

if (alpha==0 && beta==0)
    % treat special case alpha=beta=0 and assign Stokes modes
    % W(y)=cos(beta*y) and U(y)=cos(beta(y)).
    % I am avoiding the mode U(y)=cos(beta*y) to avoid mass-flux
    % fluctuation. Should probably avoid other modes.
    for i=1:2:nmodes
        bbb                     = pi/2;
        if(mod((round(i+1)/2),2)==1)
            q(1:np,i)           = cos((i+1)/2*bbb*y);
            q(2*np+1:3*np,i+1)  = cos((i+1)/2*bbb*y);
        else
            q(1:np,i)           = sin((i+1)/2*bbb*y);
            q(2*np+1:3*np,i+1)  = sin((i+1)/2*bbb*y);
        end

    end
else
    Ux                          = C*Vx;
    q                           = Ux(:,1:nmodes);
end
% Put zeros at y=\pm 1
qaux                            = zeros(3*Ny,nmodes);
qaux(2:Ny-1,:)                  = q(1:np,:);
qaux(Ny+2:2*Ny-1,:)             = q(np+1:2*np,:);
qaux(2*Ny+2:3*Ny-1,:)           = q(2*np+1:3*np,:);
q                               = qaux;
