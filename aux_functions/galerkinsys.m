function gsys = galerkinsys(~,X,L,QQ,F,Re)
% Computes right-hand side of ODE

XX          = X*X.'; XX = XX(:);
Xp          = 1/Re*L*X + QQ*XX + F; Xp = full(Xp);  
gsys        = Xp;
end