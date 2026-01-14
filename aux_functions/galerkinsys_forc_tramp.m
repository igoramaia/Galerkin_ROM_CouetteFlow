function gsys = galerkinsys_forc_tramp(t,X,L,QQ,F,Re,targ1,targ2)
% Time integration of the forced system

XX      = X*X.'; XX = XX(:);
Xp      = 1/Re*L*X + QQ*XX +F*(0.5+0.5*tanh(0.1*(t-targ1)))-F*(0.5+0.5*tanh(0.1*(t-targ2))); Xp = full(Xp);
gsys    = Xp;
end