function inp = inprod_s(f1,f2,W)
% Compute inner product
u1      = f1(:,:,:,1); u1 = u1(:);
v1      = f1(:,:,:,2); v1 = v1(:);
w1      = f1(:,:,:,3); w1 = w1(:);

u2      = f2(:,:,:,1); u2 = u2(:);
v2      = f2(:,:,:,2); v2 = v2(:);
w2      = f2(:,:,:,3); w2 = w2(:);

inp = u1'*W*u2 + v1'*W*v2 + w1'*W*w2;
