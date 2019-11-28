function [md, fa, rd, ad, l1, l2, l3] = dti_parameters(D)
% Calculates mean diffusion (md), fractional anisotropy (fa)
% radial diffusivity (rd) and axial diffusivity (ad)
% from diffusion tensor D

%D = [d(1) d(2) d(3); d(2) d(4) d(5); d(3) d(5) d(6)];

[eigvec, eigval] = eigs(D);
eigval = diag(eigval);
%eigval = sort(abs(eigval), 'descend'); %XXX ABS???
eigval = sort(eigval, 'descend');

l1 = eigval(1);
l2 = eigval(2);
l3 = eigval(3);
 
%  lambda = eig(D); %should be equal to D_diag(L)
%  l1 = lambda(1);
%  l2 = lambda(2);
%  l3 = lambda(3);

md = (l1+l2+l3)/3;
rd = (l2+l3)/2;
ad = l1;

A = (l1 - md)^2 + (l2 - md)^2 + (l3 - md)^2;
B = l1^2 + l2^2 + l3^2;
fa = sqrt(3*A/(2*B));

%fa = sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2);


end