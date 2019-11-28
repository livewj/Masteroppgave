function plot_ellipsoid(L, xrot, yrot, zrot)
%Function that takes the diffusion tensor eigenvalues and angle of
%rotation (degrees) and
%illustrates the corresponding diffusion ellipsoid

D_diag = diag(L); %Diffusion tensor

%center
xc=0; yc=0; zc=0;

%semi-axis length
xr=L(1);
yr=L(2);
zr=L(3);

n=50; %number of grids

%generate ellipsoid 
[x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,n); 


%Rotate the ellipsoid
fontsize = 22; % 14 orig

figure, clf 
s = surf(x, y, z);
rotate(s, [1 0 0], xrot); %Rotation about x
rotate(s, [0 1 0], yrot); %Rotation about y
rotate(s, [0 0 1], zrot); %Rotation about z
axis equal
set(s,'EdgeColor','none'); % xxx
set(gca, 'fontsize', fontsize);
axis equal
camlight left
lighting gouraud
shading interp
hold on
plot3([-3 3], [0 0], [0 0]), hold on, plot3([0 0], [-3 3], [0 0]), hold on, plot3([0 0], [0 0], [-3 3]);
xlabel('x', 'fontsize', fontsize)
ylabel('y', 'fontsize', fontsize)
zlabel('z', 'fontsize', fontsize)
title('Diffusion ellipsoid')
w = 7;
h = 5;
set(gcf, 'PaperSize', [w h]);
set(gcf, 'PaperPosition', [0 0 w h]);

end