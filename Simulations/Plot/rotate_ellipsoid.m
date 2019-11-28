function [D_rot] = rotate_ellipsoid(Lvec, xrot, yrot, zrot)
%Rotates a diagonal tensor D_diag
%Lvec = diganocal elements of D_diag
%xrot, yrot, zrot = rotation angle (degrees)

L = Lvec;
D_diag  = diag(L);

%Angle of rotation
x = deg2rad(xrot); %pi/6;
y = deg2rad(yrot); %pi/3;
z = deg2rad(zrot); %pi/2;

%Rotation about x
Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
%Rotation about y
Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
%Rotation about z
Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];

%Rotation about z, then y, then x
R = Rx*Ry*Rz;

D2 = R*D_diag*R'; %D_diag rotated
a = 4;
D_rot = round(10^a*D2)/10^a; %Final D

end