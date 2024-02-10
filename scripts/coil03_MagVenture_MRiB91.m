%   This script creates the mesh (both CAD surface mesh and a computational
%   wire grid) for a figure-eight elliptical coil with 1 A of total current
%   The output is saved in the binary file coil.mat and includes:
%   strcoil.Pwire(:, 3) - set of nodes for all wires 
%   strcoil.Ewire(:, 2) - set of edges or current dipoles for all wires
%   (current flows from the first edge node to the second edge node)
%   strcoil.Swire{:, 1} - current strength weight for every elementary
%   dipole asssuring that the total conductor current through any
%   cross-section is 1 A.
%
%   Copyright SNM 2018-2020

clear all; %#ok<CLALL>
if ~isunix
    s = pwd; addpath(strcat(s(1:end-5), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-5), '/Engine'));
end

%   The figure-eight coil includes elliptical turns/loops. The coil axis is
%   the z-axis. We construct one part of the coil first.
%   When crossing the xz-plane, the i   ntersection points for the loop
%   centerlines are
x0  = 1e-3*[20.2500   25.7500   31.2500   36.7500   ...
            20.2500   25.7500   31.2500   36.7500];  
%   When crossing the yz-plane, the intersection points for the loop
%   centerlines are
y0  = 1e-3*[28.7500   34.2500   39.7500   45.2500  ...
            28.7500   34.2500   39.7500   45.2500];
%   When crossing the xz- or yz-planes, the intersection points for the
%   loop centerlines are (this is a z-offset)
z0  = 1e-3*[-4.6000   -4.6000   -4.6000   -4.6000  ...
            4.6000    4.6000    4.6000    4.6000];

%   Other parameters
a    = 3.50e-3;     %   z-side, m  (for a rectangle cross-section)
b    = 2.20e-3;     %   x-side, m  (for a rectangle cross-section)
M    = 32;          %   number of cross-section subdivisions 
N    = 128;          %   number of perimeter subdivisions
flag = 2;           %   rectangular cross-section    
sk   = 1;           %   surface current distribution (skin layer)
        
%   Construct the coil mesh for one part
[Pwire, Ewire, Swire, P, t, tind] = meshcoil(x0, y0, z0, M, N, a, b, flag, sk);
  
%   Construct the entire coil as a combination of two parts
%   Separate two coils along the x-axis
offset = -39.5e-3;
Pwire1 = Pwire; Pwire1(:, 1) = Pwire1(:, 1) - offset;
Pwire2 = Pwire; Pwire2(:, 1) = Pwire2(:, 1) + offset;
Ewire1 = Ewire;
Ewire2 = Ewire+size(Pwire, 1);
 
strcoil.Pwire = [Pwire1; Pwire2]; 
strcoil.Pwire(:, 3) = strcoil.Pwire(:, 3) - min(strcoil.Pwire(:, 3));
strcoil.Ewire = [Ewire1; Ewire2];
strcoil.Swire = [+Swire; -Swire]; %  swap current direction for the second part

P1 = P; P1(:, 1) = P(:, 1) - offset;
P2 = P; P2(:, 1) = P(:, 1) + offset;
t1 = t;
t2 = t+size(P, 1);
P = [P1; P2];
P(:, 3) = P(:, 3) - min(P(:, 3)); 
t    = [t1; t2];
tind = [tind; tind+max(tind)];

%   Display CAD and wire models
bemf1_graphics_coil_CAD(P, t, 0);
%bemf1_graphics_coil_wire(strcoil, [0 1 0]); 
view(-4, 24);
 
save('coil', 'strcoil');
save('coilCAD', 'P', 't', 'tind');  %   optional, slow
