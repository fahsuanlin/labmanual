%   This script computes and plots the Lorentz force (N/m^2) just inside/outside
%   any brain compartment surface (plots the surface field + optionally
%   coil geometry).
%
%   Copyright SNM/WAW 2017-2020

%%   Find the E-field just inside or just outside any model surface
par = +1;    %      par=-1 -> E-field just inside surface; par=+1 -> E-field just outside surface     
E = Einc + Eadd + par/(2)*normals.*repmat(c, 1, 3);    %   full field
%   Select surface/interface and compute field magnitude (tangential or normal or total)

tissue_to_plot  = 'GM';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
condnumber      = find(strcmp(tissue, 'CSF'));
E               = E(Indicator==objectnumber, :);
B               = Binc(Indicator==objectnumber, :);
F               = cond(condnumber)* cross(E, B, 2);
Normals         = normals(Indicator==objectnumber, :);
Fnormal         = sum(F.*Normals, 2); % this is a projection onto normal vector (directed outside!)
temp            = Normals.*repmat(Fnormal, 1, 3);
Ftangent        = F - temp;
Ftangent        = sqrt(dot(Ftangent, Ftangent, 2));
Ftotal          = sqrt(dot(F, F, 2));

f.MAXFtotal     = max(Ftotal);
f.MAXFnormal    = max(abs(Fnormal));
f.MAXFtangent   = max(abs(Ftangent));
f

%%   Graphics
figure;
temp = Ftotal;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
if par == +1; string = ' just outside:'; end
if par == -1; string = ' just inside:'; end
title(strcat('Solution: Lorentz force-field (total, normal, or tang.) in N/m^2: ', string, tissue{objectnumber}));

% Coil centerline graphics 
bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
hold on;
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

% General
axis tight; view(30, 35), camzoom(1);