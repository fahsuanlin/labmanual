
%   This script computes and plots the electric field just inside/outside
%   any brain compartment surface (plots the surface field + optionally
%   coil geometry).
%
%   Copyright SNM/WAW 2017-2020

%%   Find the E-field just inside or just outside any model surface
par = -1;    %      par=-1 -> E-field just inside surface; par=+1 -> E-field just outside surface     
E = Einc + Eadd + par/(2)*normals.*repmat(c, 1, 3);    %   full field
%   Select surface/interface and compute field magnitude (tangential or normal or total)

tissue_to_plot = 'GM_LH';

objectnumber= find(strcmp(tissue, tissue_to_plot));
E          = E(Indicator==objectnumber, :);
Normals     = normals(Indicator==objectnumber, :);
Enormal     = sum(E.*Normals, 2); % this is a projection onto normal vector (directed outside!)
temp        = Normals.*repmat(Enormal, 1, 3);
Etangent    = E - temp;
Etangent    = sqrt(dot(Etangent, Etangent, 2));
Etotal      = sqrt(dot(E, E, 2));
e.MAXEtotal     = max(Etotal);
e.MAXEnormal    = max(abs(Enormal));
e.MAXEtangent   = max(abs(Etangent));
%e

%%   Graphics
temp = Etotal;
figure;
bemf2_graphics_surf_field(P, t, etc_threshold(temp,0.999), Indicator, objectnumber);
if par == +1; string = ' just outside:'; end
if par == -1; string = ' just inside:'; end
title(strcat('Solution: E-field (total, normal, or tang.) in V/m: ', string, tissue{objectnumber}));

% Coil centerline graphics 
bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
hold on;
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

% General
axis tight; view(30, 35), camzoom(1);

brighten(0.4);

save field_e_fhlin.mat P t Etotal Indicator objectnumber