%   This script plots a thresholded normal Lorentz force just
%   inside/outside any brain compartment surface (plots the surface field +
%   optionally coil geometry)
%
%   Copyright SNM/WAW 2017-2020

%%   Find the E-field just inside or just outside any model surface
par = +1;    %      par=-1 -> E-field just inside surface; par=+1 -> E-field just outside surface     
E = Einc + Eadd + par/(2)*normals.*repmat(c, 1, 3);    %   full field
%   Select surface/interface and compute field magnitude (tangential or normal or total)

%%  Identify the field
tissue_to_plot  = 'GM';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
condnumber      = find(strcmp(tissue, 'CSF'));
E               = E(Indicator==objectnumber, :);
B               = Binc(Indicator==objectnumber, :);
F               = cond(condnumber)* cross(E, B, 2);
Normals         = normals(Indicator==objectnumber, :);
Fnormal         = sum(F.*Normals, 2); % this is a projection onto normal vector (directed outside!)

%%   Identify the point cloud
i            = Indicator==objectnumber;
temp         = Fnormal;
[MAX, m]     = max(abs(temp))
Points       = Center(i, :);
th          = margin*MAX;
index       = find(abs(temp)>=th);
cloud       = Points(index, :);
position    = Points(m, :)
area        = Area(i);
TotalArea   = 1e6*sum(area(index))      % in mm^2

%%   Display the point cloud
figure;
S = load('sphere');
N = length(index);
n = length(S.P);
scale = 2.5;
for m = 1:N
    p = patch('vertices', scale*S.P+repmat(cloud(m, :), n, 1), 'faces', S.t);
    p.FaceColor = 'c';
    p.EdgeColor = 'none';
    p.FaceAlpha = 1.0;
end

%   Display the shell
p = patch('vertices', P, 'faces', t(i, :));
p.FaceColor = [1 0.75 0.65];
p.EdgeColor = 'none';
p.FaceAlpha = 1.0;
daspect([1 1 1]);
camlight; lighting phong;
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');

percentage  = num2str(margin*100);
maximum     = num2str(max(temp));
firstline   = strcat('Normal Lorentz force density is >', percentage, '%');
secondline  = strcat(' of the maximum value Fmax=', maximum, ' N/m^2');
title({firstline; secondline});

% Coil centerline graphics 
bemf1_graphics_coil_CAD(Coil.P, Coil.t, 1);
hold on;
plot3(pointsline(:, 1), pointsline(:, 2), pointsline(:, 3), '-m', 'lineWidth', 5);

% General
axis tight; view(30, 35), camzoom(1)
