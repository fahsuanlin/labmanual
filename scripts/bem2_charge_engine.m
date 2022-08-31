%   This script computes the induced surface charge density for an
%   inhomogeneous multi-tissue object in an external time-varying magnetic
%   field due to a coil via a BEM/FMM iterative solution with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Parameters of the iterative solution
iter         = 14;              %    Maximum possible number of iterations in the solution 
relres       = 1e-12;           %    Minimum acceptable relative residual 
weight       = 1/2;             %    Weight of the charge conservation law to be added (empirically found)

%%  Right-hand side b of the matrix equation Zc = b. Compute pointwise
%   Surface charge density is normalized by eps0: real charge density is eps0*c
tic
EincP    = bemf3_inc_field_electric(strcoil, P, dIdt, mu0);             %   Incident coil field
Einc     = 1/3*(EincP(t(:, 1), :) + EincP(t(:, 2), :) + EincP(t(:, 3), :));
b        = 2*(contrast.*sum(normals.*Einc, 2));                         %   Right-hand side of the matrix equation
IncFieldTime = toc

%%  GMRES iterative solution (native MATLAB GMRES is used)
h           = waitbar(0.5, 'Please wait - Running MATLAB GMRES');  
%   MATVEC is the user-defined function of c equal to the left-hand side of the matrix equation LHS(c) = b
MATVEC = @(c) bemf4_surface_field_lhs(c, Center, Area, contrast, normals, weight, EC);     
[c, flag, rres, its, resvec] = gmres(MATVEC, b, [], relres, iter, [], [], 8*b); 
close(h);

%%  Plot convergence history
figure; 
semilogy(resvec/resvec(1), '-o'); grid on;
title('Relative residual of the iterative solution');
xlabel('Iteration number');
ylabel('Relative residual');

%%  Check charge conservation law (optional)
conservation_law_error = sum(c.*Area)/sum(abs(c).*Area)

%%  Check the residual of the integral equation
solution_error = resvec(end)/resvec(1)

%%   Topological low-pass solution filtering (repeat if necessary)
c = (c.*Area + sum(c(tneighbor).*Area(tneighbor), 2))./(Area + sum(Area(tneighbor), 2));

%%  Save solution data (surface charge density, principal value of surface field)
tic
save('output_charge_solution', 'c', 'Einc', 'resvec', 'conservation_law_error', 'solution_error');
save_charge_solution_time = toc

%%   Find and save surface fields
%   (i)     total normal E-field just inside/outside any model surface; 
%   (ii)    secondary continuous E-field contribution for any model surface; 
%   (iii)   secondary continuous electric potential for any model surface; 
Eninside     = condout./(condin-condout).*c;    %   since c is normalized by eps0
Enoutside    = condin./(condin-condout).*c;     %   since c is normalized by eps0

tic
h    = waitbar(0.5, 'Please wait - computing accurate surface electric field'); 
[Pot, Eadd] = bemf4_surface_field_electric_subdiv(c, P, t, Area, 'barycentric', 3);
close(h);
Esurface_field_time = toc

tic
h    = waitbar(0.5, 'Please wait - computing coil magnetic field'); 
Binc = bemf3_inc_field_magnetic(strcoil, Center, I0, mu0);
close(h);
BincFieldTime = toc

tic
save('output_field_solution.mat', 'Eninside', 'Enoutside', 'Eadd', 'Pot', 'Binc');
save_E_solution_time = toc