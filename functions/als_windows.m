function [copt,sopt] = als_windows(windows_i,num_comp)
% Function to perform MCR-ALS analysis.
% The pure concentrations are used as initial estimates.
% This function uses no-negativity constraints on both C and ST matrices.
% The number of iterations is set to 1000.
%
% OUTPUT:
%   copt: resolved concentration profiles.
%   sopt: resolved spectra profiles.

[sp,imp]=pure_modified(windows_i,num_comp,1);
[niter,r2,copt,sopt,sdopt,ropt,areaopt,rtopt]=als_command(windows_i,sp,1,1000);
end

