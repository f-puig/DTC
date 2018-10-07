function [windows] = splitNMR(NMR,limits)
% This function splits a NMR dataset (with samples in rows, and intensity
% variables in columns) into a set of smaller spectral windows.
% The limits for these windows are defined by the 'limits' variable.
% The 'limits' variable is a matrix with two columns and as many rows as
% spectral windows to create.

% OUTPUT: A structure named 'windows' with as many fields as spectral
% windows.

name_files={};
windows=struct;
for i=1:size(limits,1)
    name_files{i} = sprintf('w_%d',i);
    windows.(name_files{i})=NMR(:,limits(i,1):limits(i,2));
end
end

