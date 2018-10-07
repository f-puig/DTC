function [supersopt,supercopt,names_sopts]=supersopt_supercopt_ppm(mat_file_path,ignored)
%This function is a variant of the first one (i got some errors sometimes
%when appending the S features but i don't remember why)

%With this function, an augmented matrix of Ck (supercopt) and an augmented
%cell of Sk are constructed.
%Each Sk contains two vectors:
%    --> The vector of ppm
%    --> The vector of intensities.

%The 'file.mat' contains the following objects:
%    --> The set of copt_wk
%    --> The set of sopt_wk
%    --> The ppm vector of the whole region: It must be named 'ppm'.

%Only Ck and Sk objects are named with the 'copt_wk' and 'sopt_wk'
%nomenclature.
if nargin == 1
    ignored= [];
end
% 1. Load .Mat file
mat_file = matfile(mat_file_path);
mat_file2=load(mat_file_path);
names = who(mat_file);

ppm=mat_file2.ppm;
length_ppm=length(ppm);

%2.Check dimensions
copts = strfind(names,'copt_w');
sopts = strfind(names,'sopt_w');
position_copts=[];
position_sopts=[];
names_sopts={};
names_sopts{1,1}='a';
l = 0;
for i=1:size(sopts,1)
    if sopts{i,1}==1
        l = l + 1;
        names_sopts{1,l} = names{i};
    end
end
size_copts=[];
ppm_sopts=[];

components_w_c=[];
components_w_s=[];
j=0;
for i=1:length(copts)
    if ~isempty(copts{i})
        j=j+1;
        position_copts(j) = i;
        size_copts(j)=size(mat_file2.(names{i}),1);
        components_w_c(j)=size(mat_file2.(names{i}),2);
    end
end

j=0;
for i=1:length(sopts)
    if ~isempty(sopts{i})
        j=j+1;
        position_sopts(j) = i;
        ppm_sopts(j)=size(mat_file2.(names{i}),2);
        components_w_s(j)=size(mat_file2.(names{i}),1);
    end
end

if length(components_w_s) ~= length(components_w_c)
    display('There is a different number of Ck and Sk');
    return
end

if sum(components_w_s) ~= sum(components_w_c)
    display('The total number of components for the C and S submatrices do not agree');
    for i=1:length(components_w_s)
        sk_size=size(mat_file2.(sprintf('sopt_w%d',i)),1);
        ck_size=size(mat_file2.(sprintf('copt_w%d',i)),2);
        if sk_size ~= ck_size
            display(sprintf('Window %d does not have the same size for Ck and Sk.', i));
        end
    end
    return
end

if length_ppm ~= sum(ppm_sopts)
    display('The augmented supercopt will not be equal to the number of ppm');
    return
end

if length(unique(size_copts)) ~= 1
   display('Not all Ck submatrices contain the same number of samples (rows)');
   return
end

% 3. Create supercopt
supercopt=[];
for i=1:length(components_w_c)
supercopt=[supercopt,mat_file2.(sprintf('copt_w%d',i))];
end
supercopt(:,ignored)=[];
% 4. create supersopt
supersopt = {};
k = 0;
m = 0;
ppm_ini = 1;
limits = zeros(length(components_w_s),2);
supersopt{1,1}=0;

for i = 1:length(components_w_s)
    %i
    a = mat_file2.(sprintf('sopt_w%d',i));
    ppm_w_size = size(mat_file2.(sprintf('sopt_w%d',i)),2);
    limits(i,1)=ppm_ini;
    limits(i,2)=ppm_ini+ppm_w_size-1;
    for j = 1:size(mat_file2.(sprintf('sopt_w%d',i)),1)
        m = m+1;
        if ~ismember(m,ignored)
            k = k+1;
            supersopt{1,k} = [ppm(ppm_ini:(ppm_ini+ppm_w_size-1));a(j,:)];
        end
    end
    ppm_ini = ppm_ini + ppm_w_size;
end

end