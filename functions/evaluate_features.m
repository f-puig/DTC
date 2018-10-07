function [ignored]=evaluate_features(mat_file_path)

%The 'file.mat' contains the following objects:
%    --> The set of copt_wk
%    --> The set of sopt_wk
%    --> The ppm vector of the whole region: It must be named 'ppm'.

%Only Ck and Sk objects are named with the 'copt_wk' and 'sopt_wk'
%nomenclature.

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

%%% PLOT
ignored=[];
k = 0;
ppm_ini = 1;
limits = zeros(length(components_w_s),2);
num_comp=0;
for i =1:sum([copts{:,:}])
    num_comp2 = size(mat_file2.(sprintf('copt_w%d',i)),2);
    legend_label=[(num_comp+1):(num_comp+num_comp2)];
    leg_char={};
    for j=1:length(legend_label)
    leg_char{j}=sprintf('Component %d', legend_label(j));
    end
    subplot(1,2,1)
    plot(mat_file2.(sprintf('copt_w%d',i)))
    title(sprintf('C matrix for window %d',i),'fontweight','bold')
    xlabel('samples')
    legend(leg_char)
    num_comp = num_comp + num_comp2;
    subplot(1,2,2)
    ppm_w_size = size(mat_file2.(sprintf('sopt_w%d',i)),2);
    plot([ppm(ppm_ini:(ppm_ini+ppm_w_size-1))],mat_file2.(sprintf('sopt_w%d',i))')
    set(gca,'xdir','reverse')
    ppm_ini = ppm_ini + ppm_w_size;
    title(sprintf('S matrix for window %d',i),'fontweight','bold')
    xlabel('ppm')
    display('Do you want to ignore any of these components?')
    ans1 = input('Type the number. If not, type 0.');
    if ans1 ~= 0
        ignored = [ignored,ans1];
    end
end
end