function [ssel, ini_est]=nmr_ssel_iniest(sopt, copt, grouped_copt_indeces, ppm)

%These three lines are added because each sopt was transposed in a previous
%verion. It might be not required.
for i=1:length(sopt)
    sopt{1,i}=sopt{1,i}';
end
%

ini_est=zeros(size(grouped_copt_indeces,2),length(ppm));
ssel=zeros(size(grouped_copt_indeces,2),length(ppm));

for i=1:size(grouped_copt_indeces,2)
    for j=1:size(grouped_copt_indeces{1,i},2)
        beg=find(ppm == sopt{grouped_copt_indeces{1,i}(j)}(1,1));
        ends=find(ppm == sopt{grouped_copt_indeces{1,i}(j)}(end,1));
        ini_est(i,beg:ends)=sopt{grouped_copt_indeces{1,i}(j)}(:,2)*mean(copt(:,grouped_copt_indeces{1,i}(j)));
    end
end

for i=1:size(grouped_copt_indeces,2)
    for j=1:size(grouped_copt_indeces{1,i},2)
        ssel_beg=find(ppm == sopt{grouped_copt_indeces{1,i}(j)}(1,1));
        ssel_end=find(ppm == sopt{grouped_copt_indeces{1,i}(j)}(end,1));
        ssel(i,ssel_beg:ssel_end)=NaN(1,length(ssel_beg:ssel_end));
    end
end

figure
plot(ppm,ini_est)
title('Initial estimates','fontweight','bold')
set(gca,'xdir','reverse')

end