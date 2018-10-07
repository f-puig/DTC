function [grouped_copt_indeces]=selNMRfeat(nmrspec, ref, ppm, copt, sopt, grouped_copt_indeces, new)
% If a new grouping is started, write:
% new = 1;
% If grouping is continued from the previous one, write:
% new = 0,

% If it is the first time running the program:
% grouped_copt_indeces={};

% nmrspec is an intensity matrix of N x M (samples x peak intensity)
% ref is the index of the reference nmr spectra to plot the SOPTs
% copt is the matrix of copt (C x N) (copt x samples)
% sopt is the cell containing sopt (1 x sopt),
% where each sopt is 2 x length(window)
% ppm is the vector of chemical shifts of nmrspec

% grouped_copt_indeces is a cell containing the indeces for each group
for i=1:length(sopt)
    sopt{1,i}=sopt{1,i}';
end

% 1) STACKED PLOT OF RAW SPECTRA
subplot(4,1,1)
plot(ppm, nmrspec(1,:),'Color',[0.72 0.82 0.96]);
title('^{1}H NMR dataset','fontweight','bold')
set(gca,'xdir','reverse')
hold on
color_sample = zeros(size(nmrspec,1),3);
color_sample(1,:)=[0.72 0.82 0.96];
for i=2:size(nmrspec,1)
    value=mod(i,50);
    tenth=mod(i,10);
    if (value<=10)
        color_sample(i,:)=[0 tenth/10 1-tenth/10];
    elseif (value>10 && value<=20)
        color_sample(i,:)=[tenth/10 1 0];
    elseif (value>20 && value<=30)
        color_sample(i,:)=[1 1-tenth/10 0];
    elseif (value>30 && value<=40)
        color_sample(i,:)=[1 0 tenth/10];
    elseif (value>40)
        color_sample(i,:)=[1-tenth/10 0 1];
    end
    plot(ppm, nmrspec(i,:),'Color',color_sample(i,:));
end
hold off

%CALCULATE GROUPINGS
assigned = cell2mat(grouped_copt_indeces);
cor_copt = corr(copt);
not_assigned = setdiff(1:1:size(cor_copt,1),assigned);
if length(not_assigned)> 1
    for i=1:size(cor_copt,1)
        cor_copt(i,i)=-1; %To avoid finding matches with themselves.
    end

    cor_copt2=cor_copt(not_assigned,not_assigned); % I create a copy of cor_copt in order to overwrite it.
    %cor_copt2 will reduce its dimensions at every iteration.

    if new == 1
        previous = [];
        match=[0 0]; % Vector that will keep the index of the two related features.
        [match(1),match(2)]=find(cor_copt==max(max(cor_copt2)),1,'first');
        not_assigned=not_assigned(not_assigned~=match(1));
        not_assigned=not_assigned(not_assigned~=match(2));
        %! I should have another way to recover the initial indeces
        fprintf('Higher correlation is found between %d and %d \n', match(1), match(2))
        fprintf('Correlation value is %f \n', max(max(cor_copt2)))
    else
        previous = grouped_copt_indeces{size(grouped_copt_indeces,2)};
        not_grouped = sort([not_assigned,previous]);
        list_correlations=[0 0 0];
        for i = 1:length(previous)
            cors1=[1-cor_copt(not_grouped,previous(i))';ones(1,length(not_grouped))*previous(i);not_grouped]'; % 1 - cor: sortrows sort in ascending order
            list_correlations=sortrows([list_correlations;cors1]);
        end
        % remove those pairs already picked-up.
        combinations = combnk(previous,2);
        for i = 1:size(combinations,1)
            remov1 = find(list_correlations(:,2)==combinations(i,1));
            if ~isempty(remov1)
                remov2 = find(list_correlations(remov1,3)==combinations(i,2));
                if ~isempty(remov2)
                    list_correlations = list_correlations(setdiff(1:size(list_correlations,1),remov1(remov2)),:);
                end
            end
        end

        for i = 1:size(combinations,1)
            remov1 = find(list_correlations(:,2)==combinations(i,2));
            if ~isempty(remov1)
                remov2 = find(list_correlations(remov1,3)==combinations(i,1));
                if ~isempty(remov2)
                    list_correlations = list_correlations(setdiff(1:size(list_correlations,1),remov1(remov2)),:);
                end
            end
        end

        match=list_correlations(2,3); % The first 'not_assigned' entity.
        fprintf('Higher correlation is found between %d and %d \n', list_correlations(2,2), list_correlations(2,3))
        fprintf('Correlation value is %f \n', 1-list_correlations(2,1))
        not_assigned=not_assigned(not_assigned~=match);
    end
    
    last = 0;
else
  fprintf('The only remaining feature is no. %d', not_assigned)
  match = not_assigned;
  if new == 1
      last = 1; 
      previous = not_assigned;
  end
end
t = size(grouped_copt_indeces,2)+1;

% **** PLOT ****
if last == 0
    features = [previous, match];
else
    features = [match];
end

subplot(4,1,2)
met_profile=zeros(1,length(ppm));
beg_met=zeros(1,length(features));
end_met=zeros(1,length(features));
for i=1:length(features)
   beg=find(ppm == sopt{features(i)}(1,1));
   ends=find(ppm == sopt{features(i)}(end,1));
   beg_met(i)=beg;
   end_met(i)=ends;
   labile=zeros(1,length(ppm));
   labile(beg:ends)= sopt{features(i)}(:,2)*copt(ref,features(i));
   met_profile=met_profile+labile;
end  
plot(ppm,met_profile)
%title(['\fontsize{16}black {\color{magenta}magenta ','\color[rgb]{0 .5 .5}teal \color{red}red} black again'],'interpreter','tex')
title({['Selected s_{k} spectral features for reference spectrum: no. ', ['\color[rgb]{0 1 0}',int2str(ref)]];['\fontsize{10}\color{blue}set of previously selected s_k      ','\fontsize{10}\color{red}new s_k']},'fontweight','bold')
set(gca,'xdir','reverse')
hold on
labile = sopt{features(end)}(:,2)*copt(ref,features(end));
plot(sopt{features(end)}(:,1),labile,'r');
hold off
        
subplot(4,1,3)
raw_sct=zeros(size(nmrspec,1),length(ppm));
 for j=1:size(nmrspec,1)
     for i=1:length(features)
        labile=zeros(1,length(ppm));
        labile(beg_met(i):end_met(i))= sopt{features(i)}(:,2)*copt(j,features(i));
        raw_sct(j,:)=raw_sct(j,:)+labile;
    end
end

plot(ppm, raw_sct(1,:),'Color',[0.72 0.82 0.96]);
title('Selected s_{k} spectral features for all samples','fontweight','bold')
set(gca,'xdir','reverse')
hold on
for i=2:size(nmrspec,1)
        plot(ppm, raw_sct(i,:),'Color',color_sample(i,:));
end
hold off

subplot(4,1,4)
if new == 1 && last == 0
    a=copt(:,match(1))./copt(:,match(2));
    scatter(1:size(nmrspec,1),copt(:,match(1))./copt(:,match(2)),10,color_sample,'filled');
    title([' Quotient between the highest correlated c_{k}. Value = ', num2str(median(a)), ' +/- ', num2str(std(a(~isinf(a))))],'fontweight','bold')
else
    a=copt(:,previous(end))./copt(:,match(1));
    scatter(1:size(nmrspec,1),copt(:,previous(end))./copt(:,match(1)),10,color_sample,'filled');
    title([' Quotient between the highest correlated c_{k}. Value = ', num2str(median(a)), ' +/- ', num2str(std(a(~isinf(a))))],'fontweight','bold')
end
xlabel('samples')
hold off
                      
% **** DECISION MAKING ****
if last == 1
    grouped_copt_indeces{t} = match(1);
elseif new == 1
    ans1 = input('Group NMR feature: 1-YES/0-NO');
        if ans1==1 % They belong to the same metabolite feature
          grouped_copt_indeces{t}=[match(1) match(2)];
        else % They belong to two different metabolite features
          grouped_copt_indeces{t}=match(1);
          grouped_copt_indeces{t+1}=match(2);  
        end
else %new == 0, there is only one index in match (the other is the already grouped features).
    ans1 = input('Group NMR feature: 1-YES/0-NO');
        if ans1==1 % They belong to the same metabolite feature
          grouped_copt_indeces{t-1} = [grouped_copt_indeces{t-1} match(1)];
        else % They belong to two different metabolite features
          grouped_copt_indeces{t} = match(1);
        end
end

not_assigned =setdiff(1:1:size(cor_copt,1),cell2mat(grouped_copt_indeces));

fprintf('Number of assigned features: %d', size(cor_copt,1)-length(not_assigned))
fprintf('Number of not-assigned features: %d', length(not_assigned))

if isempty(not_assigned)
   for i = 1:size(grouped_copt_indeces,2)
       grouped_copt_indeces{i}=sort(grouped_copt_indeces{i});
   end
disp('All features were grouped')
end    
% **** END ****
