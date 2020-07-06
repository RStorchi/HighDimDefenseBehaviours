function[pval_latency,pval_rate] = estimate_rate_latency_primitives_cl(which_stim)
%Shows all results obtained with primitive analyses. 
%INPUT: 
%which_stim: 0, 1 or 2 (flash, loom or sound
%OUTPUT:
%Retun pvalues, calculated with kruskalwallis for changes in latency and 
%rate across the 7 clusters shown in Figure.5d 

%load the final clustering
load('clusters_VMM.mat','groups');
Npre = [5 4 2 2 1 1];
L = [2 3 5 6 10 15];
frame_rate = 15;
T = 2; 
%load the decoding results
[num,txt] = xlsread('export_vomm_results.xlsx');
time_steps = num(:,2);
correct = num(:,4); 
vomm_mem = num(:,5);
num_clus = num(:,6);

%show the performances
correct2 = reshape(correct,[10 6]);
fig = figure;
set(fig,'Position',[100 100 300 250]);
errorbar(L/frame_rate, mean(correct2),std(correct2),'LineWidth',2);
xlabel('#WindowDuration');
xlim([0 1.2]);

%find the best models
figure;
correct_sort = sort(unique(correct),'descend');
ind_best = [];
for n = 1:2
    temp = find(correct==correct_sort(n));
    ind_best = [ind_best; temp];
end
Nbest = numel(ind_best);
best_correct = 100*correct(ind_best);
best_num_clus = num_clus(ind_best);
best_time_steps = 2./time_steps(ind_best);
best_vomm_mem = vomm_mem(ind_best);
T = table(best_correct,best_time_steps,best_num_clus,best_vomm_mem);
uitable('Data',T{:,:},'ColumnName',{' %Accuracy ',' Duration ',' NumClus ',' VMM Order '})

%load measures of posture and movements & original groups & stimuli list
load('clustered_response','group');
load('../3Dreconstruction/measures_postures_movements.mat','X');
load('../3Dreconstruction/data_3D_all','Yrefined','flag_val');
load('../3Dreconstruction/SSM_trained','mean_pose');
ind_flag = find(flag_val(:,1) == 1);
X = X(ind_flag,:,:); Ndim = size(X,3);
flag_val = flag_val(ind_flag,:);
Yrefined = Yrefined(ind_flag);
ind_stim_vec0 = flag_val(:,5);
ind_stim_vec = ind_stim_vec0;
ind_stim_vec((ind_stim_vec0==1)|(ind_stim_vec0==2)) = 1;
ind_stim_vec((ind_stim_vec0==3)|(ind_stim_vec0==4)) = 2;


%select the best clustering results according to decoding
group_ts = groups{1,7}; Lseq = 15; Dframe = 2; Nts_pre = 5;

%divide primitives into pre-stimulus and response 
gval_ts = unique(group_ts(:));
Ng = numel(gval_ts);
group_ts_pre = group_ts(:,1:Nts_pre);
group_ts = group_ts(:,Nts_pre+1:end);
N = numel(group);

%find the most common sequence for each stimulus
temp = group_ts(ind_stim_vec==which_stim,:);
primitive_rate = hist(temp(:),gval_ts);
which_seq = gval_ts(find(primitive_rate == max(primitive_rate)));
primitive_rate_all = hist(group_ts(:),gval_ts);

%which_seq = 5;

%map clusters to 3D movements
Nframe = 30; frame_num0 = 150; 
trial_ind = []; frame_ind = []; count = 0; 
Xm = cell(1,Ndim);
for n = 1:N
    for m = 1:Lseq
        if group_ts(n,m) == which_seq
            count = count+1;
            frame_ind{count} = frame_num0 + (m-1)*Dframe + [1:Dframe];
            trial_ind = [trial_ind n];
            for p = 1:Ndim
                Xm{p} = [Xm{p}; mean(squeeze(X(n,frame_ind{count},p)))];
            end
        end
    end
end


%show average activity
Xm = horzcat(Xm{:});
fig = figure;
set(fig,'Position',[100 100 600 300]);
h = subplot(1,1,1); hold on; 
bar(mean(Xm),'FaceColor',0.66*ones(1,3),'BarWidth',0.5);
errorbar(mean(Xm),std(Xm),'.k','LineWidth',2); %/size(Xm,1)
line([0 Ndim+1],[0.5 0.5],'LineStyle','--','Color','k','LineWidth',2);
set(h,'XTick',1:Ndim);
set(h,'XTickLabel',{'Re','Be','Bb','Lc','Fr','dRe','Rt','dBe','dBb'});
ylabel('Normalized','FontSize',12);
set(h,'FontSize',12);
title(num2str(which_seq));
ylim([0 1.2]);

%show the poses
for n = 1:2
    Y = Yrefined{trial_ind(n)}(:,:,frame_ind{n});
    visualize_trial_3D_cl(Y,90); 
end

%reorder group, group_ts_pre,group_ts,ind_stim_vec
gval0 = unique(group); 
Ng0 = numel(gval0);
Ng0_count = hist(group,gval0);
[~,ind_sort] = sort(Ng0_count,'descend');
group_new = []; 
group_ts_new = [];
group_ts_pre_new = [];
ind_stim_vec_new = [];
for n = 1:Ng0
    ind = find(group==gval0(ind_sort(n)));
    group_new = [group_new; n*ones(numel(ind),1)]; 
    group_ts_new = [group_ts_new; group_ts(ind,:)];
    group_ts_pre_new = [group_ts_pre_new; group_ts_pre(ind,:)];
    ind_stim_vec_new = [ind_stim_vec_new; ind_stim_vec(ind)];
end
group = group_new;
group_ts = group_ts_new;
group_ts_pre = group_ts_pre_new;
ind_stim_vec = ind_stim_vec_new;

cval = generate_equally_saturated_colours(1,1,Ng,false);
%show classes 
estimate_rate_latency_primitives_cl_sub0(group,group_ts,cval);
%

%calculate latency
latency = estimate_rate_latency_primitives_cl_sub1(group_ts,which_seq);
%calculate rate
rate = estimate_rate_latency_primitives_cl_sub3(group_ts,which_seq);
%original dataset sum squares for latencies
[sum_sq_lat,mean_lat,std_lat] = estimate_rate_latency_primitives_cl_sub2(group,latency);
%original dataset sum squares for rates
[sum_sq_rate,mean_rate,std_rate] = estimate_rate_latency_primitives_cl_sub2(group,rate);

%stats for rate and latency
pval_latency = kruskalwallis(latency,group,'off');
pval_rate = kruskalwallis(rate,group,'off');

%figure latency distribution
estimate_rate_latency_primitives_cl_sub5(mean_lat,std_lat,['latency p=' num2str(pval_latency)],cval(which_seq,:))
%figure rate distribution
estimate_rate_latency_primitives_cl_sub5(mean_rate,std_rate,['rate p=' num2str(pval_rate)],cval(which_seq,:))

%figure compare flash with ongoing
estimate_rate_latency_primitives_cl_sub7(group_ts_pre,group_ts,ind_stim_vec,gval_ts,'evoked vs spontaneous');

%%%%%%%%%%%%%%%%%%%%%%%%%%%SUB_FUNCTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure groups colormap
function[] = estimate_rate_latency_primitives_cl_sub0(group,group_ts,cval)
gval = unique(group); Ng = numel(gval);
gval_ts = unique(group_ts); Ng_ts = numel(gval_ts);
[Ntrial,Nt] = size(group_ts);
group_ts_ordered = []; 
count = zeros(1,Ng);
for n = 1:Ng
    ind = find(group==gval(n));
    group_ts_ordered = [group_ts_ordered; group_ts(ind,:)];  
    if n == 1
        count(n) = numel(ind);
    else
        count(n) = count(n-1)+numel(ind);
    end
end
fig = figure; 
set(fig,'Position',[50 50 250 600]);
h = subplot(1,1,1); hold on; 
p = pcolor(group_ts_ordered); set(p,'LineStyle','none');
caxis([min(gval_ts) max(gval_ts)]);
%colormap(0.5*[ 1 1 0.65; 0 0 0; 1 0.65 1; 0.65 1 1]);
colormap(cval);
for n = 1:Ng
    line([0 Nt],[count(n) count(n)]+1,'LineWidth',2,'Color',[0 0 0]);
end
set(h,'XTick',[],'YTick',[]);
xlim([1 Nt]); ylim([1 Ntrial]);
%xlabel('TimeStep');  ylabel('#Trial'); 

%find the latency
function[latency] = estimate_rate_latency_primitives_cl_sub1(group_ts,which_seq)
[Ntrial,Nstep] = size(group_ts);
Lseq = length(which_seq);
latency = NaN*ones(1,Ntrial);
for n = 1:Ntrial
    latency_trial = false(1,Nstep-Lseq+1);
    for m = 1:Nstep-Lseq+1
        latency_trial(m) = isequal(group_ts(n,m:m+Lseq-1),which_seq);
    end
    if sum(latency_trial)
       %latency(n) = find(latency_trial,1); %first element latency
       latency(n) = mean(find(latency_trial)); %mean latency
    end
end
    
%calculate the variance explained 
function[sum_sq,mean_latency,std_latency] = estimate_rate_latency_primitives_cl_sub2(group,latency)
ind = find(~isnan(latency));
latency = latency(ind); 
group = group(ind); 
N = numel(ind);
gval = unique(group);
Ng = numel(gval);
for n = 1:Ng
    mean_latency(n) = nanmean(latency(group==gval(n)));
    std_latency(n) = nanstd(latency(group==gval(n)));
end
sum_sq = 0;
for n = 1:N
    sum_sq = sum_sq+(latency(n)-mean_latency(find(group(n)==gval))).^2;
end

%find the rate
function[rate] = estimate_rate_latency_primitives_cl_sub3(group_ts,which_seq)
[Ntrial,Nstep] = size(group_ts);
Lseq = length(which_seq);
rate = zeros(1,Ntrial);
for n = 1:Ntrial
    rate_trial = false(1,Nstep-Lseq+1);
    for m = 1:Nstep-Lseq+1
        rate_trial(m) = isequal(group_ts(n,m:m+Lseq-1),which_seq);
    end
    if sum(rate_trial)
       rate(n) = sum(rate_trial);
    end
end

%figure rate and latency histograms
function[] = estimate_rate_latency_primitives_cl_sub5(x,y,tt,cval)
fig = figure; 
set(fig,'Position',[50 50 300 300]);
hold on; 
bar(x,'BarWidth',0.5,'FaceColor',cval);
errorbar(x, y,'.','LineWidth',2,'Color','k')
xlim([0 numel(x)+1]); ylim([min(x-2*y) max(x+2*y)]);
xlabel('#group');  title(tt); 

function[] = estimate_rate_latency_primitives_cl_sub7(group_ts_pre,group_ts,ind_stim_vec,gval_ts,tt)
fig = figure; hold on;
set(fig,'Position',[100 100 300 250]);
Ng = numel(gval_ts);
sval = unique(ind_stim_vec);
Ns = numel(sval);
for n = 1:Ns
    temp = group_ts(ind_stim_vec==sval(n),:);
    count(n,:) = hist(temp(:),gval_ts);
end
count_pre = hist(group_ts_pre(:),gval_ts);
rng(1);
cval = [0 0 0; 0 0 1; 1 0 0];
for n = 1:Ng
    for m = 1:Ns
        bar(2*n-0.25*m,count(m,n),'BarWidth',0.25,'FaceColor',cval(m,:));
    end
    bar(2*n,count_pre(n),'BarWidth',0.25,'FaceColor',0.66*ones(1,3));
end
xlabel('#primitive'); ylabel('count');
title(tt);
%to test with pearson's chi2
%[h,p] = chi2gof(count(n,:),'Expected',count_pre)
