function[mean_duration,std_duration] = validate_measures_manual_annotation_cl(action_select)
%Match manual annotation with automatic measures of posture and movements
%INPUT: type 'W' for walking, 'T' for turning, 'F' for freeze, 'R' for
%rear;

if strcmp(action_select, 'W')
    action_name{1} = 'Walking';
    pool_str{1}= action_name{1};
elseif strcmp(action_select, 'T')
    action_name = {'Turns Left','Turns Right','Bends left','Bends right'};
    pool_str{1}= 'Turns';
elseif strcmp(action_select, 'F') 
    action_name = {'Sniff + freeze','Whisk + freeze'}; 
    pool_str{1}= 'Freeze';
elseif strcmp(action_select, 'R') 
    action_name = {'Standing on hind legs','Climbing up walls'}; 
    pool_str{1}= 'Rearing';
else
    disp('select valid action');
    return;
end
 

%Note: frame numbers in excel files are from 1:76; in the original video
%they are: [120 195]. Thus I add 119 to action_frame (frame0=119)
which_action = [1:29];
Naction = numel(which_action); 
%load manual annotations
which_mouse = [1:6 8:25]; %one sheet per mouse in the xlsx file
Nmouse = numel(which_mouse);
action_frame_all = [];
frame_diff_all = [];
Amouse = []; Atrial = [];
frame0 = 119;
Naction = numel(action_name);
for n = 1:Nmouse
    [num,txt,raw] = xlsread('manual_annotations.xlsx',n);
    txt = txt(1,2:end); 
    for m = 1:Naction
        ind_action = [];
        for p = 1:size(txt,2)
            if numel(txt{p})>=numel(action_name{m})
                if strcmp(txt{p}(1:numel(action_name{m})),action_name{m})
                   ind_action = [ind_action p];
                end
            end
        end
        for p = 1:numel(ind_action)
            action_frame = num(:,ind_action(p));
            action_frame_start = num(1:3:end,ind_action(p));
            action_frame_end = num(2:3:end,ind_action(p));
            frame_diff = action_frame_end-action_frame_start;   
            frame_diff(frame_diff==0) = NaN;
            action_frame = action_frame(1:3:end,:)+frame0;    
            Ntrial = size(action_frame,1);
            action_frame_all = [action_frame_all; action_frame];
            frame_diff_all = [frame_diff_all; frame_diff];
            Amouse = [Amouse; which_mouse(n)*ones(Ntrial, 1)];
            Atrial = [Atrial; [1:Ntrial]'];
        end
    end
end
N = size(action_frame_all,1);
%%%%%%%%%%%%load images from videos%%%%%%%%%%%%%%%%%%
figure;
if strcmp(action_name{1},'Walking')
    imshow(imread('walk.tif'));
elseif strcmp(action_name{1},'Turns Left')
    imshow(imread('turn.tif'));
elseif strcmp(action_name{1},'Sniff + freeze')
    imshow(imread('freeze.tif'));
elseif strcmp(action_name{1},'Standing on hind legs')
    imshow(imread('rear.tif'));
end

%%%%%%%%%%%%load measures of postures and movements%%%%%%%%%%%%%%%%%%
%X is a 3dim matrix of Ntrial*Nframe*Ndim, where each dim is a posture or
%movement
load('measures_postures_movements.mat');
load('data_3D_all.mat','flag_val');
ind_ok = find(flag_val(:,1)==1);
X = X(ind_ok,:,:);
flag_val = flag_val(ind_ok,:);
Xmouse = flag_val(:,3);
Xtrial = flag_val(:,4);
Xstim = flag_val(:,5);
[Ntrial_all, Nframe, Ndim] = size(X);
%%%%%%%%%%%%combine action_frame with movements_quantified%%%%%%%%%%%%%%%%%%
action_val = cell(1);
L = 20; dL = -5;
frame_rate = 15;
tval = (dL:L+dL-1)/15;
mean_duration = nanmean(frame_diff_all)/frame_rate;
std_duration = nanstd(frame_diff_all)/frame_rate;
for n = 1:N
    comp = (Xtrial==Atrial(n)) + (Xmouse==Amouse(n));
    ind = find(comp == 2); 
    if numel(ind)
        if action_frame_all(n)>frame0 & frame_diff_all(n)>1 %& frame_diff_all(n,m)<11
            action_frame_all(n) = round(action_frame_all(n));
            ind_frame = action_frame_all(n)+dL:action_frame_all(n)+L-1+dL;
            Xtemp = squeeze(X(ind, ind_frame, :));
            action_val{1} = cat(3, action_val{1}, Xtemp);
            %action_val{m} = [action_val{m}; diff(Xtemp,1,1)];
        end
    end
end

%figure
dim_name = {'Re','Be','Bb','Lc','Fr','dRe','Br','dBe','dBb'};
cval = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 1 0.5 0; 0 1 0.5; 1 0 0.5];  
%
temp = cell(1,Ndim); 
for m = 1:Ndim
    temp{m} = squeeze(action_val{1}(:,m,:));
end
%
fig = figure; 
set(fig,'Position',[50 50 600 300]);
h(1) = subplot(1,2,1); hold on;
h(2) = subplot(1,2,2); hold on;
%h(3) = subplot(1,3,3); hold on;
%sub1
for m = 1:Ndim
    %if sum(m==rankval(1:2)) 
    if m < 4
        mtemp = mean(temp{m},2);        
        N = size(temp{m},2);
        subplot(h(1));
        stemp = std(temp{m},[],2)/sqrt(N); 
        errorbar(tval,mtemp,stemp,'.-','MarkerSize',18,'LineWidth',2,'Color',cval(m,:));
    end
    if ((m > 3) & (m < 6))|(m == 7)   
        mtemp = mean(temp{m},2);        
        N = size(temp{m},2);
        stemp = std(temp{m},[],2)/sqrt(N); 
        subplot(h(2));
        errorbar(tval,mtemp,stemp,'.-','MarkerSize',18,'LineWidth',2,'Color',cval(m,:));
    end
end
subplot(h(1)); set(h,'FontSize',12); 
title(pool_str{1},'FontSize',12);
%legend(dim_name([1:3]));
subplot(h(2)); set(h,'FontSize',12); 
title(['n=' num2str(N)],'FontSize',12);
%legend(dim_name([4 5 7]));   
for n = 1:2
    subplot(h(n));
    line([0 0],[0 1],'LineWidth',2,'Color',0.7*ones(1,3));
    line(mean_duration*[1 1],[0 1],'LineWidth',2,'Color',0.7*ones(1,3));
    xlim([tval(1) tval(end)]); ylim([0 1]);
end
 


