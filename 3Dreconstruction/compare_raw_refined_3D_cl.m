function[] = compare_raw_refined_3D_cl()


%threshold to detect outliers 
THoutliers = 5; 

%load raw and refined 3D data
load('data_3D_all.mat');
N = numel(Yraw);

%keep frames between -1s to +3s cut others
%stim onset at frame 150
for n = 1:N
    Yraw{n} = Yraw{n}(:,:,135:196);
    if numel(Yrefined{n})>0
        Yrefined{n} = Yrefined{n}(:,:,135:196);
    end
end
[Np,~,Nframe] = size(Yraw{1});

%load SSM
load('SSM_trained','mean_pose','P');
%reshape mean pose
mean_pose_3d = reshape(mean_pose,Np,3);

%all poses to the mean pose
Yraw_aligned = cell(1,N);
Yrefined_aligned = cell(1,N);
count = 0;
%shape parameters
proj_raw = zeros(3*Np,N*Nframe);
proj_refined = zeros(3*Np,N*Nframe);
%Reviewer 4 measures
nose_tail = zeros(1,N*Nframe);
neck_tail = zeros(1,N*Nframe);
head_angle = zeros(1,N*Nframe);
head_angle_yz = zeros(1,N*Nframe);
%distances from mean pose
dist_raw = zeros(1,N*Nframe);
dist_refined = zeros(1,N*Nframe);
%others
outliers_raw = false(1,N*Nframe);
outliers_refined = false(1,N*Nframe);
which_trial = zeros(1,N*Nframe);
which_frame = zeros(1,N*Nframe);
for n = 1:N
    Yraw_aligned{n} = NaN*zeros(Np,3,Nframe);
    Yrefined_aligned{n} = NaN*zeros(Np,3,Nframe);
    for m = 1:Nframe
        %%%%%update
        count = count+1;
        which_trial(count) = n;
        which_frame(count) = m;
        %%%%%%%raw%%%%%%%
        if size(Yraw{n})
            X = Yraw{n}(:,:,m);
            %align
            [~, ~, t] = procrustes(mean_pose_3d, X, 'Scaling', false,'Reflection',false);
            Yraw_aligned{n}(:,:,m) = t.b*X * t.T + t.c;
            X_aligned_1d = reshape(Yraw_aligned{n}(:,:,m),3*Np,1);            
            %calculate projection on eigenposes
            proj_raw(:,count) = P'*(X_aligned_1d - mean_pose);
            dist_raw(count) = norm(X_aligned_1d - mean_pose);
            %check for outliers
            if dist_raw(count) > THoutliers
                outliers_raw(count) = true;
            end
        else
           outliers_raw(count) = true;
           dist_raw(count) = Inf;
           
        end
        %%%%%%%refined%%%%%%% %   if prod(size(X))==3*Np
        if size(Yrefined{n})
            X = Yrefined{n}(:,:,m);
            %align
            [~, ~, t] = procrustes(mean_pose_3d, X, 'Scaling', false,'Reflection',false);
            Yrefined_aligned{n}(:,:,m) = t.b*X * t.T + t.c;
            X_aligned_1d = reshape(Yrefined_aligned{n}(:,:,m),3*Np,1); 
            %calculate Reviewer 4 measures
            temp_xyz = Yrefined_aligned{n}(:,:,m);
            %head angle xy
            temp_xy = temp_xyz(:,1:2); 
            nose_vec = temp_xy(1,:)-temp_xy(4,:); 
            tail_vec = temp_xy(5,:)-temp_xy(4,:); 
            dot_val = nose_vec*tail_vec';
            det_val = det([nose_vec; tail_vec]);
            head_angle(count) = atan2(det_val, dot_val);
            %head angle yz
            temp_yz = temp_xyz(:,2:3); 
            nose_vec = temp_yz(1,:)-temp_yz(4,:); 
            tail_vec = temp_yz(5,:)-temp_yz(4,:); 
            dot_val = nose_vec*tail_vec';
            det_val = det([nose_vec; tail_vec]);
            head_angle_yz(count) = atan2(det_val, dot_val);
            %nose tail
            nose_tail(count) = norm(temp_xyz(1,:)-temp_xyz(5,:));
            %neck tail
            neck_tail(count) = norm(temp_xyz(4,:)-temp_xyz(5,:));
            %calculate projection on eigenposes
            proj_refined(:,count) = P'*(X_aligned_1d - mean_pose);
            dist_refined(count) = norm(X_aligned_1d - mean_pose);
            %check for outliers
            if dist_refined(count) > THoutliers
                outliers_refined(count) = true;
            end
        else
           outliers_refined(count) = true;
           dist_refined(count) = Inf;
           
        end
    end
end
%concatenate 3D poses
Yraw_aligned = cat(3,Yraw_aligned{:});
Yrefined_aligned = cat(3,Yrefined_aligned{:});

 
%count N and outliers
N = size(proj_raw,2);  
Noutliers_raw = sum(outliers_raw);
Noutliers_refined = sum(outliers_refined);
Noutliers_trials_raw = numel(unique(which_trial(outliers_raw)));
Noutliers_trials_refined = numel(unique(which_trial(outliers_refined)));

%show figure outliers
%raw
dist_raw(dist_raw>100) = 100;
fig = figure; 
set(fig,'Position',[200 200 350 250]);
h = subplot(1,1,1); hold on; 
hist(log10(dist_raw),[-1.5:0.05:2]);
line(log10([THoutliers THoutliers]),[0 4000],'Color','k','LineWidth',2);
xlim([-1.5 2.25]); ylim([-0.1 4500]); 
title(['n = ' num2str(N) '; out = ' num2str(Noutliers_raw) '; Ntrial = ' num2str(Noutliers_trials_raw)]);
%refined
dist_refined(dist_refined>100) = 100;
fig = figure; 
set(fig,'Position',[200 200 350 250]);
h = subplot(1,1,1); hold on; 
hist(log10(dist_refined),[-1.5:0.05:2]);
line(log10([THoutliers THoutliers]),[0 4000],'Color','k','LineWidth',2);
xlim([-1.5 2.25]); ylim([-0.1 4500]); 
title(['n = ' num2str(N) '; out = ' num2str(Noutliers_refined) '; Ntrial = ' num2str(Noutliers_trials_refined)]);

%look at projections
proj_raw = proj_raw';
proj_refined = proj_refined';
proj_raw = proj_raw(~outliers_raw,:);
proj_refined = proj_refined(~outliers_refined,:);
dist_raw = dist_raw(~outliers_raw);
dist_refined = dist_refined(~outliers_refined);
%take principal components 
[coeff_raw, score_raw, vars_raw] = pca(proj_raw);
[coeff_refined, score_refined, vars_refined] = pca(proj_refined);
%figure distances
fig = figure; 
h = subplot(1,1,1); hold on;
set(fig,'Position',[200 200 300 250]);
bar([mean(dist_raw) mean(dist_refined)],'BarWidth',0.5,'FaceColor',0.66*ones(1,3));
errorbar([mean(dist_raw) mean(dist_refined)],[std(dist_raw) std(dist_refined)],'.k','LineWidth',2);
ylabel('Distance(cm)','FontSize',14);
set(h,'XTick',[1, 2]); set(h,'XTickLabel',{'Raw' 'Ref'},'FontSize',14);
title(['pval = ' num2str(ranksum(dist_raw,dist_refined))]);
set(h,'FontSize',12);
%figure components
fig = figure; 
h = subplot(1,1,1); hold on;
set(fig,'Position',[200 200 350 300]);
plot(cumsum(vars_raw)/sum(vars_raw),'k.-','LineWidth',2,'MarkerSize',15); 
plot(cumsum(vars_refined)/sum(vars_refined),'b.-','LineWidth',2,'MarkerSize',15); 
xlim([0 10]); ylim([0.3 1.05]); 
xlabel('PCs','FontSize',14);
ylabel('Variance Explained','FontSize',14);
set(h,'FontSize',12);

%%%%%%%%plot all 3d aligned points%%%%%%%%
Yraw_aligned = Yraw_aligned(:,:,~outliers_raw);
Yrefined_aligned = Yrefined_aligned(:,:,~outliers_refined);
validate_manual_review_v1_fig1(Yraw_aligned);
validate_manual_review_v1_fig1(Yrefined_aligned);

%show correlations with Reviewer 4 measures and the first 3 shape
%parameters
nose_tail = nose_tail(~outliers_refined);
neck_tail = neck_tail(~outliers_refined);
head_angle = head_angle(~outliers_refined);
head_angle_yz = head_angle_yz(~outliers_refined);
head_angle = 180*wrapTo2Pi(head_angle)/pi;
head_angle_yz = 180*wrapTo2Pi(head_angle_yz)/pi;
%
fig = figure; 
set(fig,'Position',[50 50 700 400]);
N = numel(nose_tail);
Nred = min(N,1000);
indred = randperm(N); 
indred = indred(1:Nred);
%nose-tail to b1-b3
for n = 1:3
    h = subplot(2,3,n); hold on;
    plot(nose_tail(indred),proj_refined(indred,n),'.'); 
    xlim([0 20]);
    R = corr2(nose_tail',proj_refined(:,n));
    title(num2str(R));
    xlim([5 15]); ylim([-2.5 2.5]);
    xlabel('nose-tail(cm)','FontSize',12);
    ylabel(['b' num2str(n)],'FontSize',12);
end
%head-angle xy to b1-b3
for n = 1:3
    h = subplot(2,3,n+3); hold on;
    plot(head_angle(indred),proj_refined(indred,n),'.'); 
    R = corr2(head_angle',proj_refined(:,n));
    title(num2str(R));
    xlim([130 230]); ylim([-2.5 2.5]);
    xlabel('head angle XY(deg)','FontSize',12);
    ylabel(['b' num2str(n)],'FontSize',12);
end
%%%%%%%%%%%%%%%%%%%%
fig = figure; 
set(fig,'Position',[50 50 700 400]);
%neck-tail xy to b1-b3
for n = 1:3
    h = subplot(2,3,n); hold on;
    plot(neck_tail(indred),proj_refined(indred,n),'.'); 
    xlim([0 20]);
    R = corr2(neck_tail',proj_refined(:,n));
    title(num2str(R));
    xlim([2 10]); ylim([-2.5 2.5]);
    xlabel('neck-tail(cm)','FontSize',12);
    ylabel(['b' num2str(n)],'FontSize',12);
end
%head-angle yz to b1-b3
for n = 1:3
    h = subplot(2,3,n+3); hold on;
    plot(head_angle_yz(indred),proj_refined(indred,n),'.'); 
    R = corr2(head_angle_yz',proj_refined(:,n));
    title(num2str(R));
    xlim([130 300]); ylim([-2.5 2.5]);
    xlabel('head angle YZ(deg)','FontSize',12);
    ylabel(['b' num2str(n)],'FontSize',12);
end
'a'
%figures
function[] = validate_manual_review_v1_fig1(poses)
[Np,~,N] = size(poses);
ind = randperm(N); 
N = min(N,10000);
poses = poses(:,:,ind(1:N));
%
fig0 = figure; 
%set(fig0,'Position',[50 50 300 600]);
set(fig0,'Position',[50 50 600 300]);
h = subplot(1,1,1); hold on;
cval = [1 0 0; 0 0 1; 0 0 0.5; 0 1 0; 0.85 0.85 0];
for m = 1:Np
    %plot3(squeeze(poses(m,1,:)),squeeze(poses(m,2,:)),squeeze(poses(m,3,:)),'.','MarkerSize',10,'Color',cval(m,:));
    plot(squeeze(poses(m,2,:)),squeeze(poses(m,1,:)),'.','MarkerSize',10,'Color',cval(m,:));
end
%xlabel('X(cm)','FontSize',14);
%ylabel('Y(cm)','FontSize',14);
%zlabel('Z(cm)','FontSize',14);
ylabel('X(cm)','FontSize',14);
xlabel('Y(cm)','FontSize',14);
set(h,'FontSize',12);
%view([0 90]); 
%xlim([-4 4]); ylim([-8 8]);
ylim([-4 4]); xlim([-8 8]);
%
fig1 = figure; 
set(fig1,'Position',[50 50 600 300]);
h = subplot(1,1,1); hold on;
cval = [1 0 0; 0 0 1; 0 0 0.5; 0 1 0; 0.85 0.85 0];
for m = 1:Np
    %plot3(squeeze(poses(m,1,:)),squeeze(poses(m,2,:)),squeeze(poses(m,3,:)),'.','MarkerSize',10,'Color',cval(m,:));
    plot(squeeze(poses(m,2,:)),squeeze(poses(m,3,:)),'.','MarkerSize',10,'Color',cval(m,:));
end
%xlabel('X(cm)','FontSize',14);
%ylabel('Y(cm)','FontSize',14);
%zlabel('Z(cm)','FontSize',14);
xlabel('Y(cm)','FontSize',14);
ylabel('Z(cm)','FontSize',14);
set(h,'FontSize',12);
%view([90 0]); 
%zlim([-3 3]); ylim([-8 8]);
xlim([-8 8]); ylim([-3 3]); 
