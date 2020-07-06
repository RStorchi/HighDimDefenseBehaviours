function[Y,Y_final] = refine_3D_cl(which_mouse, which_trial, which_frames, showfigs)
%This function performs a refined 3D reconstruction based on SSM
%INPUT: 
%which_mouse: mouse number (1-30, except 7)
%which_trial: trial number (1-18)
%which_frames: to select frame interval (1 to 320)
%show_figs: if true makes a movie of raw and refined reconstruction
%OUTPUT:
%Y: raw 3D (dimensions: dim1 = body landmark(5); dim2 = X,Y,Z coordinates;
%dim3 = frame index)
%Y_final: refined 3D (dimensions: dim1 = body landmark(5); dim2 = X,Y,Z coordinates;
%dim3 = frame index)
%EXAMPLE
%to reconstruct frames 140:180 from mouse 1 at trial 3 use:
%refine_3D_cl(1, 3, 140:180, true);
%to reconstruct all frames in a trial use []:
%refine_3D_cl(1, 3, [], true);

%%%%%%%%%%%%%INIT PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%
%smoothing window
wind = [0.2 0.6 0.2];
%regularization parameter
alpha_reg = 0.001;

%load the 3d camera calibration
load('../3Dcalbration/Pcal.mat');
Ncam = numel(P);

%load the 3D data
load('data_3D_all.mat','Yraw','flag_val');
%extract the trial of interest
ind = find((flag_val(:,3)==which_mouse)&(flag_val(:,4)==which_trial));
if numel(which_frames)>0
    Yraw = Yraw{ind}(:,:,which_frames);
else
    Yraw = Yraw{ind};
end
Y = Yraw;
[Np,~,Nframe] = size(Y);

%load shape parameters
[b,V,mean_pose1,lambda] = train_shape_parameters_cl(0 ,0);
Nshape = length(lambda);
mean_pose = reshape(mean_pose1,Np,3);

%perform first robust reconstruction without shape constraints
alpha_reg_init = 0;
[b0, dist0, R0, t0] = reconstruct_procrustes_fminunc_cl(Y,mean_pose,V,lambda,alpha_reg_init);

%detect outliers 
THdist = 2;
THb = 6;
d_outlier_flag = dist0>THdist;
b_outlier_flag = false(1,Nframe);
for n = 1:Nshape
    b_outlier_flag = b_outlier_flag | ((b0(n,:).^2)/lambda(n))>THb;
end
outlier_flag = d_outlier_flag | b_outlier_flag;
Noutlier = sum(outlier_flag);
b1 = b0;
b1(:,outlier_flag) = NaN;
ind_outlier = find(outlier_flag);

%correct for outliers
for n = 1:Nshape
    b1(n,:) = fillmissing_cl(b1(n,:),5,'pchip',1);
end

%identify the outlier points by leaving 2 points out
d = zeros(Np,Noutlier);
ind_comb = combnk([1:Np],3)';
Ncomb = size(ind_comb,2);
dist_outlier = zeros(Ncomb,Noutlier);
for n = 1:Ncomb
    for m = 1:Noutlier
        [~, Z] = procrustes(mean_pose(ind_comb(:,n),:),Y(ind_comb(:,n),:,ind_outlier(m)),'Reflection',false, 'Scaling',false);
        dist_outlier(n,m) = norm(mean_pose(ind_comb(:,n),:)-Z);
    end
end

%substitute the outlier points
for n = 1:Noutlier
    best_comb = ind_comb(:,find(dist_outlier(:,n) == min(dist_outlier(:,n))));
    temp1 = mean_pose1 + V*b1(:,ind_outlier(n));
    temp = reshape(temp1,Np,3);
    %
    [~, ~, tr] = procrustes(Y(best_comb,:,ind_outlier(n)),temp(best_comb,:),'Reflection',false, 'Scaling',false);
    Y(:,:,ind_outlier(n)) = tr.b * temp * tr.T + repmat(tr.c(1,:),Np,1);
    %
    [~, Z, tr] = procrustes(mean_pose,Y(:,:,ind_outlier(n)),'Reflection',false, 'Scaling',false);
    %
    Z1 = reshape(Z,3*Np,1);
    b1(:,ind_outlier(n)) = V'*(Z1-mean_pose1);
end

%perform first robust reconstruction with shape constraints
[b2, dist2, R2, t2] = reconstruct_procrustes_fminunc_cl(Y,mean_pose,V,lambda,alpha_reg);

%smooth the translations
for n = 1:3
    t_final(n,:) = filtfilt(wind,1,t2(n,:));
end
%smooth the shape parameters
for n = 1:Nshape
    b_final(n,:) = filtfilt(wind,1,b2(n,:));
end
%smooth the rotation matrices
R_final = zeros(3,3,Nframe);
for n = 1:3
    for m = 1:3
        R_final(n,m,:) = filtfilt(wind,1,R2(n,m,:));
    end
end
for n = 1:Nframe
    [Ur,Sr,Vr] = svd(R_final(:,:,n));
    R_final(:,:,n) = Ur*Vr'; 
end

%generate the final 3D reconstruction
Y_final = zeros(Np,3,Nframe);
for n = 1:Nframe
    Y_final(:,:,n) = reshape(mean_pose1+V*b_final(:,n),Np,3)*R_final(:,:,n)+repmat(t_final(:,n)',Np,1);
end

%show initial and final reconstructions
if showfigs
   show_trial2(Yraw,Y_final);
end



%2 Plots
function[] = show_trial2(Y,Yfinal)
[Np,~,Nframe] = size(Y);
fig = figure; 
set(fig,'Position',[100 100 1200 600]);
h(1) = subplot(1,2,1);
h(2) = subplot(1,2,2);
cval = 'rbbgykc';
for n = 1:Nframe
    %
    subplot(h(1)); cla; hold on;
    for m = 1:Np
        plot3(Y(m,1,n),Y(m,2,n),Y(m,3,n),'.','MarkerSize',15,'Color',cval(m));
    end
    line([Y(1,1,n) Y(2,1,n)],[Y(1,2,n) Y(2,2,n)],[Y(1,3,n) Y(2,3,n)],'Color','k','LineWidth',2);
    line([Y(1,1,n) Y(3,1,n)],[Y(1,2,n) Y(3,2,n)],[Y(1,3,n) Y(3,3,n)],'Color','k','LineWidth',2);
    line([Y(2,1,n) Y(4,1,n)],[Y(2,2,n) Y(4,2,n)],[Y(2,3,n) Y(4,3,n)],'Color','k','LineWidth',2);
    line([Y(3,1,n) Y(4,1,n)],[Y(3,2,n) Y(4,2,n)],[Y(3,3,n) Y(4,3,n)],'Color','k','LineWidth',2);
    line([Y(4,1,n) Y(5,1,n)],[Y(4,2,n) Y(5,2,n)],[Y(4,3,n) Y(5,3,n)],'Color','k','LineWidth',2);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view([220 10]); 
    xlim([-30 30]); ylim([-30 30]); zlim([-15 15]); 
    line([-20 20],[20 20],[0 0]);
    line([-20 20],-[20 20],[0 0]);
    line(-[20 20],[-20 20],[0 0]);
    line([20 20],[-20 20],[0 0]);
    title('Raw 3D','FontSize',14);
    %
    subplot(h(2)); cla; hold on;
    for m = 1:Np
        plot3(Yfinal(m,1,n),Yfinal(m,2,n),Yfinal(m,3,n),'.','MarkerSize',15,'Color',cval(m));
    end
    line([Yfinal(1,1,n) Yfinal(2,1,n)],[Yfinal(1,2,n) Yfinal(2,2,n)],[Yfinal(1,3,n) Yfinal(2,3,n)],'Color','k','LineWidth',2);
    line([Yfinal(1,1,n) Yfinal(3,1,n)],[Yfinal(1,2,n) Yfinal(3,2,n)],[Yfinal(1,3,n) Yfinal(3,3,n)],'Color','k','LineWidth',2);
    line([Yfinal(2,1,n) Yfinal(4,1,n)],[Yfinal(2,2,n) Yfinal(4,2,n)],[Yfinal(2,3,n) Yfinal(4,3,n)],'Color','k','LineWidth',2);
    line([Yfinal(3,1,n) Yfinal(4,1,n)],[Yfinal(3,2,n) Yfinal(4,2,n)],[Yfinal(3,3,n) Yfinal(4,3,n)],'Color','k','LineWidth',2);
    line([Yfinal(4,1,n) Yfinal(5,1,n)],[Yfinal(4,2,n) Yfinal(5,2,n)],[Yfinal(4,3,n) Yfinal(5,3,n)],'Color','k','LineWidth',2);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view([220 10]); 
    xlim([-30 30]); ylim([-30 30]); zlim([-15 15]); 
    line([-20 20],[20 20],[0 0]);
    line([-20 20],-[20 20],[0 0]);
    line(-[20 20],[-20 20],[0 0]);
    line([20 20],[-20 20],[0 0]);
    title('Refined 3D','FontSize',14);
    %
    pause(1/15);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[y] = fillmissing_cl(y,L,which_method,Npoly)
%use this function to fill missing values in shape parameters
N = numel(y);
ind_miss = find(isnan(y),1);
ncycle = 0;
while length(ind_miss) > 0
    ncycle = ncycle+1;
    if ncycle > 100000
       break
    end
    %randomise the order of ind_miss
    ind_miss = ind_miss(randperm(numel(ind_miss)));
    %
    ind_num = find(~isnan(y));
    ind_pre = ind_num(max(find(ind_num<ind_miss)));
    ind_post = ind_num(min(find(ind_num>ind_miss)));
    %generate the buffer 
    ind_pre0 = max(ind_pre-L,1);
    ind_post0 = min(ind_post+L,N);
    if length(ind_pre0)==0
       ind_pre0 = 1;
    end
    if length(ind_post0)==0
       ind_post0 = N;
    end
    ind_vec = ind_pre0:ind_post0;
    ind_num_vec = ind_vec(~isnan(y(ind_vec)));
    ind_miss_vec = ind_vec(isnan(y(ind_vec)));
    %interpolate
    if numel(ind_num_vec) >= 2
        %spline
        if strcmp(which_method,'spline')
            y(ind_miss) = spline(ind_num_vec, y(ind_num_vec), ind_miss);
        %pchip
        elseif strcmp(which_method,'pchip')
            y(ind_miss) = pchip(ind_num_vec, y(ind_num_vec), ind_miss);
        %poly
        elseif strcmp(which_method,'poly')
            p = polyfit(ind_num_vec, y(ind_num_vec),Npoly);
            y(ind_miss) = polyval(p,ind_miss);
        else 
            disp('provided methods: spline, pchip, poly');
            return
        end
    end
    %updated ind_miss
    ind_miss = find(isnan(y),1);
end



