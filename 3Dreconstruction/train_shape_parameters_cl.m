function[b,P,mean_pose,lambda,poses] = train_shape_parameters_cl(make_fig,make_movie)
%INPUT: 
%make_fig: if true two figures showing the 3D aligned training points are displayed 
%make_movie: if true the script plays a short movie showing changes in body shape along directions of the eigenposes (Supplementary Movie 2 in the paper) 
%OUTPUT:
%b: shape parameters
%P: eigenposes
%mean_pose: mean of the aligned poses
%lambda: eigenvalues associated with eigenposes
%poses: all aligned 3D poses

%load the template poses 
load('SSM_train_data.mat','good_pose','nose_3d','right_ear_3d','left_ear_3d','neck_3d','tail_3d');
ind = find(good_pose); %use only poses approved by human inspection for correctness of body-landmark detection 
N = length(ind); 
Np = 5;
nose_3d = nose_3d(ind); right_ear_3d = right_ear_3d(ind); left_ear_3d = left_ear_3d(ind); 
neck_3d = neck_3d(ind); tail_3d = tail_3d(ind);
poses = zeros(Np,3,N);
for n = 1:N
    poses(1,:,n) = nose_3d{n}(1:3)';
    poses(2,:,n) = left_ear_3d{n}(1:3)';
    poses(3,:,n) = right_ear_3d{n}(1:3)';
    poses(4,:,n) = neck_3d{n}(1:3)';
    poses(5,:,n) = tail_3d{n}(1:3)';
end

%align all poses to the first pose
X = [0 3.5 0; -1 2 0; 1 2 0; 0 1 0; 0 -3.5 0];
for n = 1:N
    Y = squeeze(poses(:,1:3,n)); 
    [~, ~, t] = procrustes(X, Y, 'Scaling', false,'Reflection',false);
    poses(:,1:3,n) = t.b*Y * t.T + t.c;
end


%plot all the poses in the dataset
if make_fig
    %
    fig0 = figure; 
    set(fig0,'Position',[200 200 300 600]);
    hold on;
    cval = [1 0 0; 0 0 1; 0 0 0.5; 0 1 0; 0.85 0.85 0];
    for n = 1:N
        for m = 1:Np
            plot3(poses(m,1,n),poses(m,2,n),poses(m,3,n),'.','MarkerSize',10,'Color',cval(m,:));
        end
    end
    view([0 90]); xlim([-4 4]); ylim([-8 8]);
    %
    fig1 = figure; 
    set(fig1,'Position',[200 200 600 300]);
    hold on;
    cval = [1 0 0; 0 0 1; 0 0 0.5; 0 1 0; 0.85 0.85 0];
    for n = 1:N
        for m = 1:Np
            plot3(poses(m,1,n),poses(m,2,n),poses(m,3,n),'.','MarkerSize',10,'Color',cval(m,:));
        end
    end
    view([90 0]); zlim([-3 3]); ylim([-8 8]);
end

%get the shape parameters
%reshape poses
poses = reshape(poses,Np*3,N);
%calculate covariance matrix
mean_pose = mean(poses')';
S = cov(poses'-repmat(mean_pose',N,1));
%get eigenvectors/values
[P,lambda] = eig(S,'vector');
%sort eigenvalues/vector, keep the relevant ones
[lambda,indsort] = sort(lambda,'descend');
P = P(:,indsort);
%save trained SSM
save('SSM_trained','mean_pose','P','lambda');
%reduce the number of eigenposes 
Nshape = min(find((cumsum(lambda)/sum(lambda))>0.9));
P = P(:,1:Nshape);
lambda = lambda(1:Nshape);
%calculate shape parameters
Dposes = poses - repmat(mean_pose,1,N);
b = P'*Dposes;

%make a movie for change in body shape along directions of the eigenposes
if make_movie
    Nmovie = 100;
    fig = figure; 
    set(fig,'Position',[200 200 800 400]);
    h(1) = subplot(1,2,1);
    h(2) = subplot(1,2,2);
    cval = [1 0 0; 0 0 1; 0 0 0.5; 0 1 0; 0.85 0.85 0];
    for nn = 1:Nshape %first 3 eigenposes
        ind0_image = (nn-1)*Nmovie;
        b_movie = zeros(Nshape,Nmovie);
        b_movie(nn,:) = sqrt(6)*lambda(nn)*sin(2*pi*[1:Nmovie]/(50));
        poses_movie = zeros(Np,3,Nmovie);
        for n = 1:Nmovie
            temp = mean_pose + P*b_movie(:,n);
            poses_movie(:,:,n) = reshape(temp,Np,3);
        end
        for n = 1:Nmovie
            %
            subplot(h(1)); hold on;
            title(['Eigenpose: ' num2str(nn)],'FontSize',16);
            for m = 1:Np
                plot3(poses_movie(m,1,n),poses_movie(m,2,n),poses_movie(m,3,n),'.','MarkerSize',18,'Color',cval(m,:));
            end
            xlim([min(min(poses_movie(:,1,:))) max(max(poses_movie(:,1,:)))]);
            ylim([min(min(poses_movie(:,2,:))) max(max(poses_movie(:,2,:)))]);
            zlim([min(min(poses_movie(:,3,:))) max(max(poses_movie(:,3,:)))]);
            xlabel('X','FontSize',16);ylabel('Y','FontSize',16);zlabel('Z','FontSize',16);
            view([0 90]); xlim([-8 8]); ylim([-8 8]);
            %
            subplot(h(2)); hold on;
            title(['Eigenpose: ' num2str(nn)],'FontSize',16);
            for m = 1:Np
                plot3(poses_movie(m,1,n),poses_movie(m,2,n),poses_movie(m,3,n),'.','MarkerSize',18,'Color',cval(m,:));
            end
            xlim([min(min(poses_movie(:,1,:))) max(max(poses_movie(:,1,:)))]);
            ylim([min(min(poses_movie(:,2,:))) max(max(poses_movie(:,2,:)))]);
            zlim([min(min(poses_movie(:,3,:))) max(max(poses_movie(:,3,:)))]);
            xlabel('X','FontSize',16);ylabel('Y','FontSize',16);zlabel('Z','FontSize',16);
            view([90 0]); zlim([-8 8]); ylim([-8 8]);
            %
            pause(0.025); 
            if n < Nmovie
                subplot(h(1)); cla;
                subplot(h(2)); cla;
            end
        end
    end
end


