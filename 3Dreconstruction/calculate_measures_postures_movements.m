function[] = calculate_measures_postures_movements()
%This script use the refined 3D reconstruction to generate measures of
%postures and movements. Data are concatenated in a 3D array where dim1 =
%trials, dim2 = time points, dim3 = measures of postures/movements.
%Such measures are:
%POSTURES
%1) rear
%2) body elongation
%3) body bending
%MOVEMENTS:
%4) locomotion
%5) freeze
%6) Drear
%7) Drotation
%8) Dbody_elongation
%9) Dbody_bending

%load data
load('data_3D_all.mat');
load('data_b_t_R_all.mat');
Ntrial = numel(Yraw);
Nt = 320;

%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate measures%%%%%%%%%%%%%%%%%%%%%%%%%%%

for count = 1:Ntrial
    if flag_val(count,1) == 1
        %postures
        rear(count,:) = dist_rear(Yrefined{count});
        body_elongation(count,:) = b_refined{count}(1,:);
        head_body_rotation(count,:) = abs(b_refined{count}(2,:));
        %movements
        locomotion(count,:) = dist_locomotion(t_refined{count});
        freeze(count,:) = dist_freeze(Yrefined{count});
        Drear(count,:) = [0 abs(diff(rear(count,:)))];
        general_rotation(count,:) = dist_general_rotation(R_refined{count});
        Dbody_elongation(count,:) = [0 abs(diff(b_refined{count}(1,:)))];
        Dhead_body_rotation(count,:) = [0 abs(diff(b_refined{count}(2,:)))];
    else
        %postures
        rear(count,:) = zeros(1,Nt);
        body_elongation(count,:) = zeros(1,Nt);
        head_body_rotation(count,:) = zeros(1,Nt);
        %movements
        locomotion(count,:) = zeros(1,Nt);
        freeze(count,:) = zeros(1,Nt);
        Drear(count,:) = zeros(1,Nt);
        general_rotation(count,:) = zeros(1,Nt);
        Dbody_elongation(count,:) = zeros(1,Nt);
        Dhead_body_rotation(count,:) = zeros(1,Nt);
    end
end

%good trials
ind_ok = flag_val(:,1)==1;
%dimension name
dim_name{1} = 'rear';
dim_name{2} = 'body elongation';
dim_name{3} = 'body bend';
dim_name{4} = 'locomotion';
dim_name{5} = 'freeze';
dim_name{6} = '\Delta rear';
dim_name{7} = 'body rotation'; %
dim_name{8} = '\Delta body elongation'; %
dim_name{9} = '\Delta body bend'; 

% %get un-normlized measures
% X0 = zeros(size(locomotion,1),size(locomotion,2),9);
% %postures
% X0(ind_ok,:,1) = rear(ind_ok,:);
% X0(ind_ok,:,2) = body_elongation(ind_ok,:);
% X0(ind_ok,:,3) = head_body_rotation(ind_ok,:);
% %movements    
% X0(ind_ok,:,4) = locomotion(ind_ok,:);
% X0(ind_ok,:,5) = freeze(ind_ok,:);
% X0(ind_ok,:,6) = Drear(ind_ok,:);
% X0(ind_ok,:,7) = general_rotation(ind_ok,:);
% X0(ind_ok,:,8) = Dbody_elongation(ind_ok,:);
% X0(ind_ok,:,9) = Dhead_body_rotation(ind_ok,:);
% 
% save('measures_postures_movements_unnormalized','X0','dim_name');

%%%%%%%%%%%%%%%%%%%%%%%%quantile normalizations%%%%%%%%%%%%%%%%%%%%%%%%

X = zeros(size(locomotion,1),size(locomotion,2),9);
%transform data into quantile values
Nq = 256;
Q = [1:Nq]/Nq;
%postures
X(ind_ok,:,1) = transform_data(rear(ind_ok,:),Q,false)/Nq;
X(ind_ok,:,2) = transform_data(body_elongation(ind_ok,:),Q,false)/Nq;
X(ind_ok,:,3) = transform_data(head_body_rotation(ind_ok,:),Q,false)/Nq;
%movements    
X(ind_ok,:,4) = transform_data(locomotion(ind_ok,:),Q,true)/Nq;
X(ind_ok,:,5) = transform_data(freeze(ind_ok,:),Q,true)/Nq;
X(ind_ok,:,6) = transform_data(Drear(ind_ok,:),Q,true)/Nq;
X(ind_ok,:,7) = transform_data(general_rotation(ind_ok,:),Q,true)/Nq;
X(ind_ok,:,8) = transform_data(Dbody_elongation(ind_ok,:),Q,true)/Nq;
X(ind_ok,:,9) = transform_data(Dhead_body_rotation(ind_ok,:),Q,true)/Nq;

%save data
save('measures_postures_movements','X','dim_name');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[rear] = dist_rear(X)
rear = squeeze(X(4,3,:)-X(5,3,:))';        

function[locomotion] = dist_locomotion(X)
X = X(1:2,:);
N = size(X,2);
locomotion = zeros(1,N);
for n = 2:N
    locomotion(n) = norm(X(:,n)-X(:,n-1));
end

function[freeze] = dist_freeze(X)
N = size(X,3);
freeze = zeros(1,N);
for n = 2:N
    freeze(n) = -norm(X(:,:,n)-X(:,:,n-1),'fro');
end

function[general_rotation] = dist_general_rotation(R)
N = size(R,3);
general_rotation = zeros(1,N);
for n = 2:N
    general_rotation(n) = norm(R(:,:,n)-R(:,:,n-1),'fro');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sub-fun%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[xt] = transform_data(x,Q,do_filt)
[Ntrial,Nframe] = size(x);
qval = quantile(x(:),Q);
Nq = numel(qval);
xt = zeros(Ntrial,Nframe);
for n = 1:Ntrial
    xt(n,x(n,:)<qval(1))=0;
    for m = 1:Nq-1
        xt(n,((x(n,:)>=qval(m))&(x(n,:)<qval(m+1)))) = m;
    end
    xt(n,x(n,:)==qval(end))=Nq-1;
end
%filter
if do_filt
    for n = 1:Ntrial
        x(n,:) = filtfilt([0.2 0.6 0.2],1,x(n,:));
    end
end


