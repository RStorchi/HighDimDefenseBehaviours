function[b_opt, dist_opt, R_opt, t_opt] = reconstruct_procrustes_fminunc_cl(X,mu,U,vars,alpha_reg)
%Minimise cost function in eq.3 robust reconstruction based on shape constraints
%INPUT: 
%X: 3D coordinates
%mu: mean pose
%U: eigenposes
%vars: eigenvalues associated with eigenposes
%alpha_reg: regularization parameter
%OUTPUT:
%b_opt: shape parameters
%dist_opt: distance from mean pose (after alignment with procrustes)
%R_opt: 3D rotation matrix
%t_opt: 3D translation matrix


%%%%%%%%%%%%INITIALIZE SHAPE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%
%first I perform procrustes superimposition to each shape individually
[Np,~,N] = size(X);
Npca = numel(vars);
b = zeros(Npca,N);
dist = zeros(1,N);
mu1 = reshape(mu,3*Np,1);
for n = 1:N
    %calculate proctrustes
    [~, Z, tr] = procrustes(mu,squeeze(X(:,:,n)),'Reflection',false, 'Scaling',false);
    %center
    Z = Z-repmat(mean(Z),Np,1);
    %generate 1d 
    Z1 = reshape(Z,3*Np,1);
    %calculate shapes
    b(:,n) = U'*(Z1-mu1);
    %calculate dist
    dist(n) = norm(mu-Z);
end

%%%%%%%%%%%%MINIIMIZE COST FUNCTION EQ3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%estimate the time series shapes and translations
options = optimoptions('fminunc','Display','none');
%
b_opt = zeros(Npca,N);
dist_opt = zeros(1,N);
R_opt = zeros(3,3,N);
t_opt = zeros(3,N);
for n = 1:N
    %init b values with initial estimate above
    %minimise
    b_opt_temp = fminunc(@(b_opt_temp) test_robust_pose_estimation_v1_fun(b_opt_temp, U, mu, squeeze(X(:,:,n)), vars, alpha_reg), b(:,n), options);
    b_opt(:,n) = b_opt_temp;
    %get error
    [~,dist_opt(n), R_opt(:,:,n), t_opt(:,n)] = test_robust_pose_estimation_v1_fun(b_opt_temp, U, mu, squeeze(X(:,:,n)), vars, alpha_reg);  
end
%
d_err = abs(dist);
d_err_opt = abs(dist_opt);
std_b_opt = std(b_opt');
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[C, dist, R, t] = test_robust_pose_estimation_v1_fun(b, U, X, Y, vars, alpha_reg)
Np = size(X,1);
%
X1 = reshape(X,3*Np,1);
%
X1 = X1 + U*b;
%
X = reshape(X1,Np,3);
%
[~, Z, tr] = procrustes(Y,X,'Reflection',false, 'Scaling',false);
%
dist = norm(Y-Z);
%
reg = alpha_reg*sum((b.^2)./vars);
%
C = dist+reg;
%
R = tr.T;
t = tr.c(1,:)';


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[] = show_3d(x,fig)
figure(fig); 
Np = size(x,1);
for n = 1:Np
    plot3(x(n,1),x(n,2),x(n,3),'.b','MarkerSize',20);
end

