function[P,Xest0,X,err] = calibrate_cameras_cl()
%estimate the camera matrices P by using the direct linear method  
%OUTOUT: 
%P: camera matrices (one cell per camera)
%X0est: reconstructed 3D coordinates 
%X: real world coordinates


Ny = 1280; Nx = 1024;

%load world coordinates
load('world_coordinates');
Np = size(X,2);

%load all coordinates
which_camera = [1:4]; Nc = length(which_camera); 
which_height = [0:4]; Nh = length(which_height);
x = cell(1,Nc);
for n = 1:Nc
    for m = 1:Nh
        load(['camera' num2str(which_camera(n)) '_lego_h' num2str(which_height(m)) '_data']);
        temp = [idx_final; idy_final; ones(1,length(idx_final))];
        x{n} = [x{n} temp];
    end
end

%convert each multiview point into a cell
x_multi = cell(1,Np);
for n = 1:Np
    for m = 1:Nc
        x_multi{n}(:,m) = x{m}(:,n);
    end
end



%initial estimate of camera matrices
P = cell(1,Nc);
for n = 1:2
    P{n} = DLT_simple_cl(x{n},X);
end
for n = 3:4
    P{n} = DLT_simple_cl(x{n},X1);
end

%use LS triangulation (ls_triangulate) to get an  estimate of the 3D
%coordinates
for n = 1:Np
    Xest0(:,n) = ls_triangulate(x_multi{n},P);
end

%calculate error
err = sqrt(sum((Xest0-X).^2));

%fig initial estimate
fig0 = figure; 
h0 = subplot(1,1,1); hold on;
plot3(X(1,:),X(2,:),X(3,:),'.k','MarkerSize',10);
plot3(Xest0(1,:),Xest0(2,:),Xest0(3,:),'.r','MarkerSize',10);

%save camera matrices
%save([filepath 'Pcal'],'P','X','Xest0');

