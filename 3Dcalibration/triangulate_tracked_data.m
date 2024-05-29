function[Y] = triangulate_tracked_data()
%This script perform triangulation of data from 4 camera views
%OUTPUT: 
%Y: is an Nbp x 3 x Nframe matrix. Nbp is the number of body landmarks, 3
%correspond to the x,y,z coordinates of each landmark, Nframe is the number of frames 

%input here filepath and file names of the deeplabcut tracking files
filepath = 'sample_data';
filename_camera{1} = 'camera1_mouse1.csv';
filename_camera{2} = 'camera2_mouse1.csv';
filename_camera{3} = 'camera3_mouse1.csv';
filename_camera{4} = 'camera4_mouse1.csv';

%load calibration file
load('Pcal.mat');

%initial parameters
Ncam = 4;
THlike = 0.8;
bp = [];
like = [];

%%%%%%%%%%%%%%%%%%%load DLC coordinates%%%%%%%%%%%%%%%%%%%%
for n = 1:Ncam
    temp = csvread([filepath '\' filename_camera{n}],3,1);
    [Nframe,Nbp] = size(temp);
    Nbp = Nbp/3;
    for m = 1:Nbp
        bp{n,m} = temp(:,3*(m-1)+1:3*(m-1)+2); %bp coordinates
        like{n,m} = temp(:,3*m); %bp likelihood
    end
end

%THIS FUNCTION: facilitates triangulation by pooling coordinate on each
%frame. Specifically it:
%1) Reformat the dimensions: 
%   FROM bp - cell(Ncam,Nbp), each cell containing Nframe-by-2 points
%   TO bp_t - cell(Nbp,Nframe), each cell containing 3-by-Ncam points 
bp_t = transform_coord(bp);

%THIS FUNCTION: trasnforms like: 
%   FROM like - cell(Ncam,Nbp), each cell containing Nframe-by-1 points
%   TO like_t - cell(Nbp,Nframe), each cell containing 1-by-Ncam points
like_t = transform_like(like);

%THIS FUNCTION: initial raw estimates in homogeneous coordinates. Outputs:
%   Y - cell(1,Nbp), each cell containing 4-by-Nframe points (homogeneous)
%   miss - logic Nbp-by-Nframe matrix, flagging false when <2 cameras have
%          reliable landmarks (reliability defined by like_t)
[X, miss] = reconstruct_3d(bp_t,P,like_t,THlike);

%Convert the data
Y = transform_X3D_format(X);
for n = 1:Nbp
    Y(n,:,miss(n,:)==1) = NaN;
end

%save data
save('initial_3D','Y');

%%%%%%%%%%%%%%%%%%SUB FUCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[x] = transform_coord(x0)
[Ncam,Nbp] = size(x0);
Nframe = size(x0{1},1);
x = cell(Nbp,Nframe);
for n = 1:Nbp
    for m = 1:Nframe
        x{n,m} = ones(3, Ncam);
        for p = 1:Ncam
            x{n,m}(1,p) = x0{p,n}(m,2); 
            x{n,m}(2,p) = x0{p,n}(m,1); 
        end
    end
end

function[x] = transform_like(x0)
[Ncam,Nbp] = size(x0);
Nframe = size(x0{1},1);
x = cell(Nbp,Nframe);
for n = 1:Nbp
    for m = 1:Nframe
        x{n,m} = ones(1, Ncam);
        for p = 1:Ncam
            x{n,m}(1,p) = x0{p,n}(m,1); 
        end
    end
end

function[X,miss] = reconstruct_3d(x,P,like,TH)
[Nbp,Nframe] = size(x); %size of both x and like
X = cell(1,Nbp);
miss = false(Nbp,Nframe);
for n = 1:Nbp
    X{n} = zeros(4,Nframe);
    for m = 1:Nframe
        ind = find(like{n,m}>TH); %camera indexes with good likelihood
        if length(ind) >= 2
            X{n}(:,m) = ls_triangulate(x{n,m}(:,ind), P(ind));
        else
            miss(n,m) = true;
        end
        disp(sprintf(['frame: %s body_part: %s'],num2str(m),num2str(n)));
    end
end