function[correct,NN,Npca] = decode_cl(epoch, which_stim, which_dim, Nrep)
%INPUT:
%epoch: for 0-1s type 'early'; for 1-2s type 'intermediate';for 2-3s type
%'late'; 0-2s type 'all'
%which_stim: for flash vs loom = [0 1]; for flash vs sound = [0 2]; for loom vs sound = [1 2];
%which_dim: type 4 for locomotion or 1:9 for the full set
%OUTPUT
%correct: matrix of dimensions Nrep x numel(NN) x numel(Npca)
%NN: nearest neighbour k (can be changed from line 28)
%Npca: number of PCs (can be changed from lines 30-32)

%load data
load('../3Dreconstruction/measures_postures_movements_10bin','X'); 
load('../3Dreconstruction/data_3D_all.mat','flag_val');

%epoch selection
if strcmp(epoch,'early')
    indX = 151:165;   
elseif strcmp(epoch,'intermediate')
    indX = 166:180; 
elseif strcmp(epoch,'late') 
    indX = 181:195; 
elseif strcmp(epoch,'all')
    indX = 151:180; 
else
   disp('epoch not recognized, see help');      
end

%NNs & PCs
NN = 2.^[0:7];
if numel(which_dim)>1
    Npca = [1 2 4 8 16 32 64 ];
    %Npca = [2:2:100];
else
    Npca = [1 2 4 8 15];
end

%select good trials
ind_good = find(flag_val(:,1) == 1);
X = X(ind_good,:,which_dim);
flag_val = flag_val(ind_good,:);
stim0 = flag_val(:,5);
stim = stim0;
stim((stim0==1)|(stim0==2)) = 1;
stim((stim0==3)|(stim0==4)) = 2;

%generate dependent variable Y
if numel(which_stim) == 2
    indS = find((stim==which_stim(1))|(stim==which_stim(2)));
    X = X(indS,:,:);
    stim = stim(indS);
    Y = stim; 
    Y(stim==which_stim(1)) = 0;
    Y(stim==which_stim(2)) = 1;
else
    Y = stim;
end
%
[N,~,Ndim] = size(X);
%
Nind = numel(indX);
merge = zeros(N, Nind*Ndim);
for n = 1:N
    for m = 1:Ndim
        merge(n,(m-1)*Nind+1:m*Nind) = squeeze(X(n,indX,m));
    end
end

%PC reduction
[coeff, score, vars] = pca(merge);

%decode
for n = 1:numel(NN)
    for m = 1:numel(Npca)
        X = score(:, 1:Npca(m));
        for nr = 1:Nrep
            Mdl = fitcknn(X,Y,'NumNeighbors',NN(n),'DistanceWeight','inverse','CrossVal','on','KFold',10);
            correct(nr,n,m) = sum(kfoldPredict(Mdl) == Y)/N;  
        end
        disp(sprintf('K: %s Npca: %s',num2str(NN(n)),num2str(Npca(m))));
    end
end





