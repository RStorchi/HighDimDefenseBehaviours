function[pval,es] = example_locomotion_elongation_cl()
%for sound vs loom
load('../3Dreconstruction/data_3D_all.mat','flag_val');
load('../3Dreconstruction/measures_postures_movements.mat');
%load('C:\Users\mqbphrs3\Documents\offline_stuff\Behaviour\visual\movements_quantified.mat');
X = X(flag_val(:,1)==1,:,:);
flag_val = flag_val(flag_val(:,1)==1,:);
ind_stim_vec0 = flag_val(:,5);
ind_stim_vec = ind_stim_vec0;
ind_stim_vec((ind_stim_vec0==1)|(ind_stim_vec0==2)) = 1;
ind_stim_vec((ind_stim_vec0==3)|(ind_stim_vec0==4)) = 2;
%
which_stim = [1 2];
indX0 = 166:180;
mean_loc = mean(X(:,indX0,4),2);
ind1 = find(ind_stim_vec==which_stim(1));
ind2 = find(ind_stim_vec==which_stim(2));
ind_freeze = find((mean_loc>=0)&(mean_loc<0.1));
ind1 = intersect(ind1,ind_freeze);
ind2 = intersect(ind2,ind_freeze);
%locomotion & body elongation
which_dim = [4 2]; Ndim = numel(which_dim);
dim_name = {'Locomotion','Body Elongation'};
fig = figure; 
set(fig,'Position',[100 100 600 250]);
for n = 1:Ndim
    h(n) = subplot(1,2,n); hold on;
    indX = 151:196;
    tval = (indX-150)/15;
    temp{1} = squeeze(X(ind1,indX,which_dim(n)));
    temp{2} = squeeze(X(ind2,indX,which_dim(n)));
    errorbar(tval,mean(temp{1}),std(temp{1})/sqrt(numel(ind1)),'b','LineWidth',2);
    errorbar(tval,mean(temp{2}),std(temp{2})/sqrt(numel(ind2)),'r','LineWidth',2);
    xlim([tval(1)-0.1 tval(end)+0.1]); ylim([-0.1 1.1]);
    title(dim_name{n});
    xlabel('Time(s)','FontSize',12);
    %test difference
    [pval(:,:,n),es(:,:,n)] = test_mean_difference_cl(temp);
end





