function[] = draw_initial_conditions_cl(x,y)
%This script draw the joint probabilities of pre-stimulus (x) and response
%(y) clusters
%To calculate joint probabilities on a restricted dataset obtained with refine_clustering_cl
%use the following lines (use those also to generate Figure.7c,d,f):
% [group0, ind_ok0, ind_stim] = refine_clustering_cl(true,1);
% [group1, ind_ok1, ind_stim] = refine_clustering_cl(false,1);
% ind_ok_flash = (ind_ok0==1)&(ind_ok1==1)&(ind_stim'==0);
% ind_ok_loom = (ind_ok0==1)&(ind_ok1==1)&(ind_stim'==1);
% ind_ok_sound = (ind_ok0==1)&(ind_ok1==1)&(ind_stim'==2);
% draw_initial_conditions_cl(group0(ind_ok_flash),group1(ind_ok_flash))
% draw_initial_conditions_cl(group0(ind_ok_loom),group1(ind_ok_loom))
% draw_initial_conditions_cl(group0(ind_ok_sound),group1(ind_ok_sound))
% To generate Figure.7e use these lines:
% [group1, ind_ok1, ind_stim] = refine_clustering_cl(false,1);
% [group0, ind_ok0, ind_stim] = refine_clustering_cl(true,1,ind_ok1);

gx = unique(x); gy = unique(y); 
Nx = numel(gx); Ny = numel(gy);
N = numel(x);
Pxy  = zeros(Nx,Ny);
for n = 1:Nx
    for m = 1:Ny
        Pxy(n,m) = sum((x==gx(n))&(y==gy(m)));
    end
end
Pxy = Pxy/N;
%draw figure
%fig1
figure; hold on;
Pxy_line = 5*Pxy/max(Pxy(:));
plot(ones(1,Nx),gx,'O','MarkerSize',12,'LineWidth',5);
plot(2*ones(1,Ny),gy,'O','MarkerSize',12,'LineWidth',5);
for n = 1:Nx
    for m = 1:Ny
        if Pxy_line(n,m)>0
           line([1 2],[gx(n) gy(m)],'LineWidth',Pxy_line(n,m),'Color','k');
        end
    end
end
xlim([0 3]); ylim([0 max(Nx,Ny)+1]);
%fig2 - legend
figure; hold on;
xlim([0 3]); ylim([0 max(Nx,Ny)+1]);
min_val = min(Pxy(Pxy>0));
mean_val = mean(Pxy(Pxy>0));
max_val = max(Pxy(:));
line([1.3 1.6],[3 3],'LineWidth',5*min_val/max_val,'Color','k');
line([1.3 1.6],[4 4],'LineWidth',5*mean_val/max_val,'Color','k');
line([1.3 1.6],[5 5],'LineWidth',5,'Color','k');
title([num2str(min_val) ' ' num2str(mean_val) ' ' num2str(max_val)]);