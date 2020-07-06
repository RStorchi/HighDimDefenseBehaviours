function[pval] = figure_specificity_K_PC_cl()
%Plot the SI index for locomotion and full set as function of Principal
%Components and K. Also shows errorbar for SI associated with all the three
%stimuli. Returns the p-value obtained for full set vs locomotion by 
%using a sign-test

Kneigh = [1 2 4 8 16 32 64];
%flash vs loom
[SI_inv{1}, SIloc_inv{1}] = figure_specificity_K_PC_cl_sub0(Kneigh,1,0); %early
[SI_inv{2}, SIloc_inv{2}] = figure_specificity_K_PC_cl_sub0(Kneigh,4,0); %intermediate
[SI_inv{3}, SIloc_inv{3}] = figure_specificity_K_PC_cl_sub0(Kneigh,7,0); %late
Xcval1 = [0 0 1; 0 0 0.8; 0 0 0.6]; 
Xcval2 = [0 0 0; 0 0 0; 0 0 0]; 
Ycval = [0.7 0.7 0.7; 0.35 0.35 0.35; 0 0 0]; 
figure_specificity_K_PC_review_v1_sub1(SI_inv,SIloc_inv,Xcval1,Xcval2,Ycval);
figure_specificity_K_PC_review_v1_sub2(SI_inv,SIloc_inv,Xcval1,Xcval2,Ycval);

%flash vs sound
[SI_inv{1}, SIloc_inv{1}] = figure_specificity_K_PC_cl_sub0(Kneigh,2,0); %early
[SI_inv{2}, SIloc_inv{2}] = figure_specificity_K_PC_cl_sub0(Kneigh,5,0); %intermediate
[SI_inv{3}, SIloc_inv{3}] = figure_specificity_K_PC_cl_sub0(Kneigh,8,0); %late
Xcval1 = [1 0 0; 0.8 0 0; 0.6 0 0]; 
Xcval2 = [0 0 0; 0 0 0; 0 0 0];
Ycval = [0.7 0.7 0.7; 0.35 0.35 0.35; 0 0 0]; 
figure_specificity_K_PC_review_v1_sub1(SI_inv,SIloc_inv,Xcval1,Xcval2,Ycval);
figure_specificity_K_PC_review_v1_sub2(SI_inv,SIloc_inv,Xcval1,Xcval2,Ycval);

%loom vs sound
[SI_inv{1}, SIloc_inv{1}] = figure_specificity_K_PC_cl_sub0(Kneigh,3,0); %early
[SI_inv{2}, SIloc_inv{2}] = figure_specificity_K_PC_cl_sub0(Kneigh,6,0); %intermediate
[SI_inv{3}, SIloc_inv{3}] = figure_specificity_K_PC_cl_sub0(Kneigh,9,0); %late
Xcval1 = [1 0 0; 0.8 0 0; 0.6 0 0]; 
Xcval2 = [0 0 1; 0 0 0.8; 0 0 0.6]; 
Ycval = [0.7 0.7 0.7; 0.35 0.35 0.35; 0 0 0]; 
figure_specificity_K_PC_review_v1_sub1(SI_inv,SIloc_inv,Xcval1,Xcval2,Ycval);
figure_specificity_K_PC_review_v1_sub2(SI_inv,SIloc_inv,Xcval1,Xcval2,Ycval);

%early+intermediate epochs for all stimuli
[SI_inv, SIloc_inv] = figure_specificity_K_PC_cl_sub0(1,1,1);
%
[ind_best_row,ind_best_col] = find(SI_inv == max(SI_inv(:)));
load(['SIinv_review_K' num2str(Kneigh(ind_best_row)) '_v1_early_inter'],'spec');
spec = spec(:,ind_best_col);
%
[ind_best_row,ind_best_col] = find(SIloc_inv == max(SIloc_inv(:)));
load(['SIinv_review_K' num2str(Kneigh(ind_best_row)) '_v1_early_inter'],'spec_loc');
spec_loc = spec_loc(:,ind_best_col);
%
N = numel(spec);
fig = figure; set(fig,'Position',[100 100 250 300]);
h = subplot(1,1,1); hold on; 
xval = 13; 
bar(xval, mean(spec),'BarWidth',0.5,'FaceColor',0.85*[1 0 1],'EdgeColor','k','LineWidth',2);
errorbar(xval,mean(spec),std(spec)/sqrt(N),'.','LineWidth',2,'Color',[1 1 1]*0.6);
bar(xval,mean(spec_loc),'BarWidth',0.5,'FaceColor',0.8*[1 1 1],'EdgeColor','k','LineWidth',2);
errorbar(xval,mean(spec_loc),std(spec_loc)/sqrt(N),'.','LineWidth',2,'Color',[1 1 1]*0.6);
xlim([-1 1]+xval); ylim([0.5 1]);
%stats
pval = signtest(spec-spec_loc);
%title
title(['pval = ' num2str(pval)]);

%load results
function[SI,SIloc] = figure_specificity_K_PC_cl_sub0(Kneigh,ind,early_inter)
Nk = numel(Kneigh); 
for n = 1:Nk
    if early_inter
        load(['SIinv_review_K' num2str(Kneigh(n)) '_v1_early_inter.mat'],'spec','spec_loc');
    else
        load(['SIinv_review_K' num2str(Kneigh(n)) '_v1.mat'],'spec','spec_loc');    
    end
    SI(n,:) = mean(mean(spec(:,:,ind),1),3); 
    SIloc(n,:) = mean(mean(spec_loc(:,:,ind),1),3); 
end

%figure as function of neighbourhood 
function[] = figure_specificity_K_PC_review_v1_sub1(X,Y,Xcval1,Xcval2,Ycval)
Kneigh = [1 2 4 8 16 32 64 ];
fig = figure; set(fig,'Position',[100 100 350 300]);
h = subplot(1,1,1); hold on; 
for n = 1:numel(X)
    plot(log2(Kneigh),max(X{n}'),'Color',Xcval1(n,:),'LineWidth',2);
    plot(log2(Kneigh),max(X{n}'),'.','Color',Xcval2(n,:),'LineWidth',2,'MarkerSize',18);
end
for n = 1:numel(Y)
    plot(log2(Kneigh),max(Y{n}'),'.:','Color',Ycval(n,:),'LineWidth',2,'MarkerSize',18);
end
xlabel('log2(K)','FontSize',14); ylabel('SI','FontSize',14); 
set(h,'FontSize',12);
xlim([-1 log2(Kneigh(end))+1]);

%figure as function of dimensions
function[] = figure_specificity_K_PC_review_v1_sub2(X,Y,Xcval1,Xcval2,Ycval)
PCs = [1:15];
fig = figure; set(fig,'Position',[100 100 350 300]);
h = subplot(1,1,1); hold on; 
for n = 1:numel(X)
    plot(PCs,max(X{n}),'Color',Xcval1(n,:),'LineWidth',2);
    plot(PCs,max(X{n}),'.','Color',Xcval2(n,:),'LineWidth',2,'MarkerSize',18);
end
for n = 1:numel(Y)
    plot(PCs,max(Y{n}),'.:','Color',Ycval(n,:),'LineWidth',2,'MarkerSize',18);
end
xlabel('PCs','FontSize',14); ylabel('SI','FontSize',14); 
set(h,'FontSize',12);
xlim([0 PCs(end)+1]);



