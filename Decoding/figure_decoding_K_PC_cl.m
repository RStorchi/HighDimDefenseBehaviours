function[SI_inv,SIloc_inv,Npca,Npca_loc,K,Kloc] = figure_decoding_K_PC_cl()
%indexes: 
%1,2 loom vs flash, early; full set & locomotion
%3,4 sound vs flash, early; full set & locomotion
%5,6 sound vs loom, early; full set & locomotion
%7,8 loom vs flash, intermediate; full set & locomotion
%9,10 sound vs flash, intermediate; full set & locomotion
%11,12 sound vs loom, intermediate; full set & locomotion
%13,14 loom vs flash, late; full set & locomotion
%15,16 sound vs flash, late; full set & locomotion
%17,18 sound vs loom, late; full set & locomotion

%flash vs loom
[SI_inv{1},SIloc_inv{1},sSI_inv{1},sSIloc_inv{1},Npca,Npca_loc,K,Kloc] = figure_decoding_K_PC_cl_sub0(1:2);
[SI_inv{2},SIloc_inv{2},sSI_inv{2},sSIloc_inv{2}] = figure_decoding_K_PC_cl_sub0(7:8);
[SI_inv{3},SIloc_inv{3},sSI_inv{3},sSIloc_inv{3}] = figure_decoding_K_PC_cl_sub0(13:14);
Xcval2 = [0 0 0; 0 0 0; 0 0 0]; 
Xcval1 = [0 0 1; 0 0 0.8; 0 0 0.6]; 
Ycval = [0.7 0.7 0.7; 0.35 0.35 0.35; 0 0 0]; 
figure_decoding_K_PC_cl_sub1(SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,Npca,Npca_loc,Xcval1,Xcval2,Ycval);
figure_decoding_K_PC_cl_sub2(SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,K,Kloc,Xcval1,Xcval2,Ycval);

%flash vs sound
[SI_inv{1},SIloc_inv{1},sSI_inv{1},sSIloc_inv{1},Npca,Npca_loc,K,Kloc] = figure_decoding_K_PC_cl_sub0(3:4);
[SI_inv{2},SIloc_inv{2},sSI_inv{2},sSIloc_inv{2}] = figure_decoding_K_PC_cl_sub0(9:10);
[SI_inv{3},SIloc_inv{3},sSI_inv{3},sSIloc_inv{3}] = figure_decoding_K_PC_cl_sub0(15:16);
Xcval2 = [0 0 0; 0 0 0; 0 0 0]; 
Xcval1 = [1 0 0; 0.8 0 0; 0.6 0 0]; 
Ycval = [0.7 0.7 0.7; 0.35 0.35 0.35; 0 0 0];
figure_decoding_K_PC_cl_sub1(SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,Npca,Npca_loc,Xcval1,Xcval2,Ycval);
figure_decoding_K_PC_cl_sub2(SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,K,Kloc,Xcval1,Xcval2,Ycval);

%loom vs sound
[SI_inv{1},SIloc_inv{1},sSI_inv{1},sSIloc_inv{1},Npca,Npca_loc,K,Kloc] = figure_decoding_K_PC_cl_sub0(5:6);
[SI_inv{2},SIloc_inv{2},sSI_inv{2},sSIloc_inv{2}] = figure_decoding_K_PC_cl_sub0(11:12);
[SI_inv{3},SIloc_inv{3},sSI_inv{3},sSIloc_inv{3}] = figure_decoding_K_PC_cl_sub0(17:18);
Xcval2 = [0 0 1; 0 0 0.8; 0 0 0.6]; 
Xcval1 = [1 0 0; 0.8 0 0; 0.6 0 0]; 
Ycval = [0.7 0.7 0.7; 0.35 0.35 0.35; 0 0 0]; 
figure_decoding_K_PC_cl_sub1(SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,Npca,Npca_loc,Xcval1,Xcval2,Ycval);
figure_decoding_K_PC_cl_sub2(SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,K,Kloc,Xcval1,Xcval2,Ycval);

function[SI_inv,SIloc_inv,sSI_inv,sSIloc_inv,Npca,Npca_loc,K,Kloc] = figure_decoding_K_PC_cl_sub0(which_ind)
load('results_stimulus_predict_knn_review.mat');
N = numel(res); Ntrial = 344;
%
SI_inv = res(which_ind(1)).correct;
SIloc_inv = res(which_ind(2)).correct;
sSI_inv = squeeze(std(SI_inv,0,1)); 
sSIloc_inv = squeeze(std(SIloc_inv,0,1)); 
SI_inv = squeeze(mean(SI_inv,1)); 
SIloc_inv = squeeze(mean(SIloc_inv,1)); 
Npca = res(which_ind(1)).Npca;
Npca_loc = res(which_ind(2)).Npca;
K = res(which_ind(1)).NN;
Kloc = res(which_ind(2)).NN;


%figure as function of neighbours K 
function[] = figure_decoding_K_PC_cl_sub1(X,Y,sX,sY,x,xloc,Xcval1,Xcval2,Ycval)
fig = figure; set(fig,'Position',[100 100 350 300]);
h = subplot(1,1,1); hold on; 
for n = 1:numel(X)
    temp = X{n};
    stemp = sX{n};
    temp1 = []; stemp1 = [];
    for m = 1:size(temp,2)
        ind = find(temp(:,m) == max(temp(:,m)),1);
        temp1(m) = temp(ind,m);
        stemp1(m) = stemp(ind,m);
    end
    plot(log2(x),temp1,'Color',Xcval1(n,:),'LineWidth',2,'MarkerSize',18);
    errorbar(log2(x),temp1,stemp1,'.','Color',Xcval2(n,:),'LineWidth',2,'MarkerSize',18);
end
for n = 1:numel(Y)
    temp = Y{n};
    stemp = sY{n};
    temp1 = []; stemp1 = [];
    for m = 1:size(temp,2)
        ind = find(temp(:,m) == max(temp(:,m)),1);
        temp1(m) = temp(ind,m);
        stemp1(m) = stemp(ind,m);
    end
    errorbar(log2(xloc),temp1,stemp1,'.:','Color',Ycval(n,:),'LineWidth',2,'MarkerSize',18);
end
xlabel('log2(PC)','FontSize',14); ylabel('%Correct','FontSize',14); 
set(h,'FontSize',12);
xlim([log2(x(1))-1 log2(x(end))+1]);

%figure as function of PC dimensions
function[] = figure_decoding_K_PC_cl_sub2(X,Y,sX,sY,x,xloc,Xcval1,Xcval2,Ycval)
fig = figure; set(fig,'Position',[100 100 350 300]);
h = subplot(1,1,1); hold on; 
for n = 1:numel(X)
    temp = X{n}';
    stemp = sX{n}';
    temp1 = []; stemp1 = [];
    for m = 1:size(temp,2)
        ind = find(temp(:,m) == max(temp(:,m)),1);
        temp1(m) = temp(ind,m);
        stemp1(m) = stemp(ind,m);
    end
    plot(log2(x),temp1,'Color',Xcval1(n,:),'LineWidth',2,'MarkerSize',18);
    errorbar(log2(x),temp1,stemp1,'.','Color',Xcval2(n,:),'LineWidth',2,'MarkerSize',18);
end
for n = 1:numel(Y)
    temp = Y{n}';
    stemp = sY{n}';
    temp1 = []; stemp1 = [];
    for m = 1:size(temp,2)
        ind = find(temp(:,m) == max(temp(:,m)),1);
        temp1(m) = temp(ind,m);
        stemp1(m) = stemp(ind,m);
    end
    errorbar(log2(xloc),temp1,stemp1,'.:','Color',Ycval(n,:),'LineWidth',2,'MarkerSize',18);
end
xlabel('log2(#K)','FontSize',14); ylabel('%Correct','FontSize',14); 
set(h,'FontSize',12);
xlim([log2(x(1))-1 log2(x(end))+1]);
