
min_clus = 2; %minimum number of clusters
max_clus = 10; %maximum number of clusters
groups = {};
% the six time interval sequences
sequence_1 = {141:142, 143:144, 145:146, 147:148, 149:150, 151:152, 153:154, 155:156, 157:158, 159:160, 161:162, 163:164, 165:166, 167:168, 169:170, 171:172, 173:174, 175:176, 177:178, 179:180};
sequence_2 = {139:141,142:144, 145:147, 148:150, 151:153, 154:156, 157:159, 160:162, 163:165, 166:168, 169:171, 172:174, 175:177, 178:180};
sequence_3 = {141:145, 146:150, 151:155, 156:160, 161:165, 166:170, 171:175, 176:180};
sequence_4 = {139:144, 145:150, 151:156, 157:162, 163:168, 169:174, 175:180};
sequence_5 = {141:150, 151:160, 161:170, 171:180};
sequence_6 = {136:150, 151:165, 166:180};
for i=min_clus:max_clus
    [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(sequence_1, i);
    writematrix(group(:,6:end),['training_set' num2str(i) '_seq1.csv']);
    groups{1,i-1} = group;
end
for i=min_clus:max_clus
    [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(sequence_2, i);
    writematrix(group(:,5:end),['training_set' num2str(i) '_seq2.csv']);
    groups{2,i-1} = group;
end
for i=min_clus:max_clus
    [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(sequence_3, i);
    writematrix(group(:,3:end),['training_set' num2str(i) '_seq3.csv']);
    groups{3,i-1} = group;
end
for i=min_clus:max_clus
    [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(sequence_4, i);
    writematrix(group(:,3:end),['training_set' num2str(i) '_seq4.csv']);
    groups{4,i-1} = group;
end
for i=min_clus:max_clus
    [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(sequence_5, i);
    writematrix(group(:,2:end),['training_set' num2str(i) '_seq5.csv']);
    groups{5,i-1} = group;
end
for i=min_clus:max_clus
    [group,Ngroup,C,Y] = cluster_movements_multi_timestep_review_v3(sequence_6, i);
    writematrix(group(:,2:end),['training_set' num2str(i) '_seq6.csv']);
    groups{6,i-1} = group;
end
save('final_workspace.mat')