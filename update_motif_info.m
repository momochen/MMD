function [motifs new_id] = update_motif_info(motifs,merged_start,merged_id)
% Update the motif. Replace the motif from the merged_id with the lowest
% id in the group as the new id

sorted_id_list = sort(unique(merged_id));
new_id = sorted_id_list(1,1);

for i = 1:length(merged_start)
    
    motifs(find(motifs(:,2)==merged_start(i,1)),1) = new_id;
    
end