% all_starts: all starts pos that needs to be merged
% motif: the lists that contain all the motifs
function new_motif = merge_sub_motif(all_starts,win_size,motif)
% This function merge the sub-motifs and return the merged motifs

%all_starts

% Sort the startings points
sorted_starts = sort(all_starts);

% For the sorted list, find the endings for each
sorted_ends = [];
for i = 1:length(sorted_starts)
   
    row_num = find(motif(:,2)==sorted_starts(i,1));
    crs_motif = motif(row_num,:);
    crs_end = crs_motif(1,2)+win_size;
    sorted_ends = [sorted_ends;crs_end];
    
end
%sorted_ends

% for each of the starting points

new_motif = [sorted_starts(1,1) sorted_ends(1,1)];
new_merged_motif_list = [new_motif];

for i = 2:length(sorted_starts)

   % If there is overlap, then merge and starts with the second ends
   if sorted_starts(i,1)<new_motif(1,2)
      new_motif(1,2) = sorted_ends(i,1);
   else % Otherwise, this is the end of the merged motif, starts a new one
       new_merged_motif_list = [new_merged_motif_list;[new_motif(1,1) new_motif(1,2)]];
       new_motif = [sorted_starts(i,1) sorted_ends(i,1)];
   end
    
end


