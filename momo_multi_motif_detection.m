function [multi_dim_motif] = momo_multi_motif_detection(motifs_collection,motif_incidence_table,alpha)
% This function aggregate multi-dimensional motif according to the
% incidence table
multi_dim_motif = [];
size_incidence_table = size(motif_incidence_table);

for i = 1:size_incidence_table(1,1)
    
    temp_multi_id = find(motif_incidence_table(i,:)>=alpha);
    for j = 1:length(temp_multi_id)
       motif_incidence_table(temp_multi_id(1,j),i) = 0;
    end
    temp_multi_id = [i temp_multi_id];
    % Show the multi-dimensional motif
    if length(temp_multi_id)>1
        multi_dim_motif = temp_multi_id;
        multi_dim_motif
    end
    temp_multi_id = [];
end