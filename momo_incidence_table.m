function [motif_collection weighted_incidence_table] = momo_incidence_table(motifs)
% This function reassigned id to each motif and return a incidence table
% between motifs of different dimension 

% Firstly read all motifs and reassign them new ids
motif_collection = [];
[sb2 num_dimension] = size(motifs);
new_id = 0;
current_id = 0;

% For each dimension
for i = 1:num_dimension
    
    [single_dim_motif_row single_dim_motif_col] = size(motifs(i).motif);
    current_dim_motif = motifs(i).motif;
    
    for j = 1:single_dim_motif_row
        
        if current_id == current_dim_motif(j,1)
            current_dim_motif(j,1) = new_id;
        else
            new_id = new_id + 1;
            current_id = current_dim_motif(j,1);
            current_dim_motif(j,1) = new_id;
            
        end
        motif_collection = [motif_collection;[i current_dim_motif(j,:)]];
        
    end
    
    %new_id = new_id + 1;
    current_id = 0;
    
end

% Create the incidence table
incidence_table_length = length(unique(motif_collection(:,2)));

weighted_incidence_table = zeros(incidence_table_length,incidence_table_length);
size_motif_collection = size(motif_collection);

% The rows are the concerning motifs, the cols are the corresponding motifs
% After that, the matrix will be weighted by dividing each row of the
% occurance of the motif corresponds to that row number

% Buggy:... One overlapping between multi-dimensional motif should only be
% counting once
overlapped_id_list = [];

for j = 1:size_motif_collection(1,1)
    
    for k = 1:size_motif_collection(1,1)
       
        % Don't compare the motif in the same dimension
        if motif_collection(j,1)~=motif_collection(k,1)
            % Check whether two motifs from different dimension overlaps;if
            % they do, then add 1 to the incidence table
            if (motif_collection(j,3)>=motif_collection(k,3) && motif_collection(j,3)<=motif_collection(k,4)) ...
                    || (motif_collection(j,4)>=motif_collection(k,3) && motif_collection(j,4)<=motif_collection(k,4))
                if ~ismember(motif_collection(k,1),overlapped_id_list)
                    weighted_incidence_table(motif_collection(j,2),motif_collection(k,2)) = weighted_incidence_table(motif_collection(j,2),motif_collection(k,2)) + 1;
                    overlapped_id_list = [overlapped_id_list;motif_collection(k,1)];
                end
            end
        end
    end
    overlapped_id_list = [];
    
end

% Counted the occurance of each motif again
new_motif_occ = [];

for i = 1:new_id
   
    weighted_incidence_table(i,:) = weighted_incidence_table(i,:)./sum(motif_collection(:,2)==i);
    
end



