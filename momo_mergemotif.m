% motif:            [id pos end occ_num]
% motif_occ:        [id occ_num]
% distance:         Integer that represents distance criteria
% win_size:         Window size of the original motifs
function [merged_motif merged_motif_occ] = momo_mergemotif(motif,motif_occ,distance,win_size)

% Merge the raw motif into motif of different length according to the
% occurance and position of the motif 

% Using a motif-collision table, where each cell represents the
% co-occurance of two motif with different id, i.e. two different motif

% Define the co-occurance as the occurance of two motif within a distance
% of 'd'
% This is a recursive function, it would continue to find overlapping
% motifs until all different motifs are apart from each other


d = distance;
%merged_indicator = 1;

% while the distance between different motifs are less than $distance
% and overalapped motifs are find in the PREVIOUS step(just in case the
% infinite loop will occur), keep finding the motifs
%while merged_indicator == 1 
    
%    merged_indicator = 0;
    
    % Find the overlapped motif based on the motif collision table
    size_distinct_motif = size(motif_occ);
    motif_collision_table = zeros(size_distinct_motif(1,1),size_distinct_motif(1,1));

    size_motif = size(motif);

    for i = 1:size_motif(1,1)

        for j = i+1:size_motif(1,1)

            % If two motif are within a distance of d, then add to the table
            % This parameter must be carefully selected!!!!
            % the distance between these two motifs must be within a small
            % range otherwise if the difference is too big, it doesn't make
            % much sense to merge them togeterh
            if abs(motif(i,2) - motif(j,2))<=d
                if motif(i,1) ~= motif(j,1)
                    motif_collision_table(motif(i,1),motif(j,1)) = motif_collision_table(motif(i,1),motif(j,1)) + 1;
                end
            end

        end

    end
    
    % Find the cell whose collision is the same as its own occurance,
    % then treat these two as a single motif,i.e. merge them together
    size_motif_occ = size(motif_occ);
    % Fold the matrix to its upper triangular
    folded_mat = motif_collision_table + motif_collision_table';
    % zeros out the lower-triangular
    upper_folded_mat = folded_mat.*triu(ones(size(folded_mat)));
    % new motif matrix,storing the new motifs
    new_motif_matrix = [];
    
    % Iterate each point and take a note on which motif have been clustered
    % to new motifs; if a set of motif have been clustered to new, then
    % continue looking for the same motif pattern in the rest of the time
    % series, until one reaches the end of the sequence
    % Then update the motif set by removing the clustered ones, continue finding
    % the new motifs from pos+1;keep looking until each one reach the end of the 
    % raw sequence
    
    % Arrange an array that contains all starts and ends
    % This is the search space,which is a col vector
    search_space = sort([motif(:,2);motif(:,2)+win_size-1]);
    
    % Create a new motif list
    new_merged_motifs = [];
    
    % New motif
    new_merged_motifs_id = 1;
    
    % Initialize where to starts from 
    % The first one should be the motif at the first index
    %current_top = min(search_space(:,1));
    
    
    fprintf('Recursion begins\r\n');
    % While there are still elements that had not been re-motifed, do
    % the motif detection algorithm
    while length(search_space)>1
        
        %fprintf('The current top is \r\n');
        %current_top
        %fprintf('Current space is %d\r\n',length(search_space));
        
        element_pos = min(search_space(:,1));
        %element_pos
        % Find the motif with the starting point of top_start
        motif_ind = find(motif(:,2)==element_pos);
        %length(motif_ind)
        large_motif_pos_range = [];
        
        
        if length(motif_ind)>0
        
            one_motif = motif(motif_ind(1,1),:);
            
            %{
            one_motif_occ_num = find(motif_occ(:,1)==one_motif(1,1));
            co_occur_in_row = find(upper_folded_mat(:,one_motif(1,1))==motif_occ(one_motif_occ_num,2));
            co_occur_in_col = find(upper_folded_mat(one_motif(1,1),:)==motif_occ(one_motif_occ_num,2));
            
            unique_id = unique([co_occur_in_row' co_occur_in_col]);
            %}
            
            one_motif_start = one_motif(1,2);
            one_motif_end = one_motif_start + win_size - 1;
            
            
            % The number of occurance for this motif
            one_motif_occ = motif_occ(find(motif_occ(:,1)==one_motif(1,1)),2);
           
            % The ids that had already been counted
            counted_cls = [one_motif(1,1)];
            
            % Call the recursive function 
            % The function will return the list of start and end, then we
            % know what motifs can be removed and what is the new start
            %fprintf('Recursion Start\r\n');
            [merged_start merged_end merged_id]=find_merged_motif_recursion(motif,motif_occ,search_space,one_motif_start,one_motif_end,win_size,one_motif_occ,upper_folded_mat,counted_cls);
            %fprintf('Recursion End\r\n');
            % For the list of starts and ends, re-motif them as new motif
            % and add them to the new motif list
            % According to the id list, move all the motif within the id
            % list to a new temporary list and merge the motif according to
            % the overlapping
            
            % Find the min in merged start and max in merged end as the
            % new motif
            %merged_start
            new_motif_list = merge_sub_motif(merged_start,win_size,motif);
            %fprintf('Merged Once\r\n');
            %new_motif_list
            % For all the merged motifs, replace them with the lowest id 
            [motif new_id]= update_motif_info(motif,merged_start,merged_id);
            %fprintf('Updated new motif\r\n');
            
            % For each row, store it as a new motif, count the occurance of
            % the motif by counting the number of rows
            new_motif_list_size = size(new_motif_list);
            new_merged_motifs = [new_merged_motifs;[ones(new_motif_list_size(1,1),1).*new_id new_motif_list]];
           
            % Remove the lists in starts and ends from the motifs list and
            % find all new positions of this kind of motif
            %fprintf('before\r\n');
            %size(search_space)
            
            % Trial
            %search_space = setdiff(search_space,[merged_start;merged_end]);
            %size(search_space);
            
            search_space = array_diff(search_space,sort([merged_start;merged_end]));
            %fprintf('Search space\r\n');
            %size(search_space)
            %fprintf('Removed items\r\n');
            %ismember([merged_start;merged_end],search_space)
            
            %size(search_space);
            
            %fprintf('after\r\n');
            %size(search_space)
        else
            %fprintf('Not found\r\n');
            
        end
            
        
    end
    
    %new_merged_motifs
    fprintf('All motif merged\r\n');
    
    merged_motif = new_merged_motifs;
    
    % Sort the motif according to their id, put the motif with the same id
    % together
    sorted_id_list = sort(unique(merged_motif(:,1)));
    ordered_motif = [];
    for j = 1:length(sorted_id_list)
       
        ordered_motif = [ordered_motif;merged_motif(find(merged_motif(:,1)==sorted_id_list(j,1)),:)];
        
    end
    
    merged_motif = ordered_motif;
    
    merged_motif_occ = [];
    
    %{
    for k = 1:size_motif_occ(1,1)
        
       % for each of the id-unique motif
        
       % find the collision motifs
        % This is a col vector
        co_occur_in_row = find(upper_folded_mat(:,motif_occ(k,1))==motif_occ(k,2));
        % This is a row vector
        co_occur_in_col = find(upper_folded_mat(motif_occ(k,1),:)==motif_occ(k,2));
        
       % The row or col num is the index of the motif that to be merged
        % Find the unique ids
        unique_id = unique([co_occur_in_row' co_occur_in_col]);
        % These ids are to be replaced by the current id
        
        membership_motif = ismember(motif(:,1),unique_id);
        bool_membership_motif = logical(membership_motif);
        motif_to_be_merged = motif(bool_membership_motif,:);
        
        for n = 1:size_motif(1,1)
           if  
        end
            
            
        if motif_occ(k,2) == co 
            
        end
        
            
    end
    %}
    
    % Form the new motif, update the motif[] and motif_occ[]
    
    % If there is merged motif found, then change the indicator to 1,
    % making sure in the next round the the motif-merge mechanism keeps
    % running;otherwise, keeps the indicator to 0 such that it gets out of
    % the looping
    % 
    
    
%end  

    
