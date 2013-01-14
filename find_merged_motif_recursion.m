% motifs [id pos_start pos_end occur]
% pos_list the list of all starting and ending points
% org_start starting point for the range
% org_end   ending point for the range
% num_occur the number of occur as the criteria for merging motifs
% collision_table the collision table of each motifs
% counted_cls: the collisioned id that had already been counted
function [merged_start merged_end merged_id]=find_merged_motif_recursion(motifs,motifs_occ,pos_list,org_start,org_end,win_size,num_occur,collision_table,counted_cls)
% This function finds the co-occur motif within a range start to and, with
% restriction of a certain number of occur. It use recursion to find all
% co-occur motifs
% The output are merged motifs with new starting point and ending point



% Find the collision with respect to the current motif,defined by start and
% end

id_ind = find(motifs(:,2)==org_start);
id = motifs(id_ind,1);
% Modification: as one motif may collide with the other motif multiple
% times if the other motif is not widely spaced, we count the ones wherer
% the number of collision is equal or greater than the number of occurance
% BUT!! when merging the motif, 3 values must be validated: the # of
% occurance of the original motif, the # of occurance the the target motif,
% and the value in the collision table 
% The number of occurance for the original motif and the occurance for the
% target motif must be the same while the number of the collision must be
% equal or greater than the # of occurance
cls_col = find(collision_table(id,:)>=num_occur);
cls_row = find(collision_table(:,id)>=num_occur);
cls_all = unique([cls_col';cls_row]);

% Remove ids that had been counted
cls_all = setdiff(cls_all,counted_cls);

% Initialization
merged_start = [org_start];
merged_id = [id];
merged_end = [org_end];

%fprintf('The index in now at %d\r\n',max(merged_end));

%fprintf('id %d cls with \r\n',id);
%cls_all

% Find the range the search will be conduct
pos_list_start = find(pos_list(:,1)==org_start);
pos_list_end   = find(pos_list(:,1)==org_end);


% Check if the motifs in between satisfy the condition
if length(cls_all)~=0 & ((pos_list_start+1) <= (pos_list_end-1))
    
    %for i = (pos_list_start+1):(pos_list_end-1)
        i = pos_list_start+1;
        % find if the motif with the specific position satisfy the
        % condition
        potential_motif = motifs(find(motifs(:,2)==pos_list(i,1)),:); 
        motif_endings = motifs(:,2)+win_size-1;
        
        
        if length(potential_motif) == 0 
            potential_motif = motifs(find(motif_endings(:,1)==pos_list(i,1)),:);
        end
        %fprintf('position at %d ending at %d\r\n',i,pos_list_end);
        if ismember(potential_motif(1,1),cls_all)
            % check if the overlapped are in the collsion group
            % if yes,contribute the starts and ends to the final group
                
            % One more thing: make sure the # of occruance is the same 
            occ_record = motifs_occ(find(motifs_occ(:,1)==potential_motif(1,1)),:);
            if occ_record(1,2)==num_occur
                
                counted_cls = [counted_cls;potential_motif(1,1)];
                cls_all = setdiff(cls_all,counted_cls);
                
                merged_start = [merged_start;potential_motif(1,2)];
                % Trial change
                %merged_end = [merged_end;potential_motif(1,2)+win_size-1];
                merged_end = [merged_end;potential_motif(1,3)];
                
                merged_id = [merged_id;potential_motif(1,1)];
                % Find the starts and ends recursively
                [new_start new_end new_id] = find_merged_motif_recursion(motifs,motifs_occ,pos_list,potential_motif(1,2),potential_motif(1,2)+win_size-1,win_size,num_occur,collision_table,counted_cls);
                merged_start = [merged_start;new_start];
                merged_end = [merged_end;new_end];
                merged_id = [merged_id;new_id];
                
            end
        end
        
    %end
   
else
    
    merged_start = [org_start];
    merged_end = [org_end];
    merged_id = [id];
    
end



