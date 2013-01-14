
%transform timeseries to symbol using SAX
function [motifs motifs_occ motif_id_pos_end] = momo_find_motif(str,win_size)

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert SAX string to matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
hash_limit = floor(win_size/2);
collision_limit = floor(hash_limit/2);

%SAX words construction

n = length(str);

row_num = n-win_size+1;
subsqs = zeros(row_num,win_size+1);
%fill in the matrix with SAX symbols

% momo{
  for i = 1:row_num
     subsqs(i,:) = [i str(1,i:(i+win_size-1))]; 
  end
% }momo

%subsqs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collision table via random projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
collision_table = zeros (row_num, row_num);

hash_process_time = 0;     
%random projection for hash_process_times
while hash_process_time ~= hash_limit
    i_j= randi(win_size,1,2);
    i=i_j(1,1);
    j=i_j(1,2);
    %got through sub_sequence table
    for r = 1:row_num-1
        for k = r+1:row_num
            if (subsqs(r,i+1)==subsqs(k,i+1) && subsqs(r,j+1)==subsqs(k,j+1))
                collision_table(k,r) = collision_table(k,r)+1;
            end
        end
    end
    hash_process_time = hash_process_time +1;
end
%collision_table

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find motif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% momo:
% If this part of the function is to find the motif, why not just find
% those entries with threshold X and then pick the corresponding row number
% , find the row in the subsqs matrix and output the motif? 

% momo
[row_index col_index] = find(collision_table>collision_limit);

% according to collision table, find corresponding row and column index,
% which represents the position of the subsequence

% Create matrix with id
id_pos_subsqs = [];
id = 1;

for i = 1:length(row_index)
    
    if i == 1
        % For the first time, add the subsequence on the collision table to
        % the new [id position subsqs] matrix
        id_pos_subsqs = [id_pos_subsqs;[id row_index(i,1) str(1,row_index(i,1):(row_index(i,1)+win_size-1))]];
        % If they are in the same position, i.e. same subsqs, then just add
        % one
        if row_index(i,1) ~= col_index(i,1)
            id_pos_subsqs = [id_pos_subsqs;[id col_index(i,1) str(1,col_index(i,1):(col_index(i,1)+win_size-1))]];
        end
        
    else
        
        % For the rest of time, check if there is subsequence that already
        % has an id. If there is, then add the same id; otherwise, add an
        % new id
        size_id_pos_subsqs = size(id_pos_subsqs);
        row_exists = 0;
        col_exists = 0;
        
        for j = 1:size_id_pos_subsqs(1,1)
           
           % Check duplicate motif w.r.t row index
           if row_index(i,1) == id_pos_subsqs(j,2) && row_exists == 0 
              %id_pos_subsqs = [id_pos_subsqs;[id_pos_subsqs(j,1) row_index(i,1) str(1,row_index(i,1):(row_index(i,1)+win_size-1))]];
              row_exists = 1;
           end
           % Check duplicate motif w.r.t col index
           if col_index(i,1) == id_pos_subsqs(j,2) && col_exists == 0
              %id_pos_subsqs = [id_pos_subsqs;[id_pos_subsqs(j,1) col_index(i,1) str(1,col_index(i,1):(col_index(i,1)+win_size-1))]];
              col_exists = 1;
           end
           
           if row_exists == 1 && col_exists == 1
              break; 
           end
           
        end
        
        % If no duplicate, then assign new id to them 
        dup_id = 0;
        if row_exists*col_exists == 0
            % either one of them is 0 which means there is duplicate motif,
            % and we need to find the id first
            if row_exists == 1
                % Find the id according to the row pos
                dup_id_pos = find(id_pos_subsqs(:,2)==row_index(i,1));
                dup_id = id_pos_subsqs(dup_id_pos(1,1),1);
                
            elseif col_exists == 1
                % Find the id according to the col pos
                dup_id_pos = find(id_pos_subsqs(:,2)==col_index(i,1));
                dup_id = id_pos_subsqs(dup_id_pos(1,1),1);
            else
                % This is absolutely new
                id = id + 1;
                id_pos_subsqs = [id_pos_subsqs;[id row_index(i,1) str(1,row_index(i,1):(row_index(i,1)+win_size-1))]];
                if row_index(i,1)~=col_index(i,1)
                    id_pos_subsqs = [id_pos_subsqs;[id col_index(i,1) str(1,col_index(i,1):(col_index(i,1)+win_size-1))]];
                end
            end
            
            % Add the one whose X_exists==0
            if dup_id~=0
                if row_exists == 0 
                    id_pos_subsqs = [id_pos_subsqs;[dup_id row_index(i,1) str(1,row_index(i,1):(row_index(i,1)+win_size-1))]];
                elseif col_exists == 0
                    id_pos_subsqs = [id_pos_subsqs;[dup_id col_index(i,1) str(1,col_index(i,1):(col_index(i,1)+win_size-1))]];
                end
            end
        
        end
        
    end
end


motifs = id_pos_subsqs;

% Calculate something I really need, [id pos end]
motif_id_pos_end = motifs(:,1:2);
motif_id_pos_end = [motif_id_pos_end motif_id_pos_end(:,2)+(win_size-1)*ones(size(motifs(:,2)))];

motifs_occ = [];


% momo

%{
%a new matrix that store motifs
motifs = [];
%%pick motif and store them into the vector
for i=1:row_num
    for j=i:row_num
        if collision_table(j,i)> 10
            %go through list 
            %if there already exist the motif do not put it in
            [size_of_motif size_col]=size(motifs);
            existi = 0;
            
            for k=1:size_of_motif
                if subsqs(i,:) == motifs(k,:)
                    existi = 1;
                end
            end
            
            if existi == 0
                % lu{
                temp_motif = subsqs(i,:);
                temp_motif2 = motifs;
                motifs = cat(1,temp_motif2,temp_motif);
                % }lu
                
                % momo{
                motifs = [motifs;subsqs(i,:)];
                % }momo
                
            end
            
        end
    end
end
%motifs

%}


%
%Create a new matrix store the motif,last column is number of occurences
%delete the duplicate motif
%calculte the number of occurences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % First col is the id and the second the # of occurance
  unique_motifs = [];
  new_motif_mat = [];
  counter = 1;
  motifs_mat_size = size(motifs);
  
  for counter = 1:motifs_mat_size(1,1)
     
      % Get the first line, while first col is # of occur and the rest is
      % the actual symbols
      one_motif = motifs(counter,:);
      % Find if there is any duplicate, find each one of them
      temp_occur = 0;
      unique_motifs_size = size(unique_motifs);
      motif_exists = 0;
      
      for i = 1:unique_motifs_size(1,1)
         if one_motif(1,1) == unique_motifs(i,1)
             
             unique_motifs(i,end) = unique_motifs(i,end) + 1;
             motif_exists = 1;
             break;
             
         end
      end
      
      if motif_exists == 0
          % new motif, then add 1
          unique_motifs = [unique_motifs;[one_motif(1,1) 1]];
      end
      
      % Continue
      
  end

  motifs_occ = unique_motifs;

