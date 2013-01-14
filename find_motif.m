
%transform timeseries to symbol using SAX
function [motifs,motifs_occ] = find_motif(str,n,win_size)

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert SAX string to matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SAX words construction

% Momo{
%n = length(str);
% }Momo

row_num = n-win_size+1;
subsqs = zeros(row_num,win_size+1);
%fill in the matrix with SAX symbols

% lu{
for i = 1:row_num
    for j = 1:win_size
        subsqs(i,j+1) = str (1,i+j-1);
    end
    subsqs(i,1) = i;
end
% }lu

%{
% momo{
  for i = 1:row_num
     subsqs(i,:) = [i str(1,i:(i+win_size-1))]; 
  end
% }momo
%}

%subsqs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Collision table via random projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
collision_table = zeros (row_num, row_num);

hash_process_time = 0;     
%random projection for hash_process_times
while hash_process_time ~= 40
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

%{
% momo:
% If this part of the function is to find the motif, why not just find
% those entries with threshold X and then pick the corresponding row number
% , find the row in the subsqs matrix and output the motif? 

% momo
[row_index col_index] = find(collision_table>10);
motifs = [];

for i = 1:length(row_index)
   
    motifs = [motifs;[collision_table(row_index(i,1),col_index(i,1)) subsqs(row_index(i,1),:)]];
    
end

% momo
%}


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




%
%Create a new matrix store the motif,last column is number of occurences
%delete the duplicate motif
%calculte the number of occurences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
% momo {
  unique_motifs = [];
  new_motif_mat = [];
  
  while length(motifs)>0
     
      % Get the first line, while first col is # of occur and the rest is
      % the actual symbols
      one_motif = motifs(1,:);
      % Find if there is any duplicate, find each one of them
      motifs_mat_size = size(motifs);
      temp_occur = 0;
      
      for i = 1:length(motifs_mat_size(1,1))
         if sum(one_motif(1,2:end) - motifs(i,2:end),2) == 0
             % Sum the occur # 
             temp_occur = temp_occur + motifs(i,1);
         else
             % Remove the duplicates ones, update motifs matrix
             new_motif_mat = [new_motif_mat;motifs(i,:)];
         end
      end
      unique_motifs = [unique_motifs;[temp_occur one_motif(1,2:end)]];
      motifs = new_motif_mat;
      
      % Continue
      
  end

   motifs = unique_motifs;
   motifs_occ = [];
% }momo
%}



%initialize the first row, make it equals to the first motif, the
%occurence of this motif is temporary 1
motifs_occ = motifs(1,:);
motifs_occ(1,win_size+2)=1;
% ???? why col first and row second???
[col_of_old_motif_table row_of_old_motif_table]=size(motifs);

for i=2:col_of_old_motif_table

    exist = 0;
    
    %a temp vector store the values of existing motifs for comparision
    motifs_occ_temp = motifs_occ;
    motifs_occ_temp(:,win_size+2)=[];
    [col_of_new_motif_table row_of_new_motif_table]=size(motifs_occ_temp);
    
    %go through the new motif table to find if there already exist a same motif
    %if exist, last col ++
    for j = 1:col_of_new_motif_table
        if motifs(i,2:win_size+1) == motifs_occ_temp(j,2:win_size+1)
            motifs_occ(j,win_size+2)= motifs_occ(j,win_size+2)+1;
            exist = 1;
        end
    end 
    %if not exist, add new motif
    if exist == 0
        add_new_motif_temp = motifs(i,:);
        add_new_motif_temp(:,win_size+2)=1;
        temp_motif3 = motifs_occ;
        motifs_occ = cat(1,temp_motif3,add_new_motif_temp);
    end
end
[num_of_motifs,col_of_motifs]=size(motifs_occ);
%motifs_occ

%plot motifs in graph




%{
%give each motif a No. which stored in last column
%%%%%% the matrix is consist of motif+occurence+no. now
temp_no = 1:num_of_motifs;
temp_no = temp_no';
motifs_occ_no = cat(2,motifs_occ,temp_no);
motifs_occ_no
%}
end



