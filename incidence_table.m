function [motif_collection all_motif_length_occ incidence_table_weight] = incidence_table(motifs,motif_length_occ,n)

   
 % total motif collection
 % a matrix that store all the motifs in all dimension
 % struct of the matrix [dim_num][time][motif 1 to n][length][occurence]
 all_motif_length_occ = [];
 [sb dim_num] = size(motif_length_occ);
 for i=1:dim_num
     [num_motif_i sb2] = size(motif_length_occ(i).motif);
     for j = 1:num_motif_i
         temp = cat(2,i,motif_length_occ(i).motif(j,:));
         all_motif_length_occ = cat(1,all_motif_length_occ,temp);
     end
 end
all_motif_length_occ

%give each motif a no. 
%which is the first begining col of the matrix
% struct of the matrix now: [motif_no.][dim_num][time][motif][occurence][length]
[total_motif_num collection_col]=size(all_motif_length_occ);
temp_no = 1:total_motif_num;
temp_no = temp_no';
all_motif_length_occ = cat(2,temp_no,all_motif_length_occ);
all_motif_length_occ


% total motif collection without occurence
 % a matrix that store all the motifs in all dimension
 % struct of the matrix [dim_num][time][motif][length]
 motif_collection = [];
 [sb2 dim_num2] = size(motifs);
 for i=1:dim_num2
     [num_motif_i sb3] = size(motifs(i).motif);
     for j = 1:num_motif_i
         temp = cat(2,i,motifs(i).motif(j,:));
         motif_collection = cat(1,motif_collection,temp);
     end
 end
motif_collection

%give no. to each motif in the orginial list
%which means motif occurs more than ones has same no. value
%the structure of matrix is now [no.][dim][time][motif][length]
motif_collection_with_no = [];
[num_of_total_motif_collection sb4]=size(motif_collection);
for i=1:num_of_total_motif_collection
    for j =1:total_motif_num
        if motif_collection(i,1) == all_motif_length_occ(j,2) 
            if motif_collection(i,3:(3+n)) == all_motif_length_occ(j,4:(4+n))
            temp = [all_motif_length_occ(j,1),motif_collection(i,:)];
            motif_collection_with_no=cat(1,motif_collection_with_no,temp);
            end
        end
    end
end
motif_collection_with_no

%a new matrix only store [no.][dim][time][length] of the motif
simple_collection = [motif_collection_with_no(:,1:3) motif_collection_with_no(:,n+4)];
simple_collection


%construct a matrix of 0s and 1s to present the motif
motif_0_1 = simple_collection(:,1:2);
motif_0_1(:,3:n+2)=0;
for i=1:total_motif_num
    motif_0_1(i,2+simple_collection(i,3):(simple_collection(i,3)+simple_collection(i,4)+1))=1;
end
motif_0_1


%a collision table to check the incidence between motifs in different time
%series
%in the incidence table, each row number means motif no.
incidence_table = zeros(total_motif_num, total_motif_num);

%fill in the collision table
%for each motif in one dimension check the incidence between the other
%motifs in other dimension

%AND each motif, if the answer = 0, it means they do not overlap
for i=1:total_motif_num
    for j=1:total_motif_num
    	if motif_0_1(j,2) ~= motif_0_1(i,2)
            
            a=motif_0_1(j,3:n+2);
            b=motif_0_1(i,3:n+2);
            temp=a&b;
            
            if sum(temp) ~= 0
                incidence_table(motif_0_1(i,1),motif_0_1(j,1))=incidence_table(motif_0_1(i,1),motif_0_1(j,1))+1;
            end
        end
    
    end
end
incidence_table;

%{
%go through each motif by searching the no. from 1 to n

for i=1:num_of_total_motif_collection
    % for motif i, go throuth the list and check every motif j in other dimension
    for j = 1:num_of_total_motif_collection
            
        %if the motif is in different dimension
        if motif_collection_with_no(j,2) ~= motif_collection_with_no(i,2)
            %check whether motif i have coincident with j
            
            %three cases
            %1.fisrt begining, just check from 0
            if motif_collection_with_no(i,3) < win_size
                if motif_collection_with_no(j,3)>0 && motif_collection_with_no(j,3) < motif_collection_with_no(i,3)+win_size
                    incidence_table(motif_collection_with_no(i,1),motif_collection_with_no(j,1))=incidence_table(motif_collection_with_no(i,1),motif_collection_with_no(j,1))+1;
                end
            end 
            %2. in the middle, check both previous and future
            if motif_collection_with_no(i,3) >= win_size && motif_collection_with_no(j,3) <= n-win_size-2
                if motif_collection_with_no(j,3)>motif_collection_with_no(i,3)-win_size && motif_collection_with_no(j,3) < motif_collection_with_no(i,3)+win_size
                    incidence_table(motif_collection_with_no(i,1),motif_collection_with_no(j,1))=incidence_table(motif_collection_with_no(i,1),motif_collection_with_no(j,1))+1;
               
                end
            end
            %3. end, check up tp total time
            if motif_collection_with_no(i,3) > n-win_size-2 && motif_collection_with_no(j,3) <= n-win_size+1
                if (motif_collection_with_no(j,3)>(motif_collection_with_no(i,3)-win_size)) && (motif_collection_with_no(j,3) <= n)
                      incidence_table(motif_collection_with_no(i,1),motif_collection_with_no(j,1))=incidence_table(motif_collection_with_no(i,1),motif_collection_with_no(j,1))+1;
               
                end
            end
            
        
        end
    end
end
%}



%incidence_table with weight
incidence_table_weight = zeros(total_motif_num, total_motif_num);
for i=1:total_motif_num
    for j=1:total_motif_num
        incidence_table_weight(i,j)=(incidence_table(i,j)/all_motif_length_occ(i,collection_col+1));        
    end
end
incidence_table_weight

end