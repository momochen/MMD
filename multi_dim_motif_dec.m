function clustered_motifs = multi_dim_motif_dec(motif_collection_with_occ,incidence_table_weight,win_size)   
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%!!!Algorithm 1 Multi dimensional motif construction!!!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

alfa=1;
 
%Sort the vertices vi in G based on the number of incidents (sizei)
[num_of_motifs col]=size(motif_collection_with_occ);
sorted_motifs_occ = sortrows(motif_collection_with_occ,-col);
sorted_motifs_occ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cluster the motifs of vk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each vertex vk in the sorted list of vertices do
%incidence_table_weight
for i = 1:num_of_motifs 
    clustered_motifs(i).vertex = i;
    clustered_motifs(i).sk=[];
    for j = 1:num_of_motifs
        if incidence_table_weight(i,j) >= alfa
            clustered_motifs(i).sk = cat(1,clustered_motifs(i).sk,j);
        end
    end
    vertex = clustered_motifs(i).vertex
    sk = clustered_motifs(i).sk
end
    

%clustered_motifs

end





%{
    vertex_vk = vertex(vk).motif_occ_no;
    neighbours_temp = find_neighbour(vertex(vk).motif_occ_no,win_size,sorted_motifs_occ);
    %get vk neighbours and no
    vertex(vk).neighbours = neighbours_temp;
    vertex(vk).neighbours(:,col_of_motifs+2)=[];

    [vk_neighbour_num,vk_neighbour_col]=size(vertex(vk).neighbours);
    
    %neighbournum_weight is a n*2 vector which left col is no. of motif and
    %right col is weight between this motif and its neighours
    vertex(vk).neighbournum_weight = neighbours_temp(:,col_of_motifs+1:col_of_motifs+2);
    vk_nei_num_weight = vertex(vk).neighbournum_weight;
    %vk_nei_num_weight
    
    %go through the neighbour of vk
    %just go through vk.neighbournum_weight vector
    for vj = 1:vk_neighbour_num
        
        %if ek,j > ? do
        if vertex(vk).neighbournum_weight(vj,2) > alfa
            %Add vj to Sk
            temp_s_clustered_motifs = vertex(vk).s_clustered_motifs;
            vertex(vk).s_clustered_motifs = cat(1,temp_s_clustered_motifs,vertex(vk).neighbournum_weight(vj,1));
        end
        %Update the weight of edges connected to vj by removing
        %the motif occurrences of rj that has coincident with rk
        % go through the entire motif list 
       
        % find vj and change its occurnece value

        %the number of coincident between rk and rj are the
        %minimun occurence value of these two motifs,
        %because they will have coinicidence every time it
        %occurs. So just minius the min
        %{
        for find_vj = 1:num_of_motifs
            min_vk_vj = min(vertex_vk(1,col_of_motifs),vertex_vj(1,col_of_motifs));
            minaa=min_vk_vj;
            minaa
            vertex(find_vj).motif_occ_no(1,col_of_motifs) = vertex(find_vj).motif_occ_no(1,col_of_motifs) - min_vk_vj;
            vertex(find_vj).motif_occ_no(1,col_of_motifs)
        end
        %}
    end
    %Update the sorted list of vertices
    sorted_motifs_occ = sortrows(motifs_occ_no,-col_of_motifs);

    %}
    