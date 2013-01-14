

%find neighbours of vk
function vk_neighbours = find_neighbour(vertex_vk,n,motifs)
    vk_neighbours=[];
    for vj=2:n
        vertex_vj = motifs(vj,:);
        weight = weight_between_two_motifs(vertex_vk,vertex_vj,n);
        %vertex_vj
        %weight
        vertex_vj(1,n+3)=weight;
        %vertex_vj
        if weight > 0
            temp_vk = vk_neighbours;
            vk_neighbours = cat(1,temp_vk,vertex_vj);
        end
    end
    %vk_neighbours
end