function weight = weight_between_two_motifs(x,y,win_size)
coincident = 0;
for i=1:win_size
    if x(1,i)== y(1,i)
        coincident = coincident + 1;
    end
end
weight = coincident/x(1,win_size+1);
end