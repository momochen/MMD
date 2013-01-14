function C = array_diff(A,B)
% Assuming both A and B are col vectors. 
% Find the difference between two arrays C = (A-B) and return the difference
% in sorted order. Notice this is not set diff
C = [];

for i = 1:length(B)
    
       index = find(A(:,1)==B(i,1));
       if ~isempty(index)
           A = [A(1:(index-1),1);A((index+1):end,1)];
       end
end

C = A;