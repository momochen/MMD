function Y = normalization(X)
% This function normalize the feature matrix to -1 to 1

% vecN = ((vec-minVec)./(maxVec-minVec) - 0.5 ) *2;

x_size = size(X);
Y = ones(x_size(1,1),x_size(1,2));

for i = 1:x_size(1,2)
    
   %Y(:,i) = X(:,i)/norm(X(:,i));
   Y(:,i) = ((X(:,i)-min(X(:,i)))./(max(X(:,i))-min(X(:,i))) - 0.5 ) * 2;
   
end
