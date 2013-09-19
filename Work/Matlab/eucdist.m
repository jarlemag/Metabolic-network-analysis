%eucdist.m
%Computes the euclidean distance between two vectors.
%Vectors must be of equal length!

%Input: Vector 1, Vector 2
%Output: Euclidean distance

%Is called by fmindistance.m

function ds = eucdist(vector1,vector2)

ds=0;

diff = vector1-vector2;
diffsquared = diff.^2;
ds = sqrt(sum(diffsquared)); 
end