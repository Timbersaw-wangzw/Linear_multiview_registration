clc
clear
addpath('../github_repo')
addpath('../IRLS')
loadData;
lambda=5*diag(ones(48,1));
Simultaneous_linear_solve(graph);

for i=1:8
    r(i)=rank(A((i-1)*6+1:6*i,:));
end
