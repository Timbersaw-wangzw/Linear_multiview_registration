clc
clear
addpath('../github_repo')
addpath('../IRLS')

% 验证快速计算公式
a=[1;2;3];
b=[4;5;6];
skew_a=skew(a);
skew_b=skew(b);
W=magic(3);
fprintf('expand form\n')
-1*(-1*trace(b*a')*eye(3)+b*a')*(-1*trace(W)*eye(3)+W')+trace(W'*b*a')*eye(3)-W'*b*a'


fprintf('compad form\n')
skew_a*W*skew_b
% M=zeros(6,48);

% 验证秩
% load("A.mat")
% for i=1:8
%     M=M+A((i-1)*6+1:i*6,:);
% end

loadData;
lambda=5*diag(ones(48,1));
Simultaneous_linear_solve(graph);
