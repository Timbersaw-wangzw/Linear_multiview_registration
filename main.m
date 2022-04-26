clc
clear
addpath('../github_repo')
addpath('../IRLS')

% 验证快速计算公式
% a=[1;2;3];
% b=[4;5;6];
% c=[5;6;7];
% d=[8;9;10];
% 
% skew_a=skew(a);
% skew_b=skew(b);
% skew_c=skew(c);
% skew_d=skew(d);
% W=magic(3);
% fprintf('expand form\n')
% sum_qp=b*a'+d*c';
% -1*(-1*trace(sum_qp)*eye(3)+sum_qp)*(-1*trace(W)*eye(3)+W')+trace(W'*sum_qp)*eye(3)-W'*sum_qp

% 5
% 
% fprintf('compad form\n')
% skew_a*W*skew_b+skew_c*W*skew_d

% M=zeros(6,48);

% 验证秩
% load("A.mat")
% for i=1:8
%     M=M+A((i-1)*6+1:i*6,:);
% end

loadData;
lambda=5*diag(ones(48,1));
Simultaneous_GN_solve(graph);
