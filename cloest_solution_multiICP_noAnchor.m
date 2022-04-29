function graph = cloest_solution_multiICP_noAnchor(init_graph)
%CLOEST_SOLUTION_MULTIICP 此处显示有关此函数的摘要
%   此处显示详细说明
edge_num=length(init_graph.edge);
node_num=length(init_graph.node);
ee=eye(node_num);
A=zeros(3*node_num,3*node_num);
B=zeros(3*node_num,node_num);
C=zeros(node_num,node_num);
for i=1:node_num
    init_graph.node{i}.T=SE3.exp(rand(6,1));
end
for e=1:edge_num
    i=init_graph.edge{e}.idx(1);
    j=init_graph.edge{e}.idx(2);

    xi=init_graph.edge{e}.orgin_pair_points_1.Location(:,1);
    yi=init_graph.edge{e}.orgin_pair_points_1.Location(:,2);
    zi=init_graph.edge{e}.orgin_pair_points_1.Location(:,3);

    xj=init_graph.edge{e}.orgin_pair_points_2.Location(:,1);
    yj=init_graph.edge{e}.orgin_pair_points_2.Location(:,2);
    zj=init_graph.edge{e}.orgin_pair_points_2.Location(:,3);

    corre_num=init_graph.edge{e}.orgin_pair_points_1.Count;
    A_temp=cell(corre_num,1);
    B_temp=cell(corre_num,1);
    C_temp=cell(corre_num,1);
    ei=ee(:,i);
    ej=ee(:,j);
    ei_kron=kron(ei,eye(3));
    ej_kron=kron(ej,eye(3));
    eij=ei-ej;
    parfor k=1:corre_num
        p0=[xi(k);yi(k);zi(k)];
        q0=[xj(k);yj(k);zj(k)];
        aij=ei_kron*p0-ej_kron*q0;
        A_temp{k}=aij*aij';
        B_temp{k}=aij*eij';
        C_temp{k}=eij*eij';
    end
    for k=1:corre_num
        A=A+A_temp{k};
        B=B+B_temp{k};
        C=C+C_temp{k};
    end
end
% A(1,:)=zeros(3*node_num,1);
% A(:,1)=zeros(1,3*node_num);
% B(1,:)=zeros(node_num,1);
% B(:,1)=zeros(1,3*node_num);
% C(1,:)=zeros(node_num,1);
% C(:,1)=zeros(1,node_num);
M=A-B*pinv(C)*B';
[U,D,V]=svd(M);
O=zeros(3,3*(node_num-1));
M*U(:,24)
R_group=[eye(3),O];
end

