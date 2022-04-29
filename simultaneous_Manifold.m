function [outputArg1,outputArg2] = simultaneous_Manifold(init_graph)
%SIMULTANEOUS_MANIFOLD 此处显示有关此函数的摘要
%   此处显示详细说明
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;
ee=eye(node_num);
MABC=zeros(4*node_num,4*node_num);
A=zeros(3*node_num,3*node_num);
B=zeros(3*node_num,node_num);
C=zeros(node_num,node_num);
D=zeors(3*node_num,3);
E=zeors(3,node_num);
anchor_idx=1;
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
    D_temp=cell(corre_num,1);
    E_temp=cell(corre_num,1);
    if idx_i~=anchor_idx&&idx_j~=anchor_idx
        type='no anchor';
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        if idx_j>anchor_idx
            j=idx_j-1;
        end
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
    end
    if idx_i~=anchor_idx
        type='anchor i';
        if idx_j>anchor_idx
            j=idx_j-1;
        end
        ej=ee(:,j);
        ei=zeros(node_num,1);
        ej_kron=kron(ej,eye(3));
        eij=ej-ei;
        parfor k=1:corre_num
            p0=[xi(k);yi(k);zi(k)];
            q0=[xj(k);yj(k);zj(k)];
            aij=ej_kron*q0;
            B_temp{k}=aij*eij';
            C_temp{k}=eij*eij';
            D_temp{k}=ei_kron*q0*p0';
            E_temp{k}=p0*ej';
            
        end
    end
    if idx_j==anchor_idx
        type='anchor j';
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        ei=ee(:,i);
        ej=zeros(node_num,1);
        ei_kron=kron(ei,eye(3));
        eij=ei-ej;
        parfor k=1:corre_num
            p0=[xi(k);yi(k);zi(k)];
            q0=[xj(k);yj(k);zj(k)];
            aij=ei_kron*p0;
            B_temp{k}=aij*eij';
            C_temp{k}=eij*eij';
            D_temp{k}=ei_kron*p0*q0';
            E_temp{k}=q0*ei';          
        end
    end
    for k=1:corre_num
        A=A+A_temp{k};
        B=B+B_temp{k};
        C=C+C_temp{k};
        D=D+D_temp{k};
        E=E+E_temp{k};
    end
end
A1n=A[1:3,3:3*N-3];
Ann=A[3:3*N-3,3:3*N-3];
B1n=B[1:3,3:3*N-3];
end

