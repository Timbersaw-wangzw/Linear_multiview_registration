function graph = cloest_solution_multiICP_anchor(init_graph)
%SIMULTANEOUS_MANIFOLD 此处显示有关此函数的摘要
%   此处显示详细说明
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;
ee=eye(node_num);
A=zeros(3*node_num,3*node_num);
B=zeros(3*node_num,node_num);
C=zeros(node_num,node_num);
G=C;
D=B;
E=zeros(3*node_num,3);
F=zeros(3,node_num);
anchor_idx=1;
for i=1:node_num
    init_graph.node{i}.T=SE3.exp(rand(6,1));
end
for e=1:edge_num
    idx_i=init_graph.edge{e}.idx(1);
    idx_j=init_graph.edge{e}.idx(2);
    idx_noAnchor=init_graph.edge{e}.idx(2);
    corre_num=init_graph.edge{e}.orgin_pair_points_1.Count;
    A_temp=cell(corre_num,1);
    B_temp=cell(corre_num,1);
    C_temp=cell(corre_num,1);
    D_temp=cell(corre_num,1);
    E_temp=cell(corre_num,1);
    F_temp=cell(corre_num,1);
    G_temp=cell(corre_num,1);
    if idx_i~=anchor_idx&&idx_j~=anchor_idx
        type='no anchor';
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        if idx_j>anchor_idx
            idx_noAnchor=idx_j-1;
        end
        xi=init_graph.edge{e}.orgin_pair_points_1.Location(:,1);
        yi=init_graph.edge{e}.orgin_pair_points_1.Location(:,2);
        zi=init_graph.edge{e}.orgin_pair_points_1.Location(:,3);

        xj=init_graph.edge{e}.orgin_pair_points_2.Location(:,1);
        yj=init_graph.edge{e}.orgin_pair_points_2.Location(:,2);
        zj=init_graph.edge{e}.orgin_pair_points_2.Location(:,3);
        ei=ee(:,i);
        ej=ee(:,idx_noAnchor);
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
    if idx_i==anchor_idx||idx_j==anchor_idx
        type='anchor';
        if idx_j>anchor_idx
            idx_noAnchor=idx_j-1;
            xi=init_graph.edge{e}.orgin_pair_points_2.Location(:,1);
            yi=init_graph.edge{e}.orgin_pair_points_2.Location(:,2);
            zi=init_graph.edge{e}.orgin_pair_points_2.Location(:,3);

            xj=init_graph.edge{e}.orgin_pair_points_1.Location(:,1);
            yj=init_graph.edge{e}.orgin_pair_points_1.Location(:,2);
            zj=init_graph.edge{e}.orgin_pair_points_1.Location(:,3);
        end
        if idx_i>anchor_idx
            idx_noAnchor=idx_i-1;
            xi=init_graph.edge{e}.orgin_pair_points_1.Location(:,1);
            yi=init_graph.edge{e}.orgin_pair_points_1.Location(:,2);
            zi=init_graph.edge{e}.orgin_pair_points_1.Location(:,3);

            xj=init_graph.edge{e}.orgin_pair_points_2.Location(:,1);
            yj=init_graph.edge{e}.orgin_pair_points_2.Location(:,2);
            zj=init_graph.edge{e}.orgin_pair_points_2.Location(:,3);
        end
        ei=ee(:,idx_noAnchor);
        ei_kron=kron(ei,eye(3));
        parfor k=1:corre_num
            p0=[xi(k);yi(k);zi(k)];
            q0=[xj(k);yj(k);zj(k)];
            D_temp{k}=ei_kron*q0*ei';
            E_temp{k}=ei_kron*p0*q0';
            F_temp{k}=q0*ei';
            G_temp{k}=ei*ei';
        end
    end
    for k=1:corre_num
        switch type
            case 'no anchor'
                A=A+A_temp{k};
                B=B+B_temp{k};
                C=C+C_temp{k};
            case 'anchor'
                D=D+D_temp{k};
                E=E+E_temp{k};
                F=F+F_temp{k};
                G=G+G_temp{k};
        end
    end
end
mathcal_C=C+G;
tilde_B=B+D;
inv_mathcal_C=mathcal_C^-1;
AA=kron(A-tilde_B*inv_mathcal_C*tilde_B',eye(3));
temp=tilde_B*inv_mathcal_C*F'-E;
bb=-1*reshape(temp',[9*node_num,1]);
x=AA\bb;
mathcal_G=reshape(x,[3,3*node_num]);
mathcal_R=zeros(3,3*node_num);
for i=1:node_num
    G=mathcal_G(1:3,3*(i-1)+1:3*i);
    [U,D,V]=svd(G);
    s=eye(3);
    s(3,3)=det(U*V);
    mathcal_R(1:3,3*(i-1)+1:3*i)=U*s*V;
end
end



