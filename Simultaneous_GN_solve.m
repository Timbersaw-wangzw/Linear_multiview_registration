function update_grpah= Simultaneous_GN_solve(init_graph)
%SIMULTANEOUS simultaneous solve r and t by linear
%   init_graph stors all the infos about node
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;

max_multi_icp=10;
anchor_idx=1;
init_graph.anchor_idx=1;
err_list=zeros(max_multi_icp+1,1);
err_list(1)=init_graph.rmse;

fprintf('start inital matrix \n');
graph=initalMatrixInfo(init_graph);
fprintf('end inital matrix \n');
fprintf('iteration: 0/%d, error: %.3f\n',max_multi_icp,init_graph.rmse);
for iter=1:max_multi_icp
    [A,b]=updateMatrix(graph);
    x=A\(-b);
    % update pose
    for v=2:node_num
        w=x((v-1)*6+1:6*v);
        T=SE3.exp(w);
        graph.node{v}.T=T;
    end
    % calculate error
    graph=getMultiRegError(graph);
    err_list(i+1)=graph.rmse;
    fprintf('iteration: %d/%d, error: %.3f\n',iter,max_multi_icp,graph.rmse);
end
end
function init_graph=initalMatrixInfo(init_graph)
anchor_idx=init_graph.anchor_idx;
edge_num=length(init_graph.edge);


for e=1:edge_num
    %     corre_num=floor(length(init_graph.edge{i}.pair_points)/100);
    corre_num=init_graph.edge{e}.orgin_pair_points_1.Count;
    fprintf('edge index:%d, correspondence points number: %d\n',[e,init_graph.edge{e}.orgin_pair_points_1.Count]);
    idx_i=init_graph.edge{e}.idx(1);
    idx_j=init_graph.edge{e}.idx(2);

    init_graph.edge{e}.matrix_p0=zeros(6,6);
    init_graph.edge{e}.matrix_q0=zeros(6,6);
    init_graph.edge{e}.sum_skew_p0=zeros(3,3);
    init_graph.edge{e}.sum_skew_q0=zeros(3,3);
    init_graph.edge{e}.sum_skew_qp0=zeros(3,3);
    init_graph.edge{e}.sum_p0=zeros(3,1);
    init_graph.edge{e}.sum_q0=zeros(3,1);
    init_graph.edge{e}.sum_qp0=zeros(3,3);

    xi=init_graph.edge{e}.orgin_pair_points_1.Location(:,1);
    yi=init_graph.edge{e}.orgin_pair_points_1.Location(:,2);
    zi=init_graph.edge{e}.orgin_pair_points_1.Location(:,3);

    xj=init_graph.edge{e}.orgin_pair_points_2.Location(:,1);
    yj=init_graph.edge{e}.orgin_pair_points_2.Location(:,2);
    zj=init_graph.edge{e}.orgin_pair_points_2.Location(:,3);


    matrix_p0_temp=cell(corre_num,1);
    matrix_q0_temp=cell(corre_num,1);
    skew_p0_temp=cell(corre_num,1);
    skew_q0_temp=cell(corre_num,1);
    qp0_temp=cell(corre_num,1);
    p0_temp=cell(corre_num,1);
    q0_temp=cell(corre_num,1);

    if idx_i~=anchor_idx&&idx_j~=anchor_idx
        type='no anchor';
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        if idx_j>anchor_idx
            j=idx_j-1;
        end
        init_graph.edge{e}.anchor_graph_idx=[i;j];
        parfor k=1:corre_num
            p0=[xi(k);yi(k);zi(k);1];
            q0=[xj(k);yj(k);zj(k);1];
            Mi0=dotVec(p0);
            Mj0=-1*dotVec(q0);
            matrix_p0_temp{k}=Mi0'*Mi0;
            matrix_q0_temp{k}=Mj0'*Mj0;
            skew_p0_temp{k}=skew(p0(1:3));
            skew_q0_temp{k}=skew(q0(1:3));
            qp0_temp{k}=q0(1:3)*p0(1:3)';
            p0_temp{k}=p0(1:3);
            q0_temp{k}=q0(1:3);
        end
    end
    if idx_i==anchor_idx
        type='anchor i';
        if idx_j>anchor_idx
            j=idx_j-1;
        end
        init_graph.edge{e}.anchor_graph_idx=[0;j];
        parfor k=1:corre_num
            p0=[xi(k);yi(k);zi(k);1];
            q0=[xj(k);yj(k);zj(k);1];
            Mj0=-1*dotVec(q0);
            matrix_q0_temp{k}=Mj0'*Mj0;
            p0_temp{k}=p0(1:3);
            q0_temp{k}=q0(1:3);
            qp0_temp{k}=q0(1:3)*p0(1:3)';
            skew_q0_temp{k}=skew(q0(1:3));
        end
    end
    if idx_j==anchor_idx
        type='anchor j';
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        init_graph.edge{e}.anchor_graph_idx=[i;0];
        parfor k=1:corre_num
            p0=[xi(k);yi(k);zi(k);1];
            q0=[xj(k);yj(k);zj(k);1];
            Mi0=dotVec(p0);
            matrix_p0_temp{k}=Mi0'*Mi0;
            p0_temp{k}=p0(1:3);
            qp0_temp{k}=q0(1:3)*p0(1:3)';
        end
    end

    for k=1:corre_num
        switch type
            case 'no anchor'
                init_graph.edge{e}.matrix_p0=init_graph.edge{e}.matrix_p0+matrix_p0_temp{k};
                init_graph.edge{e}.matrix_q0=init_graph.edge{e}.matrix_q0+matrix_q0_temp{k};
                init_graph.edge{e}.sum_skew_p0=init_graph.edge{e}.sum_skew_p0+skew_p0_temp{k};
                init_graph.edge{e}.sum_skew_q0=init_graph.edge{e}.sum_skew_q0+skew_q0_temp{k};
                init_graph.edge{e}.sum_qp0=init_graph.edge{e}.sum_qp0+qp0_temp{k};
                init_graph.edge{e}.sum_q0=init_graph.edge{e}.sum_q0+q0_temp{k};
                init_graph.edge{e}.sum_p0=init_graph.edge{e}.sum_p0+p0_temp{k};
            case 'anchor i'
                init_graph.edge{e}.matrix_q0=init_graph.edge{e}.matrix_q0+matrix_q0_temp{k};
                init_graph.edge{e}.sum_qp0=init_graph.edge{e}.sum_qp0+qp0_temp{k};
                init_graph.edge{e}.sum_skew_q0=init_graph.edge{e}.sum_skew_q0+skew_q0_temp{k};
                init_graph.edge{e}.sum_q0=init_graph.edge{e}.sum_q0+q0_temp{k};
                init_graph.edge{e}.sum_p0=init_graph.edge{e}.sum_p0+p0_temp{k};

            case 'anchor j'
                init_graph.edge{e}.matrix_p0=init_graph.edge{e}.matrix_p0+matrix_p0_temp{k};
                init_graph.edge{e}.sum_skew_p0=init_graph.edge{e}.sum_skew_p0+skew_p0_temp{k};
                init_graph.edge{e}.sum_qp0=init_graph.edge{e}.sum_qp0+qp0_temp{k};
                init_graph.edge{e}.sum_q0=init_graph.edge{e}.sum_q0+q0_temp{k};
                init_graph.edge{e}.sum_p0=init_graph.edge{e}.sum_p0+p0_temp{k};
        end
    end
    init_graph.edge{e}.edge_type=type;
end

end
function [A,b]=updateMatrix(graph)
edge_num=length(graph.edge);
node_num=length(graph.node)-1;
A=zeros(6*node_num,6*node_num);
b=zeros(6*node_num,1);
for e=1:edge_num
    i=graph.edge{e}.anchor_graph_idx(1);
    j=graph.edge{e}.anchor_graph_idx(2);
    type=graph.edge{e}.edge_type;
    idx_i=graph.edge{e}.idx(1);
    idx_j=graph.edge{e}.idx(2);
    Ti=graph.node{idx_i}.T;
    Tj=graph.node{idx_j}.T;
    corre_num=graph.edge{e}.orgin_pair_points_1.Count;
    switch type
        case 'no anchor'
            Mi=getDiagnoalBlock(Ti,graph.edge{e}.matrix_p0);
            Mj=getDiagnoalBlock(Tj,graph.edge{e}.matrix_q0);
            Mij=getNonDiagnoalBlock(graph.edge{e}.sum_skew_p0,graph.edge{e}.sum_skew_q0,...
                graph.edge{e}.sum_qp0,Ti,Tj,corre_num);
            bi=getNonAnchorb(graph.edge{e}.sum_p0,graph.edge{e}.sum_q0,...
                graph.edge{e}.sum_qp0,Ti,Tj,corre_num);
            A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)=A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)+Mi;
            A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)=A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)+Mj;
            A(6*(i-1)+1:6*i,6*(j-1)+1:6*j)=A(6*(i-1)+1:6*i,6*(j-1)+1:6*j)+Mij;
            A(6*(j-1)+1:6*j,6*(i-1)+1:6*i)=A(6*(j-1)+1:6*j,6*(i-1)+1:6*i)+Mij';
            b(6*(i-1)+1:6*i)=b(6*(i-1)+1:6*i)+bi;
            b(6*(j-1)+1:6*j)=b(6*(j-1)+1:6*j)-bi;

        case 'anchor i'
            Mj=getDiagnoalBlock(Tj,graph.edge{e}.matrix_q0);
            bi=getAnchorb(graph.edge{e}.sum_q0,graph.edge{e}.sum_p0,...
                graph.edge{e}.sum_qp0,graph.edge{e}.sum_skew_q0,Tj,corre_num);
            A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)=A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)+Mj;
            b(6*(j-1)+1:6*j)=b(6*(j-1)+1:6*j)+bi;
        case 'anchor j'
            Mi=getDiagnoalBlock(Tj,graph.edge{e}.matrix_p0);
            bi=getAnchorb(graph.edge{e}.sum_p0,graph.edge{e}.graph.edge{e}.sum_q0,...
                graph.edge{e}.sum_qp0,Ti,corre_num);
            A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)=A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)+Mi;
            b(6*(i-1)+1:6*i)=b(6*(i-1)+1:6*i)+bi;
    end
end
end
function M=getDiagnoalBlock(T,sum_m)
inv_Ad_T=invAd(T);
M=inv_Ad_T'*sum_m*inv_Ad_T;
end
function b=getAnchorb(sum_p,sum_q,sum_qp,sum_skew_p,Ti,N)
b=zeros(6,1);
Ri=double(SO3(Ti));
ti=double(transl(Ti)');
b(1:3)=Ri*sum_p+N*ti-sum_q;

b(4:6)=skew(ti)*Ri*sum_p+Ri*sum_skew_p*Ri'*ti-skew(ti)*sum_q;
i1=skew([1,0,0]);
i2=skew([0,1,0]);
i3=skew([0,0,1]);
a=-1*(sum_qp'*Ri');
bb(1,:)=trace(i1*a);
bb(2,:)=trace(i2*a);
bb(3,:)=trace(i3*a);
b(4:6,:)=b(4:6,:)-bb;
end
function b=getNonAnchorb(sum_p,sum_q,sum_qp,Ti,Tj,N)
b=zeros(6,1);
Ri=double(SO3(Ti));
ti=double(transl(Ti)');
Rj=double(SO3(Tj));
tj=double(transl(Tj)');
b(1:3)=Ri*sum_p+N*ti-Rj*sum_q-N*tj;
i1=skew([1,0,0]);
i2=skew([0,1,0]);
i3=skew([0,0,1]);
% Ri*p0(1:3)*q0(1:3)'*Rj'+Ri*p0(1:3)*tj'+ti*q0(1:3)'*Rj'+ti*tj'
a=-1*(Ri*sum_qp'*Rj'+Ri*sum_p*tj'+ti*sum_q'*Rj'+N*ti*tj');
b(4)=trace(i1*a);
b(5)=trace(i2*a);
b(6)=trace(i3*a);
end
function out=getNonDiagnoalBlock(sum_skew_p,sum_skew_q,sum_qp,Ti,Tj,N)
Ri=double(SO3(Ti));
Rj=double(SO3(Tj));
C=Ri'*Rj;
M=zeros(6,6);
M(1:3,1:3)=-1*C*N;
M(1:3,4:6)=C*sum_skew_q;
M(4:6,1:3)=-1*sum_skew_p*C;
M(4:6,4:6)=-1*(-1*trace(sum_qp)*eye(3)+sum_qp)*(-1*trace(C)*eye(3)+C')+trace(C'*sum_qp)*eye(3)-C'*sum_qp;
inv_Ad_Ti=invAd(Ti);
inv_Ad_Tj=invAd(Tj);
out=inv_Ad_Ti'*M*inv_Ad_Tj;
end
function M=invAd(T)
M=zeros(6,6);
rot=double(SO3(T));
t=double(transl(T));
M(1:3,1:3)=rot';
M(1:3,4:6)=-1*rot'*skew(t);
M(4:6,4:6)=rot';
end