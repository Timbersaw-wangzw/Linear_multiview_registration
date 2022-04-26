function update_grpah= Simultaneous_GN_solve(init_graph)
%SIMULTANEOUS simultaneous solve r and t by linear
%   init_graph stors all the infos about node
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;
% A=zeros(4,6*edge_num);
j=zeros(6*node_num,1);
anchor_idx=1;
init_graph.anchor_idx=1;
% 测试快速运算是否正确，设定每个位姿有一个变换
% for v=1:node_num+1
%     if v==anchor_idx
%         vec=zeros(6,1);
%     else
%         vec=rand(6,1);
%     end
%     T=SE3.exp(vec);
%     init_graph.node{v}.T=T;
% end
t1=cputime;
[graph,~,~,~,~]=initalMatrixInfo(init_graph,false);
% ans=A_compare-A;
% lambda=5*diag(ones(6*node_num,1));
for i=1:max_iter
    [A0,b0]=updateMatrix(graph);
end
t2=cputime;
fprintf('parallel time:%.3f \n',t2-t1);
for iter=1:max_multi_icp


end
for e=1:node_num
    T=SE3.exp(x((e-1)*6+1:6*e));
    rot=double(SO3(T));
    trans=double(transl(T));
    tfrom=rigid3d(rot,trans);
end
end
function [init_graph,A,b,A_compare,b_compare]=initalMatrixInfo(init_graph,test)
anchor_idx=init_graph.anchor_idx;
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;
A=zeros(6*node_num,6*node_num);
b=zeros(6*node_num,1);
A_compare=zeros(6*node_num,6*node_num);
b_compare=zeros(6*node_num,1);
for e=1:edge_num
    %     corre_num=floor(length(init_graph.edge{i}.pair_points)/100);
    corre_num=init_graph.edge{e}.orgin_pair_points_1.Count;
    fprintf('edge index:%d, correspondence points number: %d\n',[e,init_graph.edge{e}.orgin_pair_points_1.Count]);
    idx_i=init_graph.edge{e}.idx(1);
    idx_j=init_graph.edge{e}.idx(2);
    Ti=init_graph.node{idx_i}.T;
    Tj=init_graph.node{idx_j}.T;

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

    b_temp_compare=cell(corre_num,1);
    A_temp_compare=cell(corre_num,1);
    if idx_i~=anchor_idx&&idx_j~=anchor_idx
        type='no anchor';
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        if idx_j>anchor_idx
            j=idx_j-1;
        end
        init_graph.edge{e}.anchor_graph_idx=[i;j];
        for k=1:corre_num
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


%             Ri=double(SO3(Ti));
%             ti=double(transl(Ti)');
%             Rj=double(SO3(Tj));
%             tj=double(transl(Tj)');
%             C=Ri'*Rj;
%             a0(1:3,:)=p0(1:3)+Ri'*ti-C*q0(1:3)-Ri'*tj;
%             i1=skew([1,0,0]);
%             i2=skew([0,1,0]);
%             i3=skew([0,0,1]);
%             aa0(1,:)=trace(i1*C*q0(1:3)*p0(1:3)');
%             aa0(2,:)=trace(i2*C*q0(1:3)*p0(1:3)');
%             aa0(3,:)=trace(i3*C*q0(1:3)*p0(1:3)');
%             aaa=skew(p0(1:3))*C*q0(1:3);
%             a0(4:6,:)=skew(p0(1:3))*Ri'*ti+aa0-skew(q0(1:3))*Ri'*tj;
%             x123=Ri*p0(1:3)*q0(1:3)'*Rj'+Ri*p0(1:3)*tj'+ti*q0(1:3)'*Rj'+ti*tj';
%             aa1(1,:)=trace(i1*x123);
%             aa1(2,:)=trace(i2*x123);
%             aa1(3,:)=trace(i3*x123);
%             b=getNonAnchorb(p0(1:3),q0(1:3),qp0_temp{k},Ti,Tj,1);

            %             M=invAd(Ti)';
            %             test=M*a0;
            %             test0(1:3,:)=p(1:3)-q(1:3);
            %             test0(4:6,:)=skew(q(1:3))*p(1:3);
            %             test1=Mi'*(p-q);
            %             test2(1:3,:)=-1*(p(1:3)-q(1:3));
            %             test2(4:6,:)=-1*skew(q(1:3))*p(1:3);
            %             test3=Mj'*(p-q);
            %             Ak_compare=zeros(4,6*node_num);
            %             Ak_compare(:,(i-1)*6+1:i*6)=Mi;
            %             Ak_compare(:,(j-1)*6+1:j*6)=Mj;
            %             aa=Ak_compare'*(p-q);
%                         Mij=getNonDiagnoalBlock(skew_p0_temp{k},skew_q0_temp{k},qp0_temp{k},Ti,Tj,1);
            if test
                Ak_compare=zeros(4,6*node_num);
                p=double(Ti)*p0;
                q=double(Tj)*q0;
                Mi=dotVec(p);
                Mj=-1*dotVec(q);
                Ak_compare(:,(i-1)*6+1:i*6)=Mi;
                Ak_compare(:,(j-1)*6+1:j*6)=Mj;
                qwe=Mi'*Mj;
                A_temp_compare{k}=Ak_compare'*Ak_compare;
                b_temp_compare{k}=Ak_compare'*(p-q);
            end
            %             
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
            %             Rj=double(SO3(Tj));
            %             tj=double(transl(Tj)');
            %
            %             b(4:6)=skew(tj)*Rj*q0(1:3)+Rj*skew_q0_temp{k}*Rj'*tj...
            %                 -skew(tj)*p0(1:3)-Rj*skew(q0(1:3))*Rj'*p0(1:3);
            if test
                p=[xi(k);yi(k);zi(k);1];
                q=double(Tj)*q0;
                Mj=-1*dotVec(q);
                Ak_compare=zeros(4,6*node_num);
                Ak_compare(:,(j-1)*6+1:j*6)=Mj;
                A_temp_compare{k}=Ak_compare'*Ak_compare;
                b_temp_compare{k}=Ak_compare'*(p-q);
            end
            %             bb=Ak_compare'*(p-q);
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
            if test
                Ak_compare=zeros(4,6*node_num);
                q=[xj(k);yj(k);zj(k);1];
                p=double(Ti)*p0;
                Mi=dotVec(p);
                Ak_compare(:,(i-1)*6+1:i*6)=Mi;
                A_temp_compare{k}=Ak_compare'*Ak_compare;
                b_temp_compare{k}=Ak_compare'*(p-q);
            end

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
        if test
            A_compare=A_compare+A_temp_compare{k};
            b_compare=b_compare+b_temp_compare{k};
        end
    end

    init_graph.edge{e}.edge_type=type;
    if test
        switch type
        case 'no anchor'
            Mi=getDiagnoalBlock(Ti,init_graph.edge{e}.matrix_p0);
            Mj=getDiagnoalBlock(Tj,init_graph.edge{e}.matrix_q0);
            Mij=getNonDiagnoalBlock(init_graph.edge{e}.sum_skew_p0,init_graph.edge{e}.sum_skew_q0,...
                init_graph.edge{e}.sum_qp0,Ti,Tj,corre_num);
            bi=getNonAnchorb(init_graph.edge{e}.sum_p0,init_graph.edge{e}.sum_q0,...
                init_graph.edge{e}.sum_qp0,Ti,Tj,corre_num);
            A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)=A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)+Mi;
            A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)=A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)+Mj;
            A(6*(i-1)+1:6*i,6*(j-1)+1:6*j)=A(6*(i-1)+1:6*i,6*(j-1)+1:6*j)+Mij;
            A(6*(j-1)+1:6*j,6*(i-1)+1:6*i)=A(6*(j-1)+1:6*j,6*(i-1)+1:6*i)+Mij';
            b(6*(i-1)+1:6*i)=b(6*(i-1)+1:6*i)+bi;
            b(6*(j-1)+1:6*j)=b(6*(j-1)+1:6*j)-bi;

        case 'anchor i'
            Mj=getDiagnoalBlock(Tj,init_graph.edge{e}.matrix_q0);
            bi=getAnchorb(init_graph.edge{e}.sum_q0,init_graph.edge{e}.sum_p0,...
                init_graph.edge{e}.sum_qp0,init_graph.edge{e}.sum_skew_q0,Tj,corre_num);
            A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)=A(6*(j-1)+1:6*j,6*(j-1)+1:6*j)+Mj;
            b(6*(j-1)+1:6*j)=b(6*(j-1)+1:6*j)+bi;
        case 'anchor j'
            Mi=getDiagnoalBlock(Tj,init_graph.edge{e}.matrix_p0);
            bi=getAnchorb(init_graph.edge{e}.sum_p0,init_graph.edge{e}.init_graph.edge{e}.sum_q0,...
                init_graph.edge{e}.sum_qp0,Ti,corre_num);
            A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)=A(6*(i-1)+1:6*i,6*(i-1)+1:6*i)+Mi;
            b(6*(i-1)+1:6*i)=b(6*(i-1)+1:6*i)+bi;
        end
    end
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