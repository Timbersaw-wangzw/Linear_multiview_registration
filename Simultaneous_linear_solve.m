function update_grpah= Simultaneous_linear_solve(init_graph)
%SIMULTANEOUS simultaneous solve r and t by linear
%   init_graph stors all the infos about node
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;
% A=zeros(4,6*edge_num);
j=zeros(6*node_num,1);
A=zeros(6*node_num,6*node_num);
A_compare=zeros(6*node_num,6*node_num);

anchor_idx=1;
% 测试快速运算是否正确，设定每个位姿有一个变换
for v=1:node_num+1
    if v==anchor_idx
        vec=zeros(6,1);
    else
        vec=rand(6,1);
    end
    T=SE3.exp(vec);
    init_graph.node{v}.T=T;
end
for e=3:edge_num
%     corre_num=floor(length(init_graph.edge{i}.pair_points)/100);
    corre_num=init_graph.edge{e}.pair_points_1.Count;
    fprintf('edge index:%d, correspondence points number: %d\n',[e,init_graph.edge{e}.pair_points_1.Count]);
    idx_i=init_graph.edge{e}.idx(1);
    idx_j=init_graph.edge{e}.idx(2);
    Ti=init_graph.node{idx_i}.T;
    Tj=init_graph.node{idx_j}.T;


    xi=init_graph.edge{e}.pair_points_1.Location(:,1);
    yi=init_graph.edge{e}.pair_points_1.Location(:,2);
    zi=init_graph.edge{e}.pair_points_1.Location(:,3);

    xj=init_graph.edge{e}.pair_points_2.Location(:,1);
    yj=init_graph.edge{e}.pair_points_2.Location(:,2);
    zj=init_graph.edge{e}.pair_points_2.Location(:,3);

    temp_points=zeros(6*node_num,1);
    A_temp={};
    A_temp_compare={};
    a_temp={};
    skew_p_temp={};
    skew_q_temp={};
    qp_temp={};
    if idx_i~=anchor_idx&&idx_j~=anchor_idx
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        if idx_j>anchor_idx
            j=idx_j-1;
        end
        parfor k=1:corre_num
            Ak=zeros(4,6*node_num);
            p0=[xi(k);yi(k);zi(k);1];
            q0=[xj(k);yj(k);zj(k);1];
            Mi0=dotVec(p0);
            Mj0=-1*dotVec(q0);
            Ak(:,(i-1)*6+1:i*6)=Mi0;
            Ak(:,(j-1)*6+1:j*6)=Mj0;
            A_temp{k}=Ak'*Ak;
            temp_points(:,k)=Ak'*(p0-q0);
            skew_p_temp{k}=skew(p0(1:3));
            skew_q_temp{k}=skew(q0(1:3));
            qp_temp{k}=q0(1:3)*p0(1:3)';


            Ak_compare=zeros(4,6*node_num);
            p=double(Ti)*[xi(k);yi(k);zi(k);1];
            q=double(Tj)*[xi(k);yi(k);zi(k);1];
            Mi=dotVec(p);
            Mj=-1*dotVec(q);
            Ak_compare(:,(i-1)*6+1:i*6)=Mi;
            Ak_compare(:,(j-1)*6+1:j*6)=Mj;
            A_temp_compare{k}=Ak_compare'*Ak_compare;
        end
    end
    if idx_i==anchor_idx
        if idx_j>anchor_idx
            j=idx_j-1;
        end
        parfor k=1:corre_num
            Ak=zeros(4,6*node_num);
            p0=[xi(k);yi(k);zi(k);1];
            q0=[xj(k);yj(k);zj(k);1];
            Mj0=-1*dotVec(q0);
            Ak(:,(j-1)*6+1:j*6)=Mj0;
            a_temp{k}=Mj0'*Mj0;
            A_temp{k}=Ak'*Ak;
            temp_points(:,k)=Ak'*(p0-q0);

            Ak_compare=zeros(4,6*node_num);
            q=double(Tj)*[xi(k);yi(k);zi(k);1];
            Mj=-1*dotVec(q);
            Ak_compare(:,(j-1)*6+1:j*6)=Mj;
            A_temp_compare{k}=Ak_compare'*Ak_compare;
        end
    end
    if idx_j==anchor_idx
        if idx_i>anchor_idx
            i=idx_i-1;
        end
        parfor k=1:corre_num
            Ak=zeros(4,6*node_num);
            p0=[xi(k);yi(k);zi(k);1];
            q0=[xj(k);yj(k);zj(k);1];
            Mi0=dotVec(p0);
            Ak(:,(i-1)*6+1:i*6)=Mi0;
            A_temp{k}=Ak'*Ak;
            temp_points(:,k)=Ak'*(p0-q0);

            Ak_compare=zeros(4,6*node_num);
            p=double(Ti)*[xi(k);yi(k);zi(k);1];
            q=double(Tj)*[xi(k);yi(k);zi(k);1];
            Mi=dotVec(p);
            Ak_compare(:,(i-1)*6+1:i*6)=Mi;
            A_temp_compare{k}=Ak_compare'*Ak_compare;
        end
    end
    a=zeros(6,6);
    sum_skew_p=zeros(3,3);
    sum_skew_q=zeros(3,3);
    sum_qp=zeros(3,3);
    for k=1:corre_num
        A=A+A_temp{k};
        j=j-temp_points(:,k);
%         a=a+a_temp{k};
        sum_skew_p=sum_skew_p+skew_p_temp{k};
        sum_skew_q=sum_skew_q+skew_q_temp{k};
        sum_qp=sum_qp+qp_temp{k};
        A_compare=A_compare+A_temp_compare{k};
    end
    M=getDiagnoalBlock(Tj,A(1:6,1:6));
    M1=getNonDiagnoalBlock(sum_skew_p,sum_skew_q,sum_qp,Ti,Tj,corre_num);
end
% lambda=5*diag(ones(6*node_num,1));
x=A\j;
for e=1:node_num
    T=SE3.exp(x((e-1)*6+1:6*e));
    rot=double(SO3(T));
    trans=double(transl(T));
    tfrom=rigid3d(rot,trans);
end
end
function M=getDiagnoalBlock(T,sum_m)
inv_Ad_T=Ad(T)^-1;
M=inv_Ad_T'*sum_m*inv_Ad_T;
end
function M=getNonDiagnoalBlock(sum_skew_p,sum_skew_q,sum_qp,Ti,Tj,N)
rot1=double(SO3(Ti));
rot2=double(SO3(Tj));
W=rot1'*rot2;
M=zeros(6,6);
M(1:3,1:3)=-1*W*N;
M(1:3,4:6)=W*sum_skew_q;
M(4:6,1:3)=-1*W*sum_skew_p;
M(4:6,4:6)=-1*(-1*trace(sum_qp)*eye(3)+sum_qp)*(-1*trace(W)*eye(3)+W')+trace(W'*sum_qp)*eye(3)-W'*sum_qp;
inv_Ad_Ti=Ad(Ti)^-1;
inv_Ad_Tj=Ad(Tj)^-1;
M=inv_Ad_Ti'*M*inv_Ad_Tj;
end


function [L,D,U]=LDU(pair_points)
corre_num=length(pair_points);
N=floor(corre_num/2000);
sum_x=[0;0;0];
sum_y=[0;0;0];
sum_Mi=zeros(6,6);
sum_mb=zeros(6,6);
ax=pair_points(1,:);
ay=pair_points(2,:);
az=pair_points(3,:);
bx=pair_points(4,:);
by=pair_points(5,:);
bz=pair_points(6,:);
Aj=zeros(4,12);
A=zeros(12,12);
for j=1:N
    va=[ax(j);ay(j);az(j);1];
    vb=[bx(j);by(j);bz(j);1];
    Mi=dotVec(va);
    Mb=-1*dotVec(vb);
    Aj(:,1:6)=Mi;
    Aj(:,7:12)=Mb;
    A=A+Aj'*Aj;
    sum_x=sum_x+va(1:3);
    sum_y=sum_y+vb(1:3);
end
La=eye(6);
Lb=eye(6);
ave_x=sum_x/N;
ave_y=sum_y/N;
La(4:6,1:3)=skew(ave_x);
Lb(4:6,1:3)=skew(ave_y);
Daa=eye(6);
Dbb=eye(6);
Dab=eye(6);
Dba=eye(6);
Daa(1:3,1:3)=N*eye(3);
Dbb(1:3,1:3)=N*eye(3);
Dab(1:3,1:3)=-1*N*eye(3);
Dba(1:3,1:3)=-1*N*eye(3);
Ua=eye(6);
Ua(1:3,4:6)=-1*skew(ave_x);
Ub=eye(6);
Ub(1:3,4:6)=-1*skew(ave_y);
for i=1:N
    xi=pair_points(1:3,i);
    yi=pair_points(4:6,i);
    mi=skew(xi-ave_x);
    ni=skew(yi-ave_y);
    Daa(4:6,4:6)=Daa(4:6,4:6)-mi*mi;
    Dbb(4:6,4:6)=Dbb(4:6,4:6)-ni*ni;
    Dab(4:6,4:6)=Dab(4:6,4:6)+mi*ni;
    Dba(4:6,4:6)=Dba(4:6,4:6)+ni*mi;
end
ta=La*Daa*Ua;
tb=Lb*Dbb*Ub;
O=zeros(6,6);
L=[La,O;O,Lb];
D=[Daa,Dab;Dba,Dbb];
U=[Ua,O;O,Ub];
t=L*D*U;
end