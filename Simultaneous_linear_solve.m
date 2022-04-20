function update_grpah= Simultaneous_linear_solve(init_graph)
%SIMULTANEOUS simultaneous solve r and t by linear
%   init_graph stors all the infos about node
edge_num=length(init_graph.edge);
node_num=length(init_graph.node)-1;
% A=zeros(4,6*edge_num);
b=zeros(6*node_num,1);
A=zeros(6*node_num,6*node_num);
anchor_idx=1;
for i=1:edge_num
%     corre_num=floor(length(init_graph.edge{i}.pair_points)/100);
    corre_num=init_graph.edge{i}.pair_points_1.Count;
    fprintf('edge index:%d, correspondence points number: %d\n',[i,init_graph.edge{i}.pair_points_1.Count]);
    idx_a=init_graph.edge{i}.idx(1);
    idx_b=init_graph.edge{i}.idx(2);
    
    ax=init_graph.edge{i}.pair_points_1.Location(:,1);
    ay=init_graph.edge{i}.pair_points_1.Location(:,2);
    az=init_graph.edge{i}.pair_points_1.Location(:,3);

    bx=init_graph.edge{i}.pair_points_2.Location(:,1);
    by=init_graph.edge{i}.pair_points_2.Location(:,2);
    bz=init_graph.edge{i}.pair_points_2.Location(:,3);

    temp_points=zeros(6*node_num,1);
    A_temp={};
    if idx_a~=anchor_idx&&idx_b~=anchor_idx
        if idx_a>anchor_idx
            a=idx_a-1;
        end
        if idx_b>anchor_idx
            b=idx_b-1;
        end
        parfor j=1:corre_num
            Aj=zeros(4,6*node_num);
            va=[ax(j);ay(j);az(j);1];
            vb=[bx(j);by(j);bz(j);1];
            Ma=dotVec(va);
            Mb=-1*dotVec(vb);
            Aj(:,(a-1)*6+1:a*6)=Ma;
            Aj(:,(b-1)*6+1:b*6)=Mb;
            A_temp{j}=Aj'*Aj;
            temp_points(:,j)=Aj'*(va-vb);
        end
    end
    if idx_a==anchor_idx
        if idx_b>anchor_idx
            b=idx_b-1;
        end
        parfor j=1:corre_num
            Aj=zeros(4,6*node_num);
            va=[ax(j);ay(j);az(j);1];
            vb=[bx(j);by(j);bz(j);1];
            Mb=-1*dotVec(vb);
            Aj(:,(b-1)*6+1:b*6)=Mb;
            A_temp{j}=Aj'*Aj;
            temp_points(:,j)=Aj'*(va-vb);
        end
    end
    if idx_b==anchor_idx
        if idx_a>anchor_idx
            a=idx_a-1;
        end
        parfor j=1:corre_num
            Aj=zeros(4,6*node_num);
            va=[ax(j);ay(j);az(j);1];
            vb=[bx(j);by(j);bz(j);1];
            Ma=dotVec(va);
            Aj(:,(a-1)*6+1:a*6)=Ma;
            A_temp{j}=Aj'*Aj;
            temp_points(:,j)=Aj'*(va-vb);
        end
    end
    for j=1:corre_num
        A=A+A_temp{j};
        b=b-temp_points(:,j);
    end
end
lambda=5*diag(ones(6*node_num,1));
x=(A+lambda)\b;
end
function [L,D,U]=LDU(pair_points)
corre_num=length(pair_points);
N=floor(corre_num/2000);
sum_x=[0;0;0];
sum_y=[0;0;0];
sum_ma=zeros(6,6);
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
    Ma=dotVec(va);
    Mb=-1*dotVec(vb);
    Aj(:,1:6)=Ma;
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