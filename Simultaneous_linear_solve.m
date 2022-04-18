function update_grpah= Simultaneous_linear_solve(init_graph)
%SIMULTANEOUS simultaneous solve r and t by linear
%   init_graph stors all the infos about node
edge_num=length(init_graph.edge);
node_num=length(init_graph.node);
% A=zeros(4,6*edge_num);
b=zeros(6*node_num,1);
A=zeros(6*node_num,6*node_num);
for i=1:edge_num
    corre_num=length(init_graph.edge{i}.pair_points);
    pair_points=init_graph.edge{i}.pair_points;
    idx_a=init_graph.edge{i}.idx(1);
    idx_b=init_graph.edge{i}.idx(2);
    ax=pair_points(1,:);
    ay=pair_points(2,:);
    az=pair_points(3,:);
    bx=pair_points(4,:);
    by=pair_points(5,:);
    bz=pair_points(6,:);
    temp_points=zeros(6*node_num,1);
    A_temp={};
    for j=1:corre_num/100
        Aj=zeros(4,6*node_num);
        va=[ax(j);ay(j);az(j);1];
        vb=[bx(j);by(j);bz(j);1];
        Ma=dotVec(va);
        Mb=-1*dotVec(vb);
        Aj(:,(idx_a-1)*6+1:idx_a*6)=Ma;
        Aj(:,(idx_b-1)*6+1:idx_b*6)=Mb;
        A_temp{j}=Aj'*Aj;
        temp_points(:,j)=Aj'*(va-vb);
    end
    for j=1:corre_num
        A=A+A_temp{j};
        b=b-temp_points(:,j);
    end
    %% test LDU factorization
end
x=A\b;
end

