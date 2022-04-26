function graph = getMultiRegError(graph)
%GETMULTIREGERROR 此处显示有关此函数的摘要
%   此处显示详细说明
edge_num=length(graph.edge);
all_err=0;
all_num=0;
for e=1:edge_num
    idx_i=graph.edge{e}.idx(1);
    idx_j=graph.edge{e}.idx(2);
    Ti=graph.node{idx_i}.T;
    Tj=graph.node{idx_j}.T;
    Ri=double(SO3(Ti));
    ti=double(transl(Ti));
    Rj=double(SO3(Tj));
    tj=double(transl(Tj));
    tform1 = rigid3d(Ri,ti);
    tform2 = rigid3d(Rj,tj);
    pts1=graph.edge{e}.orgin_pair_points_1;
    pts2=graph.edge{e}.orgin_pair_points_2;
    ptsOut1 = pctransform(pts1,tform1);
    ptsOut2 = pctransform(pts2,tform2);
    err_temp=cell(ptsOut2.Count,1);
    r=ptsOut1.Location-ptsOut2.Location;
    parfor i=1:ptsOut2.Count
        err_temp{i}=r(i,:)*r(i,:)';
    end
    err=0;
    for i=1:ptsOut2.Count
        err=err+err_temp{i};
    end
    graph.edge{e}.rmse=sqrt(err/ptsOut2.Count);
    all_err=all_err+err;
    all_num=all_num+ptsOut2.Count;
end
graph.rmse=sqrt(all_err/all_num);
end

