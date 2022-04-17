clc
clear
inputdata='.\highspeedtrain';
addpath(inputdata);
subdir  = dir( fullfile(inputdata) );
nodeNum=1;
for i = 1 : length( subdir )
    if ~isempty(strfind(subdir(i).name,'.mat'))
        correspond_file=subdir(i).name;
        graph.node(nodeNum).filePath=correspond_file;
        graph.edge{nodeNum}.idx=str2double(regexp(correspond,'\d','match'));
        load(correspond_file);
        nodeNum=nodeNum+1;
    end
end


% 提前获取边信息
% graph.edge{1}.idx=[1,2];
% graph.edge{1}.correspondences=[1,2];
% graph.edge{2}.idx=[1,8];
% graph.edge{2}.correspondences=[1,2];
% graph.edge{3}.idx=[2,3];
% graph.edge{3}.correspondences=[1,2];
% graph.edge{4}.idx=[2,7];
% graph.edge{4}.correspondences=[2,7];
% graph.edge{5}.idx=[2,8];
% graph.edge{5}.correspondences=[1,2];
% graph.edge{6}.idx=[3,4];
% graph.edge{6}.correspondences=[1,2];
% graph.edge{7}.idx=[3,5];
% graph.edge{7}.correspondences=[3,5];
% graph.edge{8}.idx=[3,6];
% graph.edge{8}.correspondences=[3,5];
% graph.edge{9}.idx=[3,7];
% graph.edge{9}.correspondences=[3,5];
% graph.edge{10}.idx=[4,5];
% graph.edge{10}.correspondences=[3,5];
% graph.edge{11}.idx=[4,6];
% graph.edge{11}.correspondences=[3,5];
% graph.edge{12}.idx=[5,6];
% graph.edge{12}.correspondences=[3,5];
% graph.edge{13}.idx=[6,7];
% graph.edge{13}.correspondences=[3,5];
% edge_num=length(graph.edge);
% for i=1:edge_num
%     idx_p=graph.edge{i}.idx(1);
%     idx_q=graph.edge{i}.idx(2);
%     AreaPt1=importdata(graph.node(idx_p).filePath);
%     AreaPt2=importdata(graph.node(idx_q).filePath);
%     AreaPt1=AreaPt1(:,1:3);
%     AreaPt2=AreaPt2(:,1:3);
%
%     AreaPt1=[AreaPt1,ones(length(AreaPt1(:,1)), 1)]';
%     AreaPt2=[AreaPt2,ones(length(AreaPt2(:,1)), 1)]';
%     [pt1,pt2]=FilterPoints(AreaPt1,AreaPt2);
%     [~,~,~,pair_points]=IRLS_ICP(pt1(1:3,:),pt2(1:3,:),0.0001,'Cauchy',50);
% end
