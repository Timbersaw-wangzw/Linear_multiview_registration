inputdata='.\highspeedtrain';
addpath(inputdata);
subdir  = dir( fullfile(inputdata) );
nodeNum=1;
edgeNum=1;
for i = 1 : length( subdir )
    if ~isempty(strfind(subdir(i).name,'.txt'))
        if isempty(strfind(subdir(i).name,'filter'))
            graph.node{nodeNum}.filePath=subdir(i).name;
            graph.node{nodeNum}.idx=nodeNum;
            nodeNum=nodeNum+1;
        end
    end
    if ~isempty(strfind(subdir(i).name,'.mat'))
        correspond_file=subdir(i).name;
        graph.edge{edgeNum}.idx=str2double(regexp(correspond_file,'\d','match'));
        load(correspond_file);
%         graph.edge{edgeNum}.pair_points=pair_points;
        temp_source_points=pointCloud(pair_points(1:3,:)');
        temp_target_points=pointCloud(pair_points(4:6,:)');
        %         figure(1);
        %         pcshowpair(source_points,target_points);
%         gridStep = 2;
%         [source_points,idx] = pcdownsample(temp_source_points,'gridAverage',gridStep);
        
%                 ptCloud2 = pcdownsample(target_points,'gridAverage',gridStep);
        
%         [ptCloud1,test] = pcdownsample(source_points,'random',0.5);
%         ptCloud2 = pcdownsample(target_points,'random',0.5);
        %         figure(2);
%                 pcshowpair(ptCloud1,ptCloud2);
%         [source_points,idx] = pcdownsample(temp_source_points,'random',0.5);
        [source_points,idx]  = pcdownsample(temp_source_points,'nonuniformGridSample',6);
        target_points=select(temp_target_points,idx);
        graph.edge{edgeNum}.orgin_pair_points_1 = source_points;
        graph.edge{edgeNum}.orgin_pair_points_2 = target_points;
        graph.edge{edgeNum}.relative_pair_points_2 = target_points;
        if source_points(:,1)==target_points(:,1)
            fprintf('test');
        end
        edgeNum=edgeNum+1;
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

% idx_p=graph.edge{2}.idx(1);
% idx_q=graph.edge{2}.idx(2);
% AreaPt1=importdata(graph.node{idx_p}.filePath);
% AreaPt2=importdata(graph.node{idx_q}.filePath);
% AreaPt1=AreaPt1(:,1:3);
% AreaPt2=AreaPt2(:,1:3);
% hold on
% axis off
% AreaPt1=[AreaPt1,ones(length(AreaPt1(:,1)), 1)]';
% AreaPt2=[AreaPt2,ones(length(AreaPt2(:,1)), 1)]';
% plot3(AreaPt1(1,:),AreaPt1(2,:),AreaPt1(3,:),'.','MarkerSize',0.5);
% plot3(AreaPt2(1,:),AreaPt2(2,:),AreaPt2(3,:),'.','MarkerSize',0.5);
% [pt1,pt2]=FilterPoints(AreaPt1,AreaPt2);
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
