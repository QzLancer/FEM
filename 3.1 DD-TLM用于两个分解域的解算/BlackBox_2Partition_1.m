%通过传输线来实现被分解的两个区域的信息交换
%BlackBox-TLM 区域内通过直接法求解，看看是否能够收敛
%Domain为COMSOL定义的区域 Partition为metis分解得到的区域
%如果是带边界条件的情况，线单元应该怎么处理？
%把所有的单元和节点分离，各个区域的单元重新编号
%By QzLancer
%2019/4/9
%-----------------------------读取文件
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
[fileID] = fopen('mesh_source.mpmetis.epart.2');
ePartTable = fscanf(fileID,'%d\n');
ePartTable = ePartTable+1;
[fileID] = fopen('mesh_source.mpmetis.npart.2');
nPartTable = fscanf(fileID,'%d\n');
nPartTable = nPartTable+1;
Y0 = 1;
%------------------------------得到每个Part包含的Element和Node编号
for i=1:2
    PartElementNum(:,i) = find(ePartTable==i);
end
% Part1ElementNum = find(ePartTable==1);
% Part2ElementNum = find(ePartTable==2);
Part1NodeNum = find(nPartTable==1);
Part2NodeNum = find(nPartTable==2);
% -----------------------------读取交界处的节点
PartEleNodeTable = nPartTable(TriElement);
k = 1; 
for i = 1:length(ePartTable)
    for j = 1:3
        if PartEleNodeTable(i,j) ~= ePartTable(i)
            JunNodeNum(k,1) = TriElement(i,j);
            k = k+1;
        end
    end
end
JunNodeNum = unique(JunNodeNum);
%------------------------------添加交界处编号，并且获得交界处节点在每个Part中的编号
Part1NodeNum = unique([Part1NodeNum;JunNodeNum]);
Part2NodeNum = unique([Part2NodeNum;JunNodeNum]);
PartNodeNum = [Part1NodeNum,Part2NodeNum];
PartJunNodeNum = zeros(length(JunNodeNum),2);
% Part1JunNodeNum = zeros(length(JunNodeNum),1);
% Part2JunNodeNum = zeros(length(JunNodeNum),1);
for i=1:length(JunNodeNum)
%     Part1JunNodeNum(i) = find(Part1NodeNum==JunNodeNum(i));
%     Part2JunNodeNum(i) = find(Part2NodeNum==JunNodeNum(i));
    for j=1:2
        PartJunNodeNum(i,j) = find(PartNodeNum(:,j)==JunNodeNum(i));
    end
end
%------------------------------每个Part中的Element中的Node重新编号为该Part中的节点编号
PartElement = zeros(length(PartElementNum),3,2);
for k = 1:2
    PartElement(:,:,k) = TriElement(PartElementNum(:,k),:); 
    for j = 1:3
        for i = 1:length(PartElementNum)
            PartElement(i,j,k) = find(PartNodeNum(:,k) == PartElement(i,j,k));
        end
    end
end
%------------------------------每个单元的节点坐标
PartCoor = zeros(length(PartNodeNum),2,2);
for i = 1:2
    PartCoor(:,:,i) = Coor(PartNodeNum(:,i),:);
end
%------------------------------绘图验证单元和节点的分离和重新编号是否正确
% plot(PartCoor(:,1,1),PartCoor(:,2,1),'.b');
% axis equal;  
% hold on;
% plot(PartCoor(:,1,2),PartCoor(:,2,2),'.r');
% hold on
% for j = 1:2
%     for i = 1:3
%         PartElementCoorX(:,i,j) = PartCoor(PartElement(:,i,j),1,j);
%         PartElementCoorY(:,i,j) = PartCoor(PartElement(:,i,j),2,j);
%     end
% end
% patch(PartElementCoorX(:,:,1)',PartElementCoorY(:,:,1)','blue','FaceAlpha',.3);
% hold on;
% patch(PartElementCoorX(:,:,2)',PartElementCoorY(:,:,2)','red','FaceAlpha',.3);

%------------------------------每个单元需要的基本几何参数
%-----全局单元几何参数求解
%导出所有三角形单元RZ轴的坐标
R = Coor(:,1);
Z = Coor(:,2);
TriR = R(TriElement);
TriZ = Z(TriElement);
%计算出所有单元的p,q,r和Area
p(:,1) = TriR(:,2).*TriZ(:,3) - TriZ(:,2).*TriR(:,3);
p(:,2) = TriR(:,3).*TriZ(:,1) - TriZ(:,3).*TriR(:,1);
p(:,3) = TriR(:,1).*TriZ(:,2) - TriZ(:,1).*TriR(:,2);
q(:,1) = TriZ(:,2) - TriZ(:,3);
q(:,2) = TriZ(:,3) - TriZ(:,1);
q(:,3) = TriZ(:,1) - TriZ(:,2);
r(:,1) = TriR(:,3) - TriR(:,2);
r(:,2) = TriR(:,1) - TriR(:,3);
r(:,3) = TriR(:,2) - TriR(:,1);
Area = (q(:,1).*r(:,2) - q(:,2).*r(:,1))/2;
% 三角形所有单元的重心处的半径
TriRadius = (TriR(:,1)+TriR(:,2)+TriR(:,3))./3;
%-----每个Part的几何参数
for i = 1:2
    PartR = R(PartNodeNum);
    PartTriR(:,:,i) = TriR(PartElementNum(:,i),:);
    Partp(:,:,i) = p(PartElementNum(:,i),:);
    Partq(:,:,i) = q(PartElementNum(:,i),:);
    Partr(:,:,i) = r(PartElementNum(:,i),:);
    PartArea = Area(PartElementNum);
    PartTriRadius = TriRadius(PartElementNum);
end
%------------------------------求出每个Part的负载和热导率单元
Cond = 52;
SourceElement = find(TriEntity==2);
SourceTable = ePartTable(SourceElement);
PartSource = zeros(length(PartElementNum),2);
for i = 1:2
    PartSourceElementNum{i} = SourceElement(SourceTable==i);
    for j = 1:length(PartSourceElementNum{i})
        PartSourceElementNum{i}(j) = find(PartElementNum(:,i) == PartSourceElementNum{i}(j));
    end
    PartSource(PartSourceElementNum{i},i) = 10000000; 
end
%------------------------------每个Part的单元分析和总体合成
PartS = zeros(length(PartNodeNum),length(PartNodeNum),2);
PartF = zeros(length(PartNodeNum),2);
for m = 1:2
    for k = 1:length(PartElementNum)
        for i = 1:3
            for j = 1:3
                Se(i,j,k)= (pi*Cond*PartTriRadius(k,m)*(Partr(k,i,m)*Partr(k,j,m) + Partq(k,i,m)*Partq(k,j,m)))/(2*PartArea(k,m));
                PartS(PartElement(k,i,m),PartElement(k,j,m),m) = PartS(PartElement(k,i,m),PartElement(k,j,m),m) + Se(i,j,k);
            end
        end
    end
end
%------------------------------TLM迭代过程，每个Part用直接法求解入射过程