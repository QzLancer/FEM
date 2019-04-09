%通过传输线来实现被分解的两个区域的信息交换
%BlackBox-TLM 区域内通过直接法求解，看看是否能够收敛
%Domain为COMSOL定义的区域 Partition为metis分解得到的区域
%如果是带边界条件的情况，线单元应该怎么处理？
%By QzLancer
%2019/4/9
%-----------------------------读取文件
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
[fileID] = fopen('mesh_source.mpmetis.epart.2');
ePartition = fscanf(fileID,'%d\n');
[fileID] = fopen('mesh_source.mpmetis.npart.2');
nPartition = fscanf(fileID,'%d\n');
Y0 = 1;
%-----------------------------求解各单元的几何参数
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
%求~R^2
RR = zeros(length(TriElement), 1);
for i = 1:length(TriElement)
    RR(i) = TriR(i, :)*([2,1,1;1,2,1;1,1,2]*TriR(i, :)') / 12;
end
%如果采用《轴对称热传导有限元格式》中的所谓正确格式，需要计算和r有关的系数
cofF1(:, 1) = (3*TriR(:, 1).^2 + TriR(:, 2).^2 + TriR(:, 3).^2 + ...
               2*TriR(:, 1).*TriR(:, 2) + TriR(:, 2).*TriR(:, 3) + 2*TriR(: ,3).*TriR(:, 1));
cofF1(:, 2) = (TriR(:, 1).^2 + 3*TriR(:, 2).^2 + TriR(:, 3).^2 + ...
               2*TriR(:, 1).*TriR(:, 2) + 2*TriR(:, 2).*TriR(:, 3) + TriR(: ,3).*TriR(:, 1));
cofF1(:, 3) = (TriR(:, 1).^2 + TriR(:, 2).^2 + 3*TriR(:, 3).^2 + ...
               TriR(:, 1).*TriR(:, 2) + 2*TriR(:, 2).*TriR(:, 3) + 2*TriR(: ,3).*TriR(:, 1));
%-------------------------------求出负载和热导率
Cond = 52;
%等效热导率
% EquCond = Cond./TriRadius;
%找出热源所在单元
Domain = find(TriEntity==2);
%给所有单元附上热源
Source0 = zeros(length(TriElement),1);
Source1 = zeros(length(TriElement),1);
DomainTable = ePartition(Domain);
DomainPart0 = Domain(DomainTable==0);
DomainPart1 = Domain(DomainTable==1);
Source0(DomainPart0) = 10000000;
Source1(DomainPart1) = 10000000;
% -----------------------------读取交界处的节点
eleNodePartition = nPartition(TriElement);
k = 1; 
for i = 1:length(ePartition)
    for j = 1:3
        if eleNodePartition(i,j) ~= ePartition(i)
            JunNode(k,1) = TriElement(i,j);
            k = k+1;
        end
    end
end
JunNode = unique(JunNode);
% %-----------------------------两个区域的分离，交界处理
Part0ElementTable = find(ePartition==0);
Part1ElementTable = find(ePartition==1);
Part0Element = TriElement(Part0ElementTable,:);
Part1Element = TriElement(Part1ElementTable,:);
% %-----------------------------Partition0面单元分析和装配
S0 = zeros(length(Coor));
F0 = zeros(length(Coor),1);
p0 = p(Part0ElementTable,:);
q0 = q(Part0ElementTable,:);
r0 = r(Part0ElementTable,:);
Area0 = Area(Part0ElementTable,:);
TriRadius0 = TriRadius(Part0ElementTable,:);
ST0 = zeros(length(Part0Element));
FT0 = zeros(length(Part0Element),1);
for i=1:length(JunNode)
    ST0(i,i) = Y0;
end
for k = 1:length(Part0Element)
    for i = 1:3
        for j = 1:3
            Se= (pi*Cond*TriRadius0(k)*(r0(k,i)*r0(k,j) + q0(k,i)*q0(k,j)))/(2*Area0(k));
            S0(Part0Element(k,i),Part0Element(k,j)) = S0(Part0Element(k,i),Part0Element(k,j)) + Se;
        end
        Fe = pi*Source0(Part0ElementTable(k))*Area0(k)*(R(Part0Element(k,i))+3*TriRadius0(k))/6;
        F0(Part0Element(k,i)) = F0(Part0Element(k,i)) + Fe;
    end
end
% %-----------------------------Partition1面单元分析和装配
S1 = zeros(length(Coor));
F1 = zeros(length(Coor),1);
p1 = p(Part1ElementTable,:);
q1 = q(Part1ElementTable,:);
r1 = r(Part1ElementTable,:);
Area1 = Area(Part1ElementTable,:);
TriRadius1 = TriRadius(Part1ElementTable,:);
for k = 1:length(Part1Element)
    for i = 1:3
        for j = 1:3
            Se= (pi*Cond*TriRadius1(k)*(r1(k,i)*r1(k,j) + q1(k,i)*q1(k,j)))/(2*Area1(k));
            S1(Part1Element(k,i),Part1Element(k,j)) = S1(Part1Element(k,i),Part1Element(k,j)) + Se;
        end
        Fe = pi*Source1(Part1ElementTable(k))*Area1(k)*(R(Part1Element(k,i))+3*TriRadius1(k))/6;
        F1(Part1Element(k,i)) = F1(Part1Element(k,i)) + Fe;
    end
end
S = S0+S1;