%ͨ����������ʵ�ֱ��ֽ�������������Ϣ����
%������ͨ��ֱ�ӷ���⣬�����Ƿ��ܹ�����
%DomainΪCOMSOL��������� PartitionΪmetis�ֽ�õ�������
%����Ǵ��߽�������������ߵ�ԪӦ����ô����
%�����еĵ�Ԫ�ͽڵ���룬��������ĵ�Ԫ���±��
%By QzLancer
%2019/4/9
%-----------------------------��ȡ�ļ�
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
e = 0.1;
Temp = 273.15;
%------------------------------�õ�ÿ��Part������Element��Node���
for i=1:2
    PartElementNum(:,i) = find(ePartTable==i);
end
% Part1ElementNum = find(ePartTable==1);
% Part2ElementNum = find(ePartTable==2);
Part1NodeNum = find(nPartTable==1);
Part2NodeNum = find(nPartTable==2);
% -----------------------------��ȡ���紦�Ľڵ�
EleNodeTable = nPartTable(TriElement);
k = 1; 
for i = 1:length(ePartTable)
    for j = 1:3
        if EleNodeTable(i,j) ~= ePartTable(i)
            JunNodeNum(k,1) = TriElement(i,j);
            k = k+1;
        end
    end
end
JunNodeNum = unique(JunNodeNum);
%------------------------------��ӽ��紦��ţ����һ�ý��紦�ڵ���ÿ��Part�еı��
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
%------------------------------ÿ��Part�е�Element�е�Node���±��Ϊ��Part�еĽڵ���
PartElement = zeros(length(PartElementNum),3,2);
for k = 1:2
    PartElement(:,:,k) = TriElement(PartElementNum(:,k),:); 
    for j = 1:3
        for i = 1:length(PartElementNum)
            PartElement(i,j,k) = find(PartNodeNum(:,k) == PartElement(i,j,k));
        end
    end
end
%------------------------------ÿ����Ԫ�Ľڵ�����
PartCoor = zeros(length(PartNodeNum),2,2);
for i = 1:2
    PartCoor(:,:,i) = Coor(PartNodeNum(:,i),:);
end
%------------------------------��ͼ��֤��Ԫ�ͽڵ�ķ�������±���Ƿ���ȷ
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

%------------------------------ÿ����Ԫ��Ҫ�Ļ������β���
%-----ȫ�ֵ�Ԫ���β������
%�������������ε�ԪRZ�������
R = Coor(:,1);
Z = Coor(:,2);
TriR = R(TriElement);
TriZ = Z(TriElement);
%��������е�Ԫ��p,q,r��Area
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
% ���������е�Ԫ�����Ĵ��İ뾶
TriRadius = (TriR(:,1)+TriR(:,2)+TriR(:,3))./3;
%-----ÿ��Part�ļ��β���
for i = 1:2
    PartR = R(PartNodeNum);
    PartZ = Z(PartNodeNum);
    PartTriR(:,:,i) = TriR(PartElementNum(:,i),:);
    Partp(:,:,i) = p(PartElementNum(:,i),:);
    Partq(:,:,i) = q(PartElementNum(:,i),:);
    Partr(:,:,i) = r(PartElementNum(:,i),:);
    PartArea = Area(PartElementNum);
    PartTriRadius = TriRadius(PartElementNum);
end
%------------------------------���ÿ��Part�ĸ��غ��ȵ��ʵ�Ԫ
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
%------------------------------ÿ��Part�ĵ�Ԫ����������ϳ�
PartS = zeros(length(PartNodeNum),length(PartNodeNum),2);
PartF = zeros(length(PartNodeNum),2);
for m = 1:2
    for k = 1:length(PartElementNum)
        for i = 1:3
            for j = 1:3
                Se = (pi*Cond*PartTriRadius(k,m)*(Partr(k,i,m)*Partr(k,j,m) + Partq(k,i,m)*Partq(k,j,m)))/(2*PartArea(k,m));
                PartS(PartElement(k,i,m),PartElement(k,j,m),m) = PartS(PartElement(k,i,m),PartElement(k,j,m),m) + Se;
            end
            Fe = pi*PartSource(k,m)*PartArea(k,m)*(PartR(PartElement(k,i,m),m)+3*PartTriRadius(k,m))/6;
            PartF(PartElement(k,i,m),m) = PartF(PartElement(k,i,m),m) + Fe;
        end
    end
end

%------------------------------TLM�������̣�ÿ��Part��ֱ�ӷ�����������
%-----�������߽�
for m = 1:2
    PartBoundary{m} = find(PartZ(:,m)==0 | PartZ(:,m)==0.14 | PartR(:,m)==0.1);
    PartFreeNodes{m} = find(~(PartZ(:,m)==0 | PartZ(:,m)==0.14 | PartR(:,m)==0.1));
end
%-----����ͷ������
PartY = zeros(length(PartS),length(PartS),2);
PartI = zeros(length(PartF),2);
PartVa = zeros(length(PartJunNodeNum),2);
PartVc = zeros(length(PartJunNodeNum),1);
PartVi = zeros(length(PartJunNodeNum),2);
PartVr = zeros(length(PartJunNodeNum),2);
for m = 1:2
    for i = 1:length(PartJunNodeNum)
        PartY(PartJunNodeNum(i,m),PartJunNodeNum(i,m),m) = Y0;%����Ҫ�Խ���һ������ֵ
    end
    PartI(PartJunNodeNum(:,m),m) = 2*PartVi(:,m)*Y0;
end
%��һ������������
PartSi = PartS+PartY;
PartFi = PartF+PartI;
PartV = zeros(length(PartNodeNum),2);
for m = 1:2
    PartV(PartBoundary{m},m) = Temp;
    PartFi(PartFreeNodes{m},m) = PartFi(PartFreeNodes{m},m)-PartSi(PartFreeNodes{m},:,m)*PartV(:,m);
    PartV(PartFreeNodes{m},m) = PartSi(PartFreeNodes{m},PartFreeNodes{m},m)\PartFi(PartFreeNodes{m},m);
    PartVa(:,m) = PartV(PartJunNodeNum(:,m),m);
end
%��ʼ����
while norm(PartVa(:,1)-PartVa(:,2))>e
    PartVr = PartVa - PartVi;
    %����������
    PartVc = PartVr(:,1)+PartVr(:,2);
    PartVi(:,1) = PartVc-PartVr(:,1);
    PartVi(:,2) = PartVc-PartVr(:,2);
    %����������
    for m = 1:2
        PartI(PartJunNodeNum(:,m),m) = 2*PartVi(:,m)*Y0;
    end
    PartFi = PartF+PartI;
    PartV = zeros(length(PartNodeNum),2);
    for m = 1:2
        PartV(PartBoundary{m},m) = Temp;
        PartFi(PartFreeNodes{m},m) = PartFi(PartFreeNodes{m},m)-PartSi(PartFreeNodes{m},:,m)*PartV(:,m);
        PartV(PartFreeNodes{m},m) = PartSi(PartFreeNodes{m},PartFreeNodes{m},m)\PartFi(PartFreeNodes{m},m);
        PartVa(:,m) = PartV(PartJunNodeNum(:,m),m);
    end
end
%------------------------------����
V = zeros(length(Coor),1);
for i = 1:length(PartNodeNum)
    for m = 1:2
        V(PartNodeNum(i,m)) = PartV(i,m);
    end
end
Interp1 = scatteredInterpolant(R,Z,V);
tx = 0.02:1e-3:0.1;
ty = 0:1e-3:0.14;
[qx,qy] = meshgrid(tx,ty);
qz = Interp1(qx,qy);
subplot(1,2,2);
contourf(qx,qy,qz,20);colorbar;