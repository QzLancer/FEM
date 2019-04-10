%ͨ����������ʵ�ֱ��ֽ�������������Ϣ����
%BlackBox-TLM ������ͨ��ֱ�ӷ���⣬�����Ƿ��ܹ�����
%DomainΪCOMSOL��������� PartitionΪmetis�ֽ�õ�������
%����Ǵ��߽�������������ߵ�ԪӦ����ô����
%���������޷����㣬��ʧ��
%By QzLancer
%2019/4/9
%-----------------------------��ȡ�ļ�
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
[fileID] = fopen('mesh_source.mpmetis.epart.2');
ePartition = fscanf(fileID,'%d\n');
[fileID] = fopen('mesh_source.mpmetis.npart.2');
nPartition = fscanf(fileID,'%d\n');
Y0 = 1;
%-----------------------------������Ԫ�ļ��β���
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
%��~R^2
RR = zeros(length(TriElement), 1);
for i = 1:length(TriElement)
    RR(i) = TriR(i, :)*([2,1,1;1,2,1;1,1,2]*TriR(i, :)') / 12;
end
%������á���Գ��ȴ�������Ԫ��ʽ���е���ν��ȷ��ʽ����Ҫ�����r�йص�ϵ��
cofF1(:, 1) = (3*TriR(:, 1).^2 + TriR(:, 2).^2 + TriR(:, 3).^2 + ...
               2*TriR(:, 1).*TriR(:, 2) + TriR(:, 2).*TriR(:, 3) + 2*TriR(: ,3).*TriR(:, 1));
cofF1(:, 2) = (TriR(:, 1).^2 + 3*TriR(:, 2).^2 + TriR(:, 3).^2 + ...
               2*TriR(:, 1).*TriR(:, 2) + 2*TriR(:, 2).*TriR(:, 3) + TriR(: ,3).*TriR(:, 1));
cofF1(:, 3) = (TriR(:, 1).^2 + TriR(:, 2).^2 + 3*TriR(:, 3).^2 + ...
               TriR(:, 1).*TriR(:, 2) + 2*TriR(:, 2).*TriR(:, 3) + 2*TriR(: ,3).*TriR(:, 1));
%-------------------------------������غ��ȵ���
Cond = 52;
%��Ч�ȵ���
% EquCond = Cond./TriRadius;
%�ҳ���Դ���ڵ�Ԫ
Domain = find(TriEntity==2);
%�����е�Ԫ������Դ
Source0 = zeros(length(TriElement),1);
Source1 = zeros(length(TriElement),1);
DomainTable = ePartition(Domain);
DomainPart0 = Domain(DomainTable==0);
DomainPart1 = Domain(DomainTable==1);
Source0(DomainPart0) = 10000000;
Source1(DomainPart1) = 10000000;
% -----------------------------��ȡ���紦�Ľڵ�
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
%-----------------------------��������ķ��룬���紦��
Part0ElementTable = find(ePartition==0);
Part1ElementTable = find(ePartition==1);
Part0Element = TriElement(Part0ElementTable,:);
Part1Element = TriElement(Part1ElementTable,:);
%-----------------------------Partition0�浥Ԫ������װ��
S0 = zeros(length(Coor));
F0 = zeros(length(Coor),1);
p0 = p(Part0ElementTable,:);
q0 = q(Part0ElementTable,:);
r0 = r(Part0ElementTable,:);
Area0 = Area(Part0ElementTable,:);
TriRadius0 = TriRadius(Part0ElementTable,:);
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
%-----------------------------Partition1�浥Ԫ������װ��
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
%-----------------------------���������ͨ�������ߵ������
%Partiton0�������
Vr0 = zeros(length(S0),1);
Vi0 = zeros(length(S0),1);
Va0 = zeros(length(S0),1);
Vc0 = zeros(length(S0),1);
Temp = zeros(length(Coor),1);
Boundary = find(Z==0 | Z==0.14 | R==0.1);
FreeNodes = find(~(Z==0 | Z==0.14 | R==0.1));
Temp(Boundary) = 273.15;
