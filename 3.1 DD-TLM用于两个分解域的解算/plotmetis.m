%����metis������������ͼ�۲죬��Ҫ�ǣ������Ľ��紦���ܲ���ͨ��metis�����������
%By QzLancer
%2019/4/9

%---------------------------------��ȡ�ļ�
close all;
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
[fileID] = fopen('mesh_source.mpmetis.epart.2');
eDomain = fscanf(fileID,'%d\n');
[fileID] = fopen('mesh_source.mpmetis.npart.2');
nDomain = fscanf(fileID,'%d\n');
%---------------------------------��ȡ��Ԫ����������
Domain0Element = TriElement(eDomain==0,:);
Domain1Element = TriElement(eDomain==1,:);
Domain0x(:,1) = Coor(Domain0Element(:,1),1);
Domain0x(:,2) = Coor(Domain0Element(:,2),1);
Domain0x(:,3) = Coor(Domain0Element(:,3),1);
Domain0y(:,1) = Coor(Domain0Element(:,1),2);
Domain0y(:,2) = Coor(Domain0Element(:,2),2);
Domain0y(:,3) = Coor(Domain0Element(:,3),2);
patch(Domain0x',Domain0y','red','FaceAlpha',.3);
hold on;
Domain1x(:,1) = Coor(Domain1Element(:,1),1);
Domain1x(:,2) = Coor(Domain1Element(:,2),1);
Domain1x(:,3) = Coor(Domain1Element(:,3),1);
Domain1y(:,1) = Coor(Domain1Element(:,1),2);
Domain1y(:,2) = Coor(Domain1Element(:,2),2);
Domain1y(:,3) = Coor(Domain1Element(:,3),2);
patch(Domain1x',Domain1y','blue','FaceAlpha',.3);
hold on;
%---------------------------------��ȡ�ڵ����������
nDomain0Coor = Coor(nDomain==0,:);
nDomain1Coor = Coor(nDomain==1,:);
plot(nDomain0Coor(:,1),nDomain0Coor(:,2),'.b');
hold on;
axis equal;
plot(nDomain1Coor(:,1),nDomain1Coor(:,2),'.r');
hold on;
%---------------------------------�ҳ������߽��ĵ�Ԫ
%�������е�Ԫ���ҳ���Щ��Ԫ�����˲�ֹһ����Ľڵ�
eleDomainNode = nDomain(TriElement);
j = 1;
k = 1
for i = 1:length(TriElement)
    if length(unique(eleDomainNode(i,:))) ~= 1
        if eDomain(i) == 0
            Bound0Element(j,1) = i;
            j = j+1;
        else
            Bound1Element(k,1) = i;
            k = k+1;
        end
    end
end
%��ȡ�����Ƴ���Щ��Ԫ
DomainB0Node = TriElement(Bound0Element,:);
eDomainB0x(:,1) = Coor(DomainB0Node(:,1),1);
eDomainB0x(:,2) = Coor(DomainB0Node(:,2),1);
eDomainB0x(:,3) = Coor(DomainB0Node(:,3),1);
eDomainB0y(:,1) = Coor(DomainB0Node(:,1),2);
eDomainB0y(:,2) = Coor(DomainB0Node(:,2),2);
eDomainB0y(:,3) = Coor(DomainB0Node(:,3),2);
% patch(eDomainB0x',eDomainB0y','blue','FaceAlpha',.3);
patch(eDomainB0x',eDomainB0y','white');
hold on;
DomainB1Node = TriElement(Bound1Element,:);
eDomainB1x(:,1) = Coor(DomainB1Node(:,1),1);
eDomainB1x(:,2) = Coor(DomainB1Node(:,2),1);
eDomainB1x(:,3) = Coor(DomainB1Node(:,3),1);
eDomainB1y(:,1) = Coor(DomainB1Node(:,1),2);
eDomainB1y(:,2) = Coor(DomainB1Node(:,2),2);
eDomainB1y(:,3) = Coor(DomainB1Node(:,3),2);
% patch(eDomainB1x',eDomainB1y','red','FaceAlpha',.3);
patch(eDomainB1x',eDomainB1y','white');
hold on;
% �ҳ����ظ������Ľڵ�
Domain0BoundNode = TriElement(Bound0Element,:);
Domain1BoundNode = TriElement(Bound1Element,:);
intersect(Domain0BoundNode,Domain1BoundNode);
%-----------------------------------�ҳ��߽�ڵ�
%����ȫ����Ԫ��ȫ���ڵ㣬������ֽڵ��������뵥Ԫ������һ�£����Ϊ���紦�ڵ�
eleDomainNode = nDomain(TriElement);
k = 1; 
for i = 1:length(eDomain)
    for j = 1:3
        if eleDomainNode(i,j) ~= eDomain(i)
            DomainBNode(k,1) = TriElement(i,j);
            k = k+1;
        end
    end
end
plot(Coor(DomainBNode,1),Coor(DomainBNode,2),'xg');