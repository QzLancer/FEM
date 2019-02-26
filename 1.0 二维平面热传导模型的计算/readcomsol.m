function [Coordinate,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol(filename)
%��ȡcomsol�еķ������ݣ��ظ�һ��poofee������
%2018/12/13
%��Ϥmatlab���ļ�IO
%by QzLancer
%------------------��ȡ�ļ��������ļ���ʶ���ʹ�����Ϣ
[fileID, Errmessage] = fopen(filename, 'r');
if fileID == -1
    disp(Errmessage);
end
%------------------���˵���Ч��Ϣ
for i = 1:18
    fgetl(fileID);
end
%------------------��ȡ�ڵ����Ŀ������
NumOfMeshpoint = fscanf(fileID, '%d # number of mesh points\n', [1, 1]);
for i = 1:3
    fgetl(fileID);
end
Coordinate = fscanf(fileID, '%f %f\n', [2, NumOfMeshpoint]);
Coordinate = Coordinate';
%------------------��ȡ���㵥Ԫ����Ŀ,�ڵ�,����ʵ��ָ�� 
for i=1:8
    fgetl(fileID);
end
NumOfVtxElement = fscanf(fileID, '%d # number of elements\n', [1, 1]);
fgetl(fileID);
VtxElement = fscanf(fileID, '%d\n', [1, NumOfVtxElement]);
VtxElement = VtxElement' + 1;
NumOfVtxEntity = fscanf(fileID, '%d # number of geometric entity indices\n', [1, 1]);
fgetl(fileID);
VtxEntity = fscanf(fileID, '%d\n', [1, NumOfVtxEntity]);
VtxEntity = VtxEntity';
%------------------��ȡ�߽絥Ԫ����Ŀ���ڵ㣬����ʵ��ָ��
for i=1:6
    fgetl(fileID);
end
NumOfEdgElement = fscanf(fileID, '%d # number of elements\n', [1, 1]);
fgetl(fileID);
EdgElement = fscanf(fileID, '%d\n', [2, NumOfEdgElement]);
EdgElement = EdgElement' + 1;
NumOfEdgEntity = fscanf(fileID, '%d # number of geometric entity indices\n', [1, 1]);
fgetl(fileID);
EdgEntity = fscanf(fileID, '%d\n', [1, NumOfEdgEntity]);
EdgEntity = EdgEntity';
%------------------��ȡ�����ε�Ԫ����Ŀ���ڵ㣬����ʵ��ָ��
for i=1:6
    fgetl(fileID);
end
NumOfTriElement = fscanf(fileID, '%d # number of elements\n', [1, 1]);
fgetl(fileID);
TriElement = fscanf(fileID, '%d\n', [3, NumOfTriElement]);
TriElement = TriElement' + 1;
NumOfTriEntity = fscanf(fileID, '%d # number of geometric entity indices\n', [1, 1]);
fgetl(fileID);
TriEntity = fscanf(fileID, '%d\n', [1, NumOfTriEntity]);
TriEntity = TriEntity';
%------------------��ȡ��ϣ��ر��ļ�
fclose(fileID);