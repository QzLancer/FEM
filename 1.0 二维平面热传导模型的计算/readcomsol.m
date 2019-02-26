function [Coordinate,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol(filename)
%读取comsol中的分网数据，重复一遍poofee的轮子
%2018/12/13
%熟悉matlab的文件IO
%by QzLancer
%------------------读取文件，返回文件标识符和错误信息
[fileID, Errmessage] = fopen(filename, 'r');
if fileID == -1
    disp(Errmessage);
end
%------------------过滤掉无效信息
for i = 1:18
    fgetl(fileID);
end
%------------------读取节点的数目和坐标
NumOfMeshpoint = fscanf(fileID, '%d # number of mesh points\n', [1, 1]);
for i = 1:3
    fgetl(fileID);
end
Coordinate = fscanf(fileID, '%f %f\n', [2, NumOfMeshpoint]);
Coordinate = Coordinate';
%------------------读取顶点单元的数目,节点,几何实体指数 
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
%------------------读取边界单元的数目，节点，几何实体指数
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
%------------------读取三角形单元的数目，节点，几何实体指数
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
%------------------读取完毕，关闭文件
fclose(fileID);