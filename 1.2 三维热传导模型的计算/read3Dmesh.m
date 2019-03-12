function [Coordinate,VtxElement,VtxEntity,TriElement,TriEntity,TetElement,TetEntity] = readcomsol(filename)
%读取三维模型的分网数据，基于之前的二维分网数据读取程序
%2019/2/28
%matlab的文件IO已经忘得差不多了（逃
%by QzLancer
%------------------读取文件，返回文件标识符和错误信息
% filename = '3D_rect_mesh.mphtxt';
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
Coordinate = fscanf(fileID, '%f %f\n', [3, NumOfMeshpoint]);
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
VtxEntity = VtxEntity' + 1;
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
TriEntity = TriEntity' + 1;
%------------------读取四面体单元的数目，节点，几何实体指数
for i=1:6
    fgetl(fileID);
end
NumOfTetElement = fscanf(fileID, '%d # number of elements\n', [1, 1]);
fgetl(fileID);
TetElement = fscanf(fileID, '%d\n', [4, NumOfTetElement]);
TetElement = TetElement' + 1;
NumOfTetEntity = fscanf(fileID, '%d # number of geometric entity indices\n', [1, 1]);
fgetl(fileID);
TetEntity = fscanf(fileID, '%d\n', [1, NumOfTetEntity]);
TetEntity = TetEntity' + 1;
%------------------读取完毕，关闭文件
fclose(fileID);
