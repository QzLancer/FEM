%将COMSOL的二维分网数据转换成metis可识别的分网数据
%输出.mpmetis文件
%边界单元应当如何处理才能被划分到合适的区域？
%By QzLancer
%2019/4/3

%-------------------------------------------原始文件读取
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
%-------------------------------------------将TriElement写入新的文件
fp=fopen('mesh_source.mpmetis','w');
a=length(TriElement);
fprintf(fp,'%d\n',a);
fprintf(fp,'%d %d %d\n',TriElement');