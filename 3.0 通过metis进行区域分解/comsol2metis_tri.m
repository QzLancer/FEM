%��COMSOL�Ķ�ά��������ת����metis��ʶ��ķ�������
%���.mpmetis�ļ�
%�߽絥ԪӦ����δ�����ܱ����ֵ����ʵ�����
%By QzLancer
%2019/4/3

%-------------------------------------------ԭʼ�ļ���ȡ
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh_source.mphtxt');
%-------------------------------------------��TriElementд���µ��ļ�
fp=fopen('mesh_source.mpmetis','w');
a=length(TriElement);
fprintf(fp,'%d\n',a);
fprintf(fp,'%d %d %d\n',TriElement');