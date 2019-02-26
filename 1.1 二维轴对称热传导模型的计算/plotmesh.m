%������֪�ķ������ݻ��Ʒ���
%2018/12/9
%By QzLancer
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh.mphtxt');
plot([0 ,0], [-0.02, 0.16], 'm--');
xlim([-0.025, 0.15]);
ylim([-0.02, 0.16]);
hold on;
%-----------------------��ʶ�����㵥Ԫ������
VtxCoor = Coor(VtxElement, :);
plot(VtxCoor(:, 1), VtxCoor(:, 2), 'x');
hold on;
%-----------------------��ʶ���߽絥Ԫ������
EdgNode1 = EdgElement(:, 1);
EdgNode2 = EdgElement(:, 2);
EdgCoor1 = Coor(EdgNode1, :);
EdgCoor2 = Coor(EdgNode2, :);
% hold on;
% plot(EdgCoor1(:, 1), EdgCoor1(:, 2), 'x');
% hold on;
% plot(EdgCoor2(:, 1), EdgCoor2(:, 2), 'x');
Edglength = length(EdgElement); 
for i = 1:Edglength
    plot([EdgCoor1(i, 1), EdgCoor2(i, 1)], [EdgCoor1(i, 2), EdgCoor2(i, 2)]);
    hold on;
end
%-----------------------��ʶ�������ε�Ԫ������
TriNode1 = TriElement(:, 1);
TriNode2 = TriElement(:, 2);
TriNode3 = TriElement(:, 3);
TriCoor1 = Coor(TriNode1, :);
TriCoor2 = Coor(TriNode2, :);
TriCoor3 = Coor(TriNode3, :);
Trilength = length(TriElement);
for i = 1:Trilength
    plot([TriCoor1(i, 1), TriCoor2(i, 1)], [TriCoor1(i, 2), TriCoor2(i, 2)]);
    hold on;
    plot([TriCoor2(i, 1), TriCoor3(i, 1)], [TriCoor2(i, 2), TriCoor3(i, 2)]);
    hold on;
    plot([TriCoor3(i, 1), TriCoor1(i, 1)], [TriCoor3(i, 2), TriCoor1(i, 2)]);
    hold on;
end