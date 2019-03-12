%ͨ��ֱ�ӷ�����άƽ����̬�ȳ���֮ǰ��Գ���̬�¶ȳ�����������⣬���Դ���������ҳ�����
%����һ���������⣬������һ��߽�����,�ڶ���߽����������ᴦ����
%���ϵ��ȵ���Ϊ52W/(m.K),�߽���¶�Ϊ273.15K�������2�ĵ�λ�������Ϊ1000W/m^3
%�ο����������������̵�ų���ֵ������ P57 ��Գƴų�����Ԫ��
%ƽ�泡�����������
%2018/12/12
%By QzLancer
%-------------------------------��ȡ�����ļ�
clear all;
[Coor,VtxElement,VtxEntity,EdgElement,EdgEntity,TriElement,TriEntity] = readcomsol('mesh.mphtxt');
%-------------------------------���ݲ�ֵҪ��������Ԫ�ļ��β�����Ŀǰֻ�����ڲ���Ԫ
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
% TriRadius = 1.5./(1./(TriR(:,1)+TriR(:,2))+1./(TriR(:,2)+TriR(:,3))+1./(TriR(:,3)+TriR(:,1)));
%-------------------------------������غ��ȵ���
Cond = 52; 
%��Ч�ȵ���
% EquCond = Cond./TriRadius;
%�ҳ���Դ���ڵ�Ԫ
Domain2 = find(TriEntity==2);
%�����е�Ԫ������Դ
Source = zeros(length(TriElement),1);
Source(Domain2) = 10000000;
%-------------------------------��Ԫ����������ϳ�
S = zeros(length(Coor));
F = zeros(length(Coor),1);
for k = 1:length(TriElement)
    for i = 1:3
        for j = 1:3
%             Se= EquCond(k) * (r(k,i)*r(k,j) + q(k,i)*q(k,j)) / (4*Area(k));
            Se= Cond .* (r(k,i)*r(k,j) + q(k,i)*q(k,j)) / (4*Area(k));
            S(TriElement(k,i),TriElement(k,j)) = S(TriElement(k,i),TriElement(k,j)) + Se;
        end
        Fe = Source(k)*Area(k)/3;
        F(TriElement(k,i)) = F(TriElement(k,i)) + Fe;
    end
end
%-------------------------------����ֱ�ӷ����
%�������߽��
% EquTemp = zeros(length(Coor),1);
Temp = zeros(length(Coor),1);
Boundary = find(Z==0 | Z==0.14 | R==0.1);
FreeNodes = find(~(Z==0 | Z==0.14 | R==0.1));

% EquTemp(Boundary) = 273.15.*R(Boundary);
Temp(Boundary) = 293.15;
% F1 = S(FreeNodes,:)*EquTemp;
F1 = S(FreeNodes,:)*Temp;
F2 = F(FreeNodes) - F1;
% EquTemp(FreeNodes) = S(FreeNodes,FreeNodes)\F2;
% Temp = EquTemp./R;
Temp(FreeNodes) = S(FreeNodes,FreeNodes)\F2;
%-------------------------------����
Interp1 = scatteredInterpolant(R,Z,Temp);
tx = 0.02:1e-3:0.1;
ty = 0:1e-3:0.14;
[qx,qy] = meshgrid(tx,ty);
qz = Interp1(qx,qy);
% mesh(qx,qy,qz);
subplot(1,2,2);
contourf(qx,qy,qz,20);colorbar;
hold on;
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
    plot([EdgCoor1(i, 1), EdgCoor2(i, 1)], [EdgCoor1(i, 2), EdgCoor2(i, 2)],'Color','k','Linewidth',1);
    hold on;
end
%-------------------------��ȡCOMSOL��������Ľ�������к���
[fileID, Errmessage] = fopen('comsolsolution.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
comsoldata = fscanf(fileID,'%lf %lf %lf\n',[3,length(Coor)]);
comsoldata = comsoldata';
fclose(fileID);
Interp2 = scatteredInterpolant(comsoldata(:,1),comsoldata(:,2),comsoldata(:,3));
qz = Interp2(qx,qy);
subplot(1,2,1);
contourf(qx,qy,qz,20);colorbar;
