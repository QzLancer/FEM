%通过直接法求解三维稳态热场
%包含第一类边界条件和第三类边界条件
%材料的热导率为50W/(m.k)，四面热绝缘，一面的热通量为5e3W/m2，一面的温度为273.15K
%参考《The finite element method in engineering》 Chapter16.4 P543-P545
%测试TriEntity与边界面的对应关系， TriEntity = 5时对应热通量面，TriEntity = 2对应恒温面
%其余面为热绝缘面
%（在读取分网程序中将Entity+1可以获得与COMSOL一致的编号）
%By QzLancer
%2019/3/1
%-------------------------------读取分网文件
clear all;
[Coor,VtxElement,VtxEntity,TriElement,TriEntity,TetElement,TetEntity] = read3Dmesh('3D_rect_mesh4.mphtxt');
tic;
%-------------------------------初始化参数
Cond = 50;
BoundTemp = 273.15;
HeatFlux = 5e3;
%-------------------------------根据插值要求求解四面体单元的几何参数
%导出四面体单元XYZ的坐标
X = Coor(:, 1);
Y = Coor(:, 2);
Z = Coor(:, 3);
TetX = X(TetElement);
TetY = Y(TetElement);
TetZ = Z(TetElement);
TetLen = length(TetElement);
%求解各单元体积Volume
Volume = zeros(TetLen, 1);
for i = 1:TetLen
    Ve = det([ones(4, 1), TetX(i, :)', TetY(i, :)', TetZ(i, :)']) / 6;
    Volume(i, 1) = Ve;
end
clear Ve;
%求解p,q,r,s
p = zeros(TetLen, 4);
q = zeros(TetLen, 4);
r = zeros(TetLen, 4);
s = zeros(TetLen, 4);
%求解p
for i = 1:TetLen
    p(i, 1) = det([TetX(i, [2,3,4])', TetY(i, [2,3,4])', TetZ(i, [2,3,4])']);
    p(i, 2) = -det([TetX(i, [1,3,4])', TetY(i, [1,3,4])', TetZ(i, [1,3,4])']);
    p(i, 3) = det([TetX(i, [1,2,4])', TetY(i, [1,2,4])', TetZ(i, [1,2,4])']);
    p(i, 4) = -det([TetX(i, [1,2,3])', TetY(i, [1,2,3])', TetZ(i, [1,2,3])']);
end
%求解q
for i = 1:TetLen
    q(i, 1) = -det([ones(3, 1), TetY(i, [2,3,4])', TetZ(i, [2,3,4])']);
    q(i, 2) = det([ones(3, 1), TetY(i, [1,3,4])', TetZ(i, [1,3,4])']);
    q(i, 3) = -det([ones(3, 1), TetY(i, [1,2,4])', TetZ(i, [1,2,4])']);
    q(i, 4) = det([ones(3, 1), TetY(i, [1,2,3])', TetZ(i, [1,2,3])']);
end
% 求解r
for i = 1:TetLen
    r(i,1) = -det([TetX(i, [2,3,4])', ones(3, 1), TetZ(i, [2,3,4])']);
    r(i,2) = det([TetX(i, [1,3,4])', ones(3, 1), TetZ(i, [1,3,4])']);
    r(i,3) = -det([TetX(i, [1,2,4])', ones(3, 1), TetZ(i, [1,2,4])']);
    r(i,4) = det([TetX(i, [1,2,3])', ones(3, 1), TetZ(i, [1,2,3])']);
end
% 求解s
for i = 1:TetLen
    s(i,1) = -det([TetX(i, [2,3,4])', TetY(i, [2,3,4])', ones(3, 1)]);
    s(i,2) = det([TetX(i, [1,3,4])', TetY(i, [1,3,4])', ones(3, 1)]);
    s(i,3) = -det([TetX(i, [1,2,4])', TetY(i, [1,2,4])', ones(3, 1)]);
    s(i,4) = det([TetX(i, [1,2,3])', TetY(i, [1,2,3])', ones(3, 1)]);
end
%-------------------------------求解边界三角形单元的几何参数
%找出热交换区域节点编号和坐标
pHeatExcgBond = find(TriEntity == 5);
TriExcgElement = TriElement(pHeatExcgBond, :);
TriExcgX = X(TriExcgElement);
TriExcgY = Y(TriExcgElement);
TriExcgZ = Z(TriExcgElement);
%求解两条向量，叉乘的模和三维三角形的面积
vec1(:, 1) = TriExcgX(:, 2) - TriExcgX(:, 1);
vec1(:, 2) = TriExcgY(:, 2) - TriExcgY(:, 1);
vec1(:, 3) = TriExcgZ(:, 2) - TriExcgZ(:, 1);
vec2(:, 1) = TriExcgX(:, 3) - TriExcgX(:, 1);
vec2(:, 2) = TriExcgY(:, 3) - TriExcgY(:, 1);
vec2(:, 3) = TriExcgZ(:, 3) - TriExcgZ(:, 1);
vec = cross(vec1,vec2);
Area = sqrt(vec(:, 1).^2 + vec(:, 2).^2 + vec(:, 3).^2) / 2;
%-------------------------------四面体单元分析和装配
S = zeros(length(Coor));
F = zeros(length(Coor), 1);
for k = 1:TetLen
    for i = 1:4
        for j = 1:4
            Se = Cond*(q(k,i)*q(k,j)+r(k,i)*r(k,j)+s(k,i)*s(k,j))/(36*Volume(k));
            S(TetElement(k,i), TetElement(k,j)) = S(TetElement(k,i), TetElement(k,j)) + Se;
        end
    end
end
%-------------------------------边界三角形单元分析和装配
for k = 1:length(TriExcgElement)
    for i = 1:3
        Fe = HeatFlux*Area(k)/3;
        F(TriExcgElement(k,i)) = F(TriExcgElement(k,i)) + Fe;
    end
end
%-------------------------------采用直接法求解
%检索温度边界条件包含的点并赋值
Temp = zeros(length(Coor), 1);
TempNodes = find(Y == 0);
FreeNodes = find(~Y == 0);
Temp(TempNodes) = BoundTemp;
%将已知部分移到等号右侧
F1 = S*Temp;
F2 = F - F1;
%求解FreeNodes
Temp(FreeNodes) = S(FreeNodes,FreeNodes)\F2(FreeNodes);
toc;
%-------------------------------后处理 先看一下y轴和z轴
Interp1 = scatteredInterpolant(Y,Z,Temp);
tx = 0:1e-2:1;
ty = 0:1e-2:0.5;
tz = 0:1e-2:1;
[qy,qz] = meshgrid(ty,tz);
Temp1 = Interp1(qy,qz);
subplot(1,2,2);
contourf(qy,qz,Temp1,20);colorbar;
hold on;
%-------------------------------读取COMSOL分网结果
[fileID, Errmessage] = fopen('solve3Drect.txt', 'r');
if fileID == -1
    disp(Errmessage);
end
for i=1:9
    fgetl(fileID);
end
comsoldata = fscanf(fileID,'%lf %lf %lf %lf\n',[4,length(Coor)]);
comsoldata = comsoldata';
fclose(fileID);
Interp2 = scatteredInterpolant(comsoldata(:,2),comsoldata(:,3),comsoldata(:,4));
Temp2 = Interp2(qy,qz);
subplot(1,2,1);
contourf(qy,qz,Temp2,20);colorbar;