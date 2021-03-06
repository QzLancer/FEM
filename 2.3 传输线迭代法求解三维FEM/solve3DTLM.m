%通过传输线法求解三维稳态热场
%包含第一类边界条件和第二类边界条件
%材料的热导率为50W/(m.k)，四面热绝缘，一面的热通量为5e3W/m2，一面的温度为273.15K
%参考《The finite element method in engineering》 Chapter16.4 P543-P545
%测试TriEntity与边界面的对应关系， TriEntity = 5时对应热通量面，TriEntity = 2对应恒温面
%其余面为热绝缘面
%（在读取分网程序中将Entity+1可以获得与COMSOL一致的编号）
%存在负电阻，所以需要将负电阻等效成电流源，且不在负电阻的两端添加传输线
%By QzLancer
%2019/3/15
%-------------------------------读取分网文件
clear all;
[Coor,VtxElement,VtxEntity,TriElement,TriEntity,TetElement,TetEntity] = read3Dmesh('3D_rect_mesh.mphtxt');
%-------------------------------初始化参数
Cond = 50;
BoundTemp = 273.15;
HeatFlux = 5e3;
Y0 = 1; 
e = 0.001;
tic
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
Se = zeros(4,4,TetLen);
Ym = zeros(length(Coor),length(Coor));
F = zeros(length(Coor), 1);
for k = 1:TetLen
    for i = 1:4
        for j = 1:4
            Se(i,j,k) = Cond*(q(k,i)*q(k,j)+r(k,i)*r(k,j)+s(k,i)*s(k,j))/(36*Volume(k));
            if i == j
                Ym(TetElement(k,i),TetElement(k,j)) = Ym(TetElement(k,i),TetElement(k,j))+3*Y0;
            else
                if Se(i,j,k) < 0
                    Ym(TetElement(k,i),TetElement(k,j)) = Ym(TetElement(k,i),TetElement(k,j))-Y0;
                else
                    Ym(TetElement(k,i),TetElement(k,i)) = Ym(TetElement(k,i),TetElement(k,i))-Y0;
                end
            end
        end
    end
end
%-------------------------------边界三角形单元分析和装配
for k = 1:length(TriExcgElement)
    for i = 1:3
        Fe = HeatFlux*Area(k)/3;
        F(TriExcgElement(k,i)) = F(TriExcgElement(k,i))+Fe;
    end
end
%-------------------------------采用TLM求解
%检索温度边界条件包含的点并赋值
TempNodes = find(Y == 0);
FreeNodes = find(~Y == 0);
%入射和反射过程
Va0 = ones(length(Coor),1);
Va1 = zeros(length(Coor),1);
Vr = zeros(4,4,length(TetElement));
Vc = zeros(4,4,length(TetElement));                              
Vi = zeros(4,4,length(TetElement));
I = zeros(length(Coor),1);
m = 1;
while norm(Va1-Va0)/norm(Va0) >= e
    Va0 = Va1;
    Va1 = zeros(length(Coor),1);
    Va1(TempNodes) = BoundTemp;
    F1 = Ym*Va1;
    F2 = F+I-F1;
    Va1(FreeNodes) = Ym(FreeNodes,FreeNodes)\F2(FreeNodes);
    I = zeros(length(Coor),1);
    for k = 1:TetLen
        for i = 1:4
            for j = 1:4
                if Se(i,j,k) < 0
                    Vr(i,j,k) = Va1(TetElement(k,i))-Va1(TetElement(k,j))-Vi(i,j,k);
                    Vc(i,j,k) = 2*Vr(i,j,k)*Y0/(Y0-Se(i,j,k));
                    Vi(i,j,k) = Vc(i,j,k) - Vr(i,j,k);
                    I(TetElement(k,i)) = I(TetElement(k,i))+2*Y0*Vi(i,j,k);         
                else
                    I(TetElement(k,i)) = I(TetElement(k,i))+Se(i,j,k)*(Va1(TetElement(k,i))-Va1(TetElement(k,j)));
                end
            end
        end
    end
    m = m+1;
end
toc
%-------------------------------后处理 先看一下y轴和z轴
Interp1 = scatteredInterpolant(Y,Z,Va1);
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