% 该函数的作用是计算一个闭合1-形式w的向量值2-形式。其中，闭合1-形式是指一个从每个边的起点到终点的赋值函数，
% 它满足对于每个面片，该面片上的边的起点和终点的赋值之和为0。
% 
% 函数的输入参数有两个，分别是网格(mesh)和闭合1-形式w。其中，mesh是一个网格结构体，包含面片(face)、
% 顶点(vertex)和边(edge)等信息；w是一个nv x 1的列向量，表示闭合1-形式在每个顶点处的取值。
% 
% 函数的输出参数是一个nf x 3的双精度数组vf，表示向量值2-形式在每个面片上的取值。其中，nf是面片的数量。
% 
% 函数的具体实现过程如下：
% 
% 1. 通过mesh中的面片(face)和顶点(vertex)计算每个面片的法向量normal。
% 
% 2. 根据闭合1-形式w构建一个nv x nv的稀疏矩阵es，其中es(i,j)表示从i到j的边上的w值。
% 
% 3. 将es转换为一个对称矩阵，即es = es - conj(es')。
% 
% 4. 对于每个面片，计算其三条边上的w值，并将其分别乘以该边所对应的向量和法向量的叉积，得到该面片上的向量值2-形式。
% 
% 5. 将所有面片上的向量值2-形式合并到一起，得到vf数组。
% 
% 6. 返回vf数组。
function vf = vector_representation(mesh, w)
% vector valued 2-form of closed 1-form w

edge = mesh.edge;
vert = mesh.vert;
face = mesh.face;
normal = face_normal(face, vert);

nv = mesh.nv;
es = sparse(edge(:,1),edge(:,2),w,nv,nv);
es = es - conj(es');

vf = - full(diag(es(face(:, 1), face(:, 2)))) .* cross(vert(face(:, 3), :), normal) ...
    - full(diag(es(face(:, 3), face(:, 1)))) .* cross(vert(face(:, 2), :), normal) ...
    - full(diag(es(face(:, 2), face(:, 3)))) .* cross(vert(face(:, 1), :), normal);

end
