% 该函数用于计算两个一形式的楔积。
% 
% 输入参数：
% - mesh: 包含网格信息的结构体。
% - w1: double类型，n x 3，表示第一个一形式。
% - w2: double类型，n x 3，表示第二个一形式。
% 
% 输出参数：
% - wp: double类型，nf x 1，表示两个一形式的楔积。
% 
% 函数的实现过程：
% 1. 使用vector_representation函数将w1和w2转换为向量表示形式。
% 2. 使用face_normal函数计算每个面的法向量。
% 3. 使用cross函数计算w1_vector和w2_vector的叉积。
% 4. 使用dot函数计算叉积结果和每个面法向量的点积，得到楔积。
function wp = wedge_product(mesh, w1, w2)
% compute wedge product of two one-forms
w1_vector = vector_representation(mesh, w1);
w2_vector = vector_representation(mesh, w2);
f_normal = face_normal(mesh.face, mesh.vert);
wp = dot(cross(w1_vector, w2_vector), f_normal, 2);
end