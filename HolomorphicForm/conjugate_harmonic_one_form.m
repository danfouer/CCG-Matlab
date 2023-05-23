% 这是一个计算给定1-形式的共轭调和1-形式的MATLAB函数。下面我们逐步解释代码：
% 
% 1. 首先定义了一个函数 `conjugate_harmonic_one_form`，它有三个输入参数 `mesh`、`w` 和 `hf`，
% 分别表示曲面的连接性和顶点坐标、给定的1-形式和共轭调和1-形式。输出参数 `cw` 是一个与 `w` 大小相同的向量，表示给定1-形式的共轭调和1-形式。
% 
% 2. 如果没有给定共轭调和1-形式 `hf`，则通过 `characteristic_one_form` 函数计算特征1-形式，
% 并通过 `harmonic_form` 函数计算特征1-形式的调和1-形式。
% 
% 3. 将调和1-形式转换为向量表示，并计算给定1-形式的叉积，得到对偶1-形式。
% 
% 4. 计算每个面的法向量和面积，并计算调和1-形式和对偶1-形式的叉积与法向量的点积的加权和，得到线性方程组的右侧向量。
% 
% 5. 构建线性方程组的系数矩阵，其中每个元素是两个调和1-形式的楔积在曲面上的积分。
% 
% 6. 通过求解线性方程组 $Ax = s$，得到系数向量 $\mu$。
% 
% 7. 将调和1-形式和系数向量 $\mu$ 相乘，并将结果相加，得到共轭调和1-形式 `cw`。
function cw = conjugate_harmonic_one_form(mesh, w, hf)
% compute conjugate harmonic 1-form of w

if ~exist('hf','var') || isempty(hf)
    if ~exist('cf', 'var')
        cf = characteristic_one_form(mesh);
    end
    hf = cell(size(cf));
    for i = 1:length(hf)
        hf{i} = harmonic_form(mesh, cf{i});
    end
end

hf_vector = cell(size(hf));
for i = 1:length(hf)
    hf_vector{i} = vector_representation(mesh, hf{i});
end

fnormal = face_normal(mesh.face, mesh.vert);
farea = face_area(mesh.face, mesh.vert);

w_vector = vector_representation(mesh, w);
dual_w = cross(fnormal, w_vector);

% build linear system
s = ones(length(hf), 1);
for i = 1:length(hf)  
    s(i) = sum(dot(cross(hf_vector{i}, dual_w), fnormal, 2).*farea);
end

A = zeros(length(hf), length(hf));
for i = 1 : length(hf)
    for j = 1 : length(hf)
        wp = wedge_product(mesh, hf{i}, hf{j});
        %integrate this 2-form
        A(i,j) = sum(wp.*farea);
    end
end
mu = A\s;

cw = zeros(size(w));
for i = 1 : length(hf)
    cw = cw + mu(i) * hf{i};
end

end