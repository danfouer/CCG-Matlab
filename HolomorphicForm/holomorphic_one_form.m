% 这是一个计算给定曲面上特征1-形式的全体全纯1-形式的MATLAB函数。下面我们逐步解释代码：
% 
% 1. 首先定义了一个函数 `holomorphic_one_form`，它有两个输入参数 `mesh` 和 `cf`，
% 分别表示曲面的连接性和顶点坐标、特征1-形式。输出参数 `holomorphic_form_basis` 
% 是一个与 `cf` 大小相同的向量，表示全体全纯1-形式。
% 
% 2. 如果没有给定特征1-形式 `cf`，则通过 `characteristic_one_form` 函数计算特征1-形式。
% 
% 3. 通过 `harmonic_form` 函数计算特征1-形式的调和1-形式。
% 
% 4. 通过 `conjugate_harmonic_one_form` 函数计算调和1-形式的共轭调和1-形式。
% 
% 5. 将调和1-形式和共轭调和1-形式组合成复数形式的1-形式，得到全纯1-形式。
% 
% 6. 将全纯1-形式存储在 `holomorphic_form_basis` 向量中，并返回。
function holomorphic_form_basis = holomorphic_one_form(mesh, cf)

if ~exist('cf','var') || isempty(cf)
    cf = characteristic_one_form(mesh);
end

hf = cell(size(cf));
for i = 1:length(cf)
    hf{i} = harmonic_form(mesh, cf{i});
end

dual_hf = cell(size(cf));
for i = 1:length(dual_hf)
    dual_hf{i} = conjugate_harmonic_one_form(mesh, hf{i}, hf);
end

holomorphic_form_basis = cell(size(dual_hf));
for i = 1:length(holomorphic_form_basis)
    holomorphic_form_basis{i} = complex(hf{i}, dual_hf{i});
end
