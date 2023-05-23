% 该函数的作用是计算一个映射f的调和能量。
% 
% 输入参数：
% - mesh: 网格数据结构，包括edge和ew等。
% - f: 双精度数组，nv x 1的数组，表示映射f的值。
% 
% 输出参数：
% - E: 双精度数值，表示映射f的调和能量。
% 
% 函数的实现过程：
% 1. 如果网格数据结构中没有边权重ew，则调用edge_weight函数计算。
% 2. 计算f在每条边上的梯度df，即f在边的终点处减去f在边的起点处。
% 3. 计算每条边的权重与df的点积的平方，并将它们相加得到调和能量E。
% 
% 该函数的作者未知，版权归原作者所有。
function E = harmonic_energy(mesh, f)
% compute harmonic energy for a map f
edge = mesh.edge;
% if edge weight is computed already, reuse it
if ~isfield(mesh,'ew')
    ew = edge_weight(mesh);
else
    ew = mesh.ew;
end
% harmonic energy
df = f(edge(:,2),:)-f(edge(:,1),:);
E = sum(ew.*dot(df,df,2))/2;
