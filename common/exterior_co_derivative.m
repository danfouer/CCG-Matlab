% 这是一个用于计算外共边缘算子的 MATLAB 函数。该函数的输入参数包括一个网格 mesh、一个差分形式 w 
% 和一个阶数 order。其中，网格 mesh 包括一个边缘数组 edge 和一个顶点数 nv，差分形式 w 的大小应该
% 与边缘数相同或者与面数相同。
% 
% 该函数的输出参数为一个 delta_w，表示差分形式 w 的外共边缘算子。
% 
% 该函数的实现过程如下：
% 
% 1. 从网格 mesh 中获取边缘数组 edge 和顶点数 nv。
% 
% 2. 如果网格 mesh 中没有边权重数组 ew，则调用 edge_weight 函数计算边权重数组 ew。
% 
% 3. 如果阶数 order 为 1，则检查差分形式 w 的大小是否与边缘数相同，如果不同则抛出错误。然后，
% 调用 accumarray 函数计算每个顶点的外共边缘算子值，得到一个 nv x 1 的数组 delta_w。
% 
% 4. 如果阶数 order 为 2，则检查差分形式 w 的大小是否与面数相同，如果不同则抛出错误。此时，
% 函数没有计算外共边缘算子的代码，需要根据具体需求进行补充。
function delta_w = exterior_co_derivative(mesh, w, order)
% exterior co-derivative of differential form w
edge = mesh.edge;
nv = mesh.nv;
if ~isfield(mesh,'ew')
    ew = edge_weight(mesh);
else
    ew = mesh.ew;
end
if order == 1
    if size(w) ~= mesh.ne
        error('differential form w has incorrect shape');
    end
    delta_w = accumarray(edge(:,1),ew.*w,[nv,1]) + accumarray(edge(:,2),-ew.*w,[nv,1]);
end
if order == 2
    if size(w) ~= mesh.nf
        error('differential form w has incorrect shape');
    end
end
    