% 这是一个用于计算外导数算子的 MATLAB 函数。该函数的输入参数包括一个网格 mesh、一个差分形式 w 和一个阶数 order。
% 其中，网格 mesh 包括一个边缘数组 edge、一个面数组 face 和一个顶点数 nv，差分形式 w 的大小应该与顶点数相同
% 或者与边缘数相同。
% 
% 该函数的输出参数为一个 dw，表示差分形式 w 的外导数算子。
% 
% 该函数的实现过程如下：
% 
% 1. 从网格 mesh 中获取边缘数组 edge。
% 
% 2. 如果阶数 order 为 0，则检查差分形式 w 的大小是否与顶点数相同，如果不同则抛出错误。然后，调用差分形式 w 
% 的起点和终点计算每条边缘的外导数值，得到一个 ne x 1 的数组 dw。
% 
% 3. 如果阶数 order 为 1，则检查差分形式 w 的大小是否与边缘数相同，如果不同则抛出错误。此时，函数需要根据面
% 和边缘的关系计算每个面的外导数值。首先，将差分形式 w 转换为一个 nv x nv 的稀疏矩阵 ws，其中 ws(i,j) 表示
% 从顶点 i 到顶点 j 的边缘上的差分形式值。然后，计算 ws 的上三角矩阵和下三角矩阵之差，得到一个 nv x nv 的矩阵，
% 表示每条边缘的外导数值。最后，根据面和边缘的关系，将每个面的外导数值相加，得到一个 nf x 1 的数组 dw。
function dw = exterior_derivative(mesh, w, order)
% exterior derivative of differential form w
edge = mesh.edge;
if order == 0
    if size(w,1) ~= mesh.nv
        error('differential form w has incorrect shape');
    end
    dw = w(edge(:,2)) - w(edge(:,1));
end
if order == 1
    if size(w,1) ~= mesh.ne
        error('differential form w has incorrect shape');
    end
    face = mesh.face;
    nv = mesh.nv;
    ws = sparse(edge(:,1),edge(:,2),w,nv,nv);
    ws = ws - ws';
    dw = zeros(mesh.nf,1);
    dw = dw + ws(face(:,1)+(face(:,2)-1)*nv);
    dw = dw + ws(face(:,2)+(face(:,3)-1)*nv);
    dw = dw + ws(face(:,3)+(face(:,1)-1)*nv);
end
