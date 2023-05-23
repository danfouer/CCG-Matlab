% 这段 MATLAB 代码实现了计算三角网格上的微分 1-形式的积分。具体实现过程如下：
% 
% 首先，从输入的 mesh 结构中获取 edge 和 nv，分别表示三角网格的边数和顶点数。同时，获取 dw 向量，
% 它是一个大小为 ne x 1 的向量，表示每条边的权重。
% 
% 然后，根据 edge 和 dw 向量构造一个大小为 nv x nv 的稀疏矩阵 es，其中 es(i,j) 表示从顶点 i 到顶点 j 的边的权重。
% 注意，由于微分 1-形式是一个复值函数，因此这里使用了 conj 函数对权重取共轭。
% 
% 接着，使用 vert_vert_ring 函数获取每个顶点的一环邻域。初始化一个大小为 nv x 1 的向量 f，其中 f(1) = 0，
% 表示从第一个顶点出发的积分值为 0。同时，初始化一个大小为 nv x 1 的逻辑向量 ind，其中 ind(i) 表示顶点 i 
% 是否已经被访问过。将 ind(1) 设为 true，表示第一个顶点已经被访问过。初始化一个空队列 qe，将第一个顶点加入队列。
% 
% 接下来，从队列中取出一个顶点 i，遍历它的一环邻域 vri。对于每个未被访问过的邻居 j，计算从 i 到 j 的边的权重
% es(i,j)，并将 f(j) 设为 f(i) 加上这个权重。将 ind(j) 设为 true，表示 j 已经被访问过。将 j 加入队列 qe，
% 以便后续访问。
% 
% 重复上述过程，直到队列 qe 为空。此时，f 向量中存储的就是从第一个顶点出发的微分 1-形式的积分值。
% 将 f 向量作为函数的输出。
function f = integration(mesh, dw)
% integration of differential one form
edge = mesh.edge;
nv = mesh.nv;
es = sparse(edge(:,1),edge(:,2),dw,nv,nv);
es = es - conj(es');
vvr = vert_vert_ring(mesh);
f = zeros(mesh.nv,1);
f(1) = 0;
ind = false(nv,1);
ind(1) = true;
qe = 1;
while ~isempty(qe)
    i = qe(end);
    qe(end) = [];
    vri = vvr{i};
    for j = vri
        if ind(j)
            continue
        end
        f(j) = f(i) + es(i,j);
        ind(j) = true;
        qe = [qe,j];
    end
end
