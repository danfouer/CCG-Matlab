% 这是一个用于计算网格边界的 MATLAB 函数。该函数接受一个网格作为输入，返回一个单元格数组，
% 其中每个单元格包含一个边界的顶点索引。
% 
% 函数实现过程：
% 
% - 首先，通过网格的 halfedge 和 edge 属性计算出边界的边 be，并使用这些边构建一个邻接矩阵 Ge。
% - 然后，通过 halfedge 和邻接矩阵 Ge 找到所有的边界半边 he，并使用这些半边构建一个有向邻接矩阵 G。
% - 接下来，找到所有的边界顶点，并将它们分组成单独的边界。
% - 最后，将每个边界的顶点索引存储在单元格数组 bds 中，并返回 bds。
% 
% 该函数的作者没有给出详细的注释，但是根据代码可以大致理解其实现过程。
function bds = boundary2(mesh)
nv = mesh.nv;
he = mesh.halfedge;
edge = mesh.edge;
eif = mesh.eif;
ind = eif(:,1)<0|eif(:,2)<0;
be = edge(ind,:);
Ge = sparse(be(:,1),be(:,2),ones(size(be,1),1),nv,nv);
Ge = Ge+Ge';
ind = full(Ge(sub2ind([nv,nv],he(:,1),he(:,2)))==1);
% find boundary halfedge
he = he(ind,:);
G = sparse(he(:,1),he(:,2),ones(size(he,1),1),nv,nv);
bds = {};
ind = false(nv,1);
ind(he) = true;
k = 1;
while true
    % find a boundary vert
    b = find(ind,1,'first');
    if isempty(b)
        break;
    end
    bd = b;
    ind(b) = false;
    while true
        bs = find(G(b,:));
        if isempty(bs)
            error(['unable to find halfedge with source vert ' num2str(b)]);
        end
        
        if length(bs)>1
            ib = ind(bs);
            i = find(ib,1);
            if isempty(i)
                break
            end
            bs = bs(i);
        end
        if ~ind(bs)
            break
        end

        ind(bs) = false;
        b = bs;
        bd = [bd;b];
    end
    bds{k} = bd;
    k = k+1;
end
