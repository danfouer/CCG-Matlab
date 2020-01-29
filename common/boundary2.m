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
