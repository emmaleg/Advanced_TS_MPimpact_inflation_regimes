function M = pick_matrix(post, k, names, n, K)
for nm = names
    if isfield(post, nm)
        v = post.(nm);
        if ndims(v)==3 && size(v,3)>=k
            M = v(:,:,k);
            return
        elseif ndims(v)==2
            % maybe stored as vec per draw: Nd x (n*K)
            if size(v,1)>=k && size(v,2)==n*K
                M = reshape(v(k,:).', [K,n]).';
                return
            end
        end
    end
end
error("Cannot find matrix field among: %s", strjoin(string(names),", "));
end
