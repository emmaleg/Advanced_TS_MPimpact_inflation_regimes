function x = pick_scalar(post, k, names)
for nm = names
    if isfield(post, nm)
        v = post.(nm);
        if isvector(v)
            x = v(k);
            return
        elseif ndims(v)==2 && size(v,1)>=k
            x = v(k,1);
            return
        end
    end
end
error("Cannot find scalar field among: %s", strjoin(string(names),", "));
end
