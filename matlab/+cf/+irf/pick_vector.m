function x = pick_vector(post, k, names)
for nm = names
    if isfield(post, nm)
        v = post.(nm);
        if ndims(v)==2
            % could be Nd x dim  OR dim x Nd
            if size(v,1)==k || size(v,1)>k
                x = v(k,:)';
                return
            elseif size(v,2)>=k
                x = v(:,k);
                return
            end
        elseif ndims(v)==3
            % dim x 1 x Nd or dim x Nd x something...
            if size(v,3)>=k
                x = v(:,:,k);
                x = x(:);
                return
            end
        end
    end
end
error("Cannot find vector field among: %s", strjoin(string(names),", "));
end
