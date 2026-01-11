function Nd = posterior_num_draws(post)
% Number of posterior draws for THIS repo format (post = res.draws).

if ~isstruct(post) || ~isfield(post,'Pstar')
    error("posterior_num_draws: expected post.Pstar (post should be res.draws).");
end

Nd = numel(post.Pstar);

if Nd < 1
    error("posterior_num_draws: post.Pstar is empty.");
end
end

