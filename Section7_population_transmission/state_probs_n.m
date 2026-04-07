function SP = state_probs_n(S,v,n)
    SP = zeros(5,n,"gpuArray");
    temp = v;
    for jj=1:n
        SP(:,jj) = S*temp;
        temp = SP(:,jj);
    end
end

