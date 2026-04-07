function EC = event_chain(e,S,num_chains,num_time_periods_per_chain)
EC = zeros(num_chains,num_time_periods_per_chain,"gpuArray");
EC(:,1) = e;
for tt=2:num_time_periods_per_chain
    
    for nn=1:num_chains
        etemp = next_event(S(:,EC(nn,tt-1)));
        EC(nn,tt) = etemp;
    end
end