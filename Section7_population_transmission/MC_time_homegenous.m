%Markov chain simulation
%%
% Time homogeneous

% Build matrix
% Order N, I, V, I', V'
% Define the transition matrix using S = transition_matrix(NtoI, NtoV,IprimetoI,IprimetoV,VprimetoI,VprimetoV,IprimetoN, VprimetoN)
S_1 = transition_matrix(0.1,0.01,0.08,0.005,0.05,0.006,0.01,0.01); %high inf, low vac
S_2 = transition_matrix(0.1,0.01,0.08,0.005,0.05,0.006,0.1,0.1); %high inf, low vac, high return to Naive
S_3 = transition_matrix(0.1,0.1,0.08,0.05,0.05,0.06,0.01,0.01); %high inf, high vac

S = S_1;


%%
% Each column of SP: state probabilities after j time steps if you start at v. 
% Plot SP' to see how each state probability changes over 50 iterations if
% you begin at v.
% Can help you identify limiting distribution for the Markov chain.

v = [1;0;0;0;0];
n = 50;
SP = state_probs_n(S,v,n);
%%

% Sample event chains

e = 1; % every event chain beginning in Naive state
n_chains = 1000; %Number of chains, start with a smaller number to test code, paper uses 1000.
n_time = 50; %Number of time steps per chain, start with a smaller number to test code, paper uses 50.
% timesteps per chain, each timestep corresponding to two weeks or so, 50 thus means 2 years

tic             
EC = event_chain(1,S,n_chains,n_time); %EC: each row is a chain.
toc

%plot(EC','LineWidth',2);
per_chain_naive = sum(EC'==1);
per_chain_infections = sum(EC'==2);
per_chain_vaccinations = sum(EC'==3);
per_chain_prev_infections = sum(EC'==4);
per_chain_prev_vaccinations = sum(EC'==5);

%%
figure()
tiledlayout(3,4);
nexttile([1,2])
histogram(per_chain_infections,'DisplayName','Infections per individual','FaceColor','red')
xlabel('Count per individual')
ylabel('Number of individuals')
legend
nexttile([1,2])
histogram(per_chain_vaccinations,'DisplayName','Vaccinations per individual','FaceColor','blue')
xlabel('Count per individual')
ylabel('Number of individuals')
legend
nexttile([1,2])
histogram(per_chain_prev_infections,'DisplayName','Previously Infected status per individual','FaceColor','red','Facealpha',0.3)
xlabel('Count per individual')
ylabel('Number of individuals')
legend
nexttile([1,2])
histogram(per_chain_prev_vaccinations,'DisplayName','Previously Vaccinated status per individual','FaceColor','blue','Facealpha',0.3)
xlabel('Count per individual')
ylabel('Number of individuals')
legend
nexttile(10,[1,2])
histogram(per_chain_naive,'DisplayName','Naive status per individual','FaceColor',[0.996 1 0.85])%'green')
xlabel('Count per individual')
ylabel('Number of individuals')
legend
hold off
fixplots(20)



%%
NtoItransitions = zeros(n_chains,1);
IprimetoItransitions = zeros(n_chains,1);
VprimetoItransitions = zeros(n_chains,1);

IprimetoVtransitions = zeros(n_chains,1);
VprimetoVtransitions = zeros(n_chains,1);


for ii=1:n_chains
    tempvec = EC(ii,:);
    indexNtoI = strfind(gather(tempvec), [1,2]);
    NtoItransitions(ii) = numel(indexNtoI);
    indexIprimetoI = strfind(gather(tempvec), [4,2]);
    IprimetoItransitions(ii) = numel(indexIprimetoI);
    indexVprimetoI = strfind(gather(tempvec), [5,2]);
    VprimetoItransitions(ii) = numel(indexVprimetoI);
    indexIprimetoV = strfind(gather(tempvec), [4,3]);
    IprimetoVtransitions(ii) = numel(indexIprimetoV);
    indexVprimetoV = strfind(gather(tempvec), [5,3]);
    VprimetoVtransitions(ii) = numel(indexVprimetoV);

end

%ismember(EC,[1,2])


%%
 UniqueSeq = unique(EC,'rows');
 [tempr,tempc] = size(UniqueSeq);
 Rcount = zeros(tempr,1);
 for ii=1:tempr
     R = EC(ii,:);
     Rcount(ii) = nnz(all(EC==R,2));
 end
% Use plot(Rcount) to identify frequent/representative sequences
% For example, %S1: 281, 412, 931 from UniqueSeq were highest for me.
% The event chain selection is stochastic; every run will get slightly
% different results.
%%
%Uncomment the following code by replacing appropriately identified indices
%from the previous few lines instead of the 281, 412, 931 in the following
%exmple.

%tempcdata = [UniqueSeq(281,:);UniqueSeq(412,:);UniqueSeq(931,:)]; %S1

%Uncomment the following code to plot the heatmap corresponding to
%representative individual sequences.

% tempcdata=gather(tempcdata);
% figure()
% h= heatmap(tempcdata);
% mycolors = [0.996 1 0.85; 1 0 0; 0 0 1; 1 0.5 0.5; 0.5 0.5 1];
% colormap(mycolors);
% h.Colormap

%%
% median(per_chain_naive)
% median(per_chain_infections)
% median(per_chain_prev_infections)
% median(per_chain_vaccinations)
% median(per_chain_prev_vaccinations)
