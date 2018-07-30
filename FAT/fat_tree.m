clc
clear
close all


K = 64; % n. pods

outfile_links = ['fat_tree-connections-K' num2str(K) '.txt'];
outfile_nodes = ['fat_tree-nodes-K' num2str(K) '.txt'];


N_core = (K/2)^2;
N_aggr = (K^2)/2;
N_edge = N_aggr;
N_serv = N_edge * K/2;
serv_nodes = (1:N_serv);
edge_nodes = N_serv + (1:N_edge);
aggr_nodes = N_serv + N_edge + (1:N_aggr);
core_nodes = N_serv + N_edge + N_aggr + (1:N_core);


N_tot = N_core + N_aggr + N_edge + N_serv;
Adj = zeros(N_tot,N_tot);

n_core_connected = zeros(1,N_aggr);

% Aggregation to edge
for kk=1:K % kk-th pod
    pod_ind = (kk-1)*K/2+1:kk*K/2;
    Adj(aggr_nodes(pod_ind),edge_nodes(pod_ind)) = ones(length(pod_ind),length(pod_ind));
end

% Edge to servers
for ee=1:N_edge
    ind_serv = (ee-1)*K/2+1:ee*K/2;
    Adj(edge_nodes(ee),serv_nodes(ind_serv)) = 1;
end

% Core to aggregation
for cc=1:N_core
    for kk=1:K % pod
        pod_ind = (kk-1)*K/2+1:kk*K/2;
        ind_aggr = find(n_core_connected(pod_ind)<K/2,1);
        Adj(core_nodes(cc),aggr_nodes(pod_ind(ind_aggr))) = 1;
        n_core_connected(pod_ind(ind_aggr)) = n_core_connected(pod_ind(ind_aggr)) + 1;
    end
end

Adj = Adj + Adj';
Adj(Adj>1) = 1;

[src,dst] = find(Adj==1);

dlmwrite(outfile_nodes,[(0:N_tot-1)', zeros(N_tot,2)],',');
dlmwrite(outfile_links,[src-1, dst-1, ones(length(src),2)],',');

