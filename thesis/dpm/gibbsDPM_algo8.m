% A Gibbs sampling algorithm that uses auxiliary variables

function c_st = gibbsDPM_algo8(y, hyperG0, alpha, niter, doPlot)

    if doPlot
        figure('name','Gibbs sampling for DPM');
        colormap('default')
        cmap = colormap;
    end

    % Have n data items
    n = size(y,2);

    % Prepare cluster statistics to return
    c_st = zeros(n, niter/2);

    % Create 200 tables, with m customers each
    m = zeros(1,200);

    % Assign to each data item a table index
    c = zeros(n, 1);

    % Initialisation
    for k=1:n
        % Assign to each data item a table index uniformly up to 30
        c(k) = ceil(30*rand);
        % Add to the table counter
        m(c(k)) = m(c(k)) + 1;
        if m(c(k))>1
            % If there are already customers assigned, update the sufficient
            % statistics with the new data item
            U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
        else
            % If this is the first customer, draw sufficient statistics from
            % the hyper prior.
            U_SS(c(k)) = update_SS(y(:,k), hyperG0);
        end
    end

    niter=2;
    n=1;

    % Neal's m, the number of auxiliary variables
    n_m=2;

    % Iterations
    for i=2:niter
        % Update cluster assignments c
        for k=1:n
            % Remove data item k from the partition
            m(c(k)) = m(c(k)) - 1;

            % Assign new table index
            c(k) = sample_c(m, alpha, z, hyperG0, U_mu, U_Sigma);

            % Add data item k to table indicated by index c
            m(c(k)) = m(c(k)) + 1;


        end
    end
end

% Different from the conjugate case, we could first establish the probability
% of sampling a new or existing clusters, only only after this decide to sample
% from a new distribution. In the nonconjugate case we have no closed-form
% description of the posterior probability, hence we actually have to sample
% from our prior and only establish then the likelihood that our observation
% originated from an existing or this new cluster.
function K = sample_c(m, alpha, z, hyperG0, U_mu, U_Sigma)

    % Find first n_m empty tables
    emptyT = find(m==0, n_m);
    % This cluster does have not a number of customers, but alpha/m as weight
    m(emptyT) = alpha/n_m;
    % Get values from prior, to be used in likelihood calculation
    for i=1:length(emptyT)
        % Sample for this empty table
        R = sample_pdf(U_SS(emptyT(i)));
        switch(hyperG0.prior)
            case 'NIW'
                U_R.mu(:, emptyT(i)) = R.mu;
                U_R.Sigma(:, :, emptyT(i)) = R.Sigma;
            case 'NIG'
                U_R.mu(:, emptyT(i)) = R.mu;
                U_R.Sigma(emptyT(i)) = R.Sigma;
            case 'DPM_Seg' % Dirichlet Process Mixture of Segments
                %
                %
            otherwise
        end
    end

    % Indices of all clusters, both existing and proposed
    c = find(m~=0);

    % Calculate likelihood for every cluster
    n = m(c).*likelihoods(hyperG0, repmat(z, 1, length(c)), U_mu(:, c)', U_Sigma(c)' )';

    % Calculate b, as b=(n-1+alpha)/sum(n), of which we can forget the nominator
    const = sum(n);
    u=rand(1);
    ind = find(cumsum(n/const)>=u, 1 );
    K = c(ind);

    % Set proposed tables to 0, except for table K
    setzero=setdiff(emptyT,K);
    m(setzero)=0;
end


