% A Gibbs sampling algorithm that uses auxiliary variables
%
% -- Function: c_st = gibbsDPM_algo8(y, hyperG0, alpha, niter, doPlot)
%     Apply Gibbs sampling with m auxiliary variables.
%
%     The data is provided in y. In case of dependent plus independent
%     variables as in a regression problem, y contains both.
%
%     The hyperparameters are provided in hyperG0.
%
%     The Dirichlet Process Mixture parameter is given through alpha.
%
%     The number of iterations is defined in niter.
%
%     There are three plotting options though doPlot. If it is 0 nothing will
%     be plotted, if it is 1 every step will be plotted, if it is 2 the plot
%     will be updated every iteration for all data points at once.
%
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

        switch (hyperG0.prior)
        case { 'NIW', 'NIG' }
            if m(c(k))>1
                % If there are already customers assigned, update the sufficient
                % statistics with the new data item
                U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
            else
                % If this is the first customer, draw sufficient statistics from
                % the hyper prior.
                U_SS(c(k)) = update_SS(y(:,k), hyperG0);
            end
        case 'DPM_Seg'
             % Nothing to do...
        otherwise
        end
    end

    % Sample parameters for the tables (unique indices in customer allocation
    % array c)
    ind = unique(c);
    for j=1:length(ind)
        switch (hyperG0.prior)
        case 'NIW'
            R = sample_pdf(hyperG0);
            U_R.mu(:, ind(j)) = R.mu;
            U_R.Sigma(:, :, ind(j)) = R.Sigma;
        case 'NIG'
            R = sample_pdf(hyperG0);
            U_R.mu(:, ind(j)) = R.mu;
            U_R.Sigma(ind(j)) = R.Sigma;
        case 'DPM_Seg'
            R = sample_pdf(hyperG0);
            U_R.mu(:, ind(j)) = R.mu;
            U_R.Sigma(ind(j)) = R.Sigma;
            U_R.a(emptyT(i)) = R.a;
            U_R.b(emptyT(i)) = R.b;
        end
    end

    % Number of sampling steps is niter
    for i=2:niter
        % Update cluster assignments c
        for k=1:n
            % Remove data item k from the partition
            m(c(k)) = m(c(k)) - 1;

            % Assign new table index
            [c(k) update, R] = sample_c(m, alpha, y(:,k), hyperG0, U_R);
            if (update)
                U_R.mu(:, c(k)) = R.mu;
                U_R.Sigma(c(k)) = R.Sigma;
                switch (hyperG0.prior)
                case 'DPM_Seg'
                    U_R.a(c(k)) = R.a;
                    U_R.b(c(k)) = R.b;
                end
            end

            % Add data item k to table indicated by index c
            m(c(k)) = m(c(k)) + 1;

            if doPlot==1
                some_plot(y, hyperG0, U_R.mu, m, c, k, i, cmap)
            end
        end
        if i>niter/2
            c_st(:, i-niter/2) = c;
        end

        if doPlot==2
            some_plot(y, hyperG0, U_R.mu, m, c, k, i, cmap)
        end

        print_clusters=true;
        if (print_clusters)
            fprintf('Iteration %d/%d\n', i, niter);
            fprintf('%d clusters\n\n', length(unique(c)));
        end
    end
end

% Different from the conjugate case! In the conjugate case we could first
% establish the probability of sampling a new or existing cluster. Only after
% that we decided to sample from the new distribution.
% Now, in the nonconjugate case we have no closed-form description of the
% posterior probability, hence we actually have to sample from our prior.
% Then we treat the proposed cluster just as an existing one and calculate the
% likelihood in one sweep.
function [K, update, R] = sample_c(m, alpha, z, hyperG0, U_R)

    % Neal's m, the number of auxiliary variables
    n_m=3;

    % Find first n_m empty tables
    emptyT = find(m==0, n_m);
    % This cluster does have not a number of customers, but alpha/m as weight
    m(emptyT) = alpha/n_m;
    % Get values from prior, to be used in likelihood calculation
    for i=1:length(emptyT)
        % Sample for this empty table
        switch(hyperG0.prior)
            case 'NIW'
                % Sample from the prior, not taken into account data items
                R = sample_pdf(hyperG0);
                U_R.mu(:, emptyT(i)) = R.mu;
                U_R.Sigma(:, :, emptyT(i)) = R.Sigma;
            case 'NIG'
                R = sample_pdf(hyperG0);
                U_R.mu(:, emptyT(i)) = R.mu;
                U_R.Sigma(emptyT(i)) = R.Sigma;
            case 'DPM_Seg' % Dirichlet Process Mixture of Segments
                R = sample_pdf(hyperG0);
                U_R.mu(:, emptyT(i)) = R.mu;
                U_R.Sigma(emptyT(i)) = R.Sigma;
                U_R.a(emptyT(i)) = R.a; % endpoint
                U_R.b(emptyT(i)) = R.b; % endpoint
            otherwise
        end
    end

    % Indices of all clusters, both existing and proposed
    c = find(m~=0);

    % Calculate likelihood for every cluster
    n = m(c).*likelihoods(hyperG0, repmat(z, 1, length(c)), U_R.mu(:, c)', U_R.Sigma(c)' )';

    % Calculate b, as b=(N-1+alpha)/sum(n), of which the nominator (N-1+alpha)
    % gets cancelled out again, so we only require 1/const = 1/sum(n)
    const = sum(n);

    % Sample cluster in n according to their weight n(c)
    u=rand(1);
    ind = find(cumsum(n/const)>=u, 1 );
    K = c(ind);

    % Set proposed tables to 0, except for table K
    setzero=setdiff(emptyT,K);
    m(setzero)=0;

    % The update flag is used to set U_mu/U_Sigma outside this function,
    % because octave/matlab doesn't support pass-by-reference
    update=true;
    if (length(setzero) == length(emptyT))
        update=false;
    end
end

% Plot mean values
function some_plot(z, hyperG0, U_mu, m, c, k, i, cmap)
    if strcmp(hyperG0.prior, 'NIG')
        z=z(2:end,:);
    end
    ind=find(m);
    hold off;
    for j=1:length(ind)
        plot(z(1,c==ind(j)),z(2,c==ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 15);
        hold on
        plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 30);
        plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'ok', 'linewidth', 2, 'markersize', 10);
    end
    plot(z(1,k),z(2,k),'or', 'linewidth', 3)
    title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
    xlabel('X');
    ylabel('Y');
    xlim([-1 1]*20);
    ylim([-1 1]*20);
    pause(.01)
end
