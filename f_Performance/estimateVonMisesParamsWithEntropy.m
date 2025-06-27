function [mu, kappa, entropy] = estimateVonMisesParamsWithEntropy(data, num_bins)
    % estimateVonMisesParamsWithEntropy estimates the parameters of a Von Mises distribution
    % and computes the entropy of the data.
    %
    % Inputs:
    %   data     - Input angles in radians (vector)
    %   num_bins - Number of bins for entropy calculation (default: 50)
    %
    % Outputs:
    %   mu       - Estimated mean direction
    %   kappa    - Estimated concentration parameter
    %   entropy  - Computed entropy of the data

    % Handle default num_bins
    if nargin < 2
        num_bins = 50;
    end

    % Number of data points
    n = length(data);
    
    % Compute sums of sin and cos
    S = sum(sin(data));
    C = sum(cos(data));
    
    % Compute mean direction (mu)
    mu = atan2(S, C);
    
    % Compute mean resultant length (R)
    R = sqrt(S^2 + C^2) / n;
    
    % Estimate concentration parameter (kappa)
    if R < 0.53
        kappa = 2 * R + R^3 + 5 * R^5 / 6;
    elseif R < 0.85
        kappa = -0.4 + 1.39 * R + 0.43 / (1 - R);
    else
        kappa = 1 / (R^3 - 4 * R^2 + 3 * R);
    end

    % Compute entropy
    % Bin the data
    bin_edges = linspace(-pi, pi, num_bins + 1);
    counts = histcounts(data, bin_edges);
    
    % Normalize counts to probabilities
    probabilities = counts / sum(counts);
    
    % Compute entropy
    entropy = -sum(probabilities .* log2(probabilities + eps), 'omitnan');
end
