function x = interp_ctr(tGrid, xGrid, t)
    % Check for the validity of input tGrid
    validateGrid(tGrid);

    % Figure out sizes
    n = floor((length(tGrid) - 1) / 2);
    m = size(xGrid, 1);
    k = length(t);
    x = zeros(m, k);

    % Find the segment for each value of t
    bin = findSegmentIndices(t, tGrid);

    % Loop over each quadratic segment
    for i = 1:n
        idx = bin == (i + 1);
        if any(idx)
            gridIdx = 2 * (i - 1) + [1, 2, 3];
            x(:, idx) = interpQuadSegment(tGrid(gridIdx), xGrid(:, gridIdx), t(idx));
        end
    end

    % Check for any points that are exactly on the upper grid point:
    x(:, t == tGrid(end)) = xGrid(:, end);
end

function validateGrid(tGrid)
    nGrid = length(tGrid);
    if mod(nGrid - 1, 2) ~= 0 || nGrid < 3
        error('The number of grid-points must be odd and at least 3');
    end
end

function bin = findSegmentIndices(t, tGrid)
    edges = [-inf, tGrid(1:2:end), inf];
    [~, bin] = histc(t, edges);
end

function x = interpQuadSegment(tGrid, xGrid, t)
    % Rescale the query points to be on the domain [-1,1]
    t = 2 * (t - tGrid(1)) / (tGrid(3) - tGrid(1)) - 1;

    % Compute the coefficients
    a = 0.5 * (xGrid(:, 3) + xGrid(:, 1)) - xGrid(:, 2);
    b = 0.5 * (xGrid(:, 3) - xGrid(:, 1));
    c = xGrid(:, 2);

    % Evaluate the polynomial for each dimension of the function
    p = length(t);
    m = size(xGrid, 1);
    x = zeros(m, p);
    tt = t.^2;
    for i = 1:m
        x(i, :) = a(i) * tt + b(i) * t + c(i);
    end
end