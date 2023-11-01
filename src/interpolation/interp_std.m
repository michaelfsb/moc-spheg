function x = interp_std(tGrid, xGrid, fGrid, t)
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
            [h, xLow, fLow, fMid, fUpp, alpha] = getSegmentInfo(i, tGrid, xGrid, fGrid, t, idx);
            x(:, idx) = interpCubicSegment(h, xLow, fLow, fMid, fUpp, alpha);
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

function [h, xLow, fLow, fMid, fUpp, alpha] = getSegmentInfo(i, tGrid, xGrid, fGrid, t, idx)
    kLow = 2 * (i - 1) + 1;
    kMid = kLow + 1;
    kUpp = kLow + 2;
    h = tGrid(kUpp) - tGrid(kLow);
    xLow = xGrid(:, kLow);
    fLow = fGrid(:, kLow);
    fMid = fGrid(:, kMid);
    fUpp = fGrid(:, kUpp);
    alpha = t(idx) - tGrid(kLow);
end

function x = interpCubicSegment(h, xLow, fLow, fMid, fUpp, alpha)
    a = (2 * (fLow - 2 * fMid + fUpp)) / (3 * h^2);
    b = -(3 * fLow - 4 * fMid + fUpp) / (2 * h);
    c = fLow;
    d = xLow;

    x = d + alpha .* (c + alpha .* (b + alpha .* a));
end