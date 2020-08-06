
function c = bisection(f, a, b)

if f(a) * f(b) > 0
    
    % Produce warning
    error('Bad starting conditions')
    
else
    % Compute center point
    c = (a + b)/2;
    
    % Loop until tolerence is met
    dist = (b - a)/2;
    while dist > 1e-6
        
        % Set new limits
        if f(a)*f(c)<0
            b = c;
        else
            a = c;
        end
        
        % Compute new center and new distance
        c = (a + b)/2;
        dist = (b - a)/2;
    end
end
end