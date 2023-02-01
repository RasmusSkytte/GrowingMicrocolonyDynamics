
function c = bisection(f, a, b)

fa = f(a);
fb = f(b);

if fa * fb > 0
    
    % Produce warning
    error('Bad starting conditions')
    
else
    
    % Compute center point
    c = (a + b)/2;
    
    % Loop until tolerence is met
    fc = f(c);
    
    while abs(fc) > 1e-2
        
        % Set new limits
        if fa*fc<0
            b = c;
        else
            a = c;
            fa = fc;
        end
        
        % Compute new center and new distance
        c = (a + b)/2;
        fc = f(c);
    end
end
end