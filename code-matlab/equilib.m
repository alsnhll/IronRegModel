function eq = equilib(fp, x0)

MAX_TIME         = 50000;
STEADY_STATE_TOL = 1e-5;
TIME_STEP        = 250;
LOOK_BACK        = ceil((1/5) * TIME_STEP);

start_t = 0;
end_t   = TIME_STEP;

while end_t < MAX_TIME

    [T, Y] = ode23(fp, start_t:end_t, x0);
    vq = Y((end - LOOK_BACK):end, :);
    fq = vq(1, :);
    eq = vq(end, :);

    if sqrt(sum((fq - eq) .* (fq - eq))) > STEADY_STATE_TOL
        start_t = end_t;
        end_t   = start_t + TIME_STEP;
        x0      = eq;
        fprintf('Extending time. Start time: %d. End time: %d.\n', start_t, end_t);
    else
        fprintf('Time to steady-state: %d.\n', end_t);
        break;
    end

end

if end_t >= MAX_TIME
     fprintf('Warning: MAX_TIME reached.\n');
     for i = (LOOK_BACK - 1):-1:1
         l2d = sqrt(sum(vq(i, :) .* vq(end, :)));
         fprintf('L2 distance to -%d vector: %f\n', i, l2d);
     end
end
