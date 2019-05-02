function [ ya_ddot ] = get_ya_ddot_lowgain_ff(ya, ya_dot, y, ydot, yddot)

    kp = 25;
    kv = 10;
    ya_ddot = kp*(y-ya) + kv*(ydot-ya_dot) + yddot;
    
end

