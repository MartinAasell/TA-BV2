function yout = SIRmodel(t, y, N, mu, beta, v)

yout = [ mu*N - mu*y(1) - beta*y(2)/N*y(1);
        beta*y(2)/N*y(1) - v*y(2) - mu*y(2);
        v*y(2) - mu*y(3)];