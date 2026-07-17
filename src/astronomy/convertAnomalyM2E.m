function E = convertAnomalyM2E(M, e)
    func = @(E) E - e*sin(E) - M;
    E = fsolve(func, M, optimset('Display','off'));
end