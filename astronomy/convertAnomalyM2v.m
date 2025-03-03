function v = convertAnomalyM2v(M, e)
    v = convertAnomalyE2v(convertAnomalyM2E(M, e), e);
end