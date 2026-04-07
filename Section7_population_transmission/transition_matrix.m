function S = transition_matrix(NtoI, NtoV,IprimetoI,IprimetoV,VprimetoI,VprimetoV,IprimetoN, VprimetoN)
    S = zeros(5,5,"gpuArray");
    S(2,1) = NtoI;
    S(3,1) = NtoV;
    S(1,1) = 1-S(2,1)-S(3,1);

    S(4,2) = 1;
    S(5,3) = 1;

    S(1,4) = IprimetoN;
    S(2,4) = IprimetoI;
    S(3,4) = IprimetoV;
    S(4,4) = 1-S(2,4)-S(3,4)-S(1,4);

    S(1,5) = VprimetoN;
    S(2,5) = VprimetoI;
    S(3,5) = VprimetoV;
    S(5,5) = 1-S(2,5)-S(3,5)-S(1,5);
    
end