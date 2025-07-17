function[resultt] = clebsch_gordan(j_1, j_2, j_3, m_1, m_2, m_3)
    res1 = (-1) ^ (j_1 - j_2 + m_3) * sqrt(2 * j_3 + 1) ;
    % first claculte the factorial list needed before hand to fast the code
    nn = max([j_1 + j_2 + j_3 + 1, j_1 + abs(m_1), j_2 + abs(m_2), j_3 + abs(m_3)]);
    Factlist = 1;
    if nn >= numel(Factlist)
        for ii = numel(Factlist):(nn + 1)
            Factlist = [Factlist, Factlist(end) * ii];
        end  % end for ii
    end % end nn
    % :) calculated all the factorial list

    m_3 = -m_3;
    % chech if j and m values given makes any sense or nor
    if mod(j_1 * 2, 1) || mod(j_2 * 2, 1) || mod(j_3 * 2, 1)
        error('j values must be integer or half integer');
    end

    if mod(m_1 * 2, 1) || mod(m_2 * 2, 1) || mod(m_3 * 2, 1)
        error('m values must be integer or half integer');
    end
    % checked all j and m values

    % Sum of all angular momentum should be zero already taken m3 = -m3 so
    % m1+m2=m3
    if m_1 + m_2 + m_3 ~= 0
        resultt = 0;
        return;
    end

    prefid = (-1) ^ (j_1 - j_2 - m_3);
    %m_3 = -m_3;
    a1 = j_1 + j_2 - j_3;
    if a1 < 0
        resultt = 0;
        return;
    end

    a2 = j_1 - j_2 + j_3;
    if a2 < 0
        resultt = 0;
        return;
    end

    a3 = -j_1 + j_2 + j_3;
    if a3 < 0
        resultt = 0;
        return;
    end

    if abs(m_1) > j_1 || abs(m_2) > j_2 || abs(m_3) > j_3
        resultt = 0;
        return;
    end

    argsqrt = Factlist(j_1 + j_2 - j_3 + 1) * Factlist(j_1 - j_2 + j_3 + 1) * Factlist(-j_1 + j_2 + j_3 + 1) * Factlist(j_1 - m_1 + 1) * ...
    Factlist(j_1 + m_1 + 1) * Factlist(j_2 - m_2 + 1) * Factlist(j_2 + m_2 + 1) * Factlist(j_3 - m_3 + 1) * Factlist(j_3 + m_3 + 1) ...
    / Factlist(j_1 + j_2 + j_3 + 2);

    ressqrt = sqrt(argsqrt);
    imin = max([-j_3 + j_1 + m_2, -j_3 + j_2 - m_1, 0]);
    imax = min([j_2 + m_2, j_1 - m_1, j_1 + j_2 - j_3]);
    sumres = 0;

    for ii = imin:imax
        den = Factlist(ii + 1) * Factlist(ii + j_3 - j_1 - m_2 + 1) * Factlist(j_2 + m_2 - ii + 1) * Factlist(j_1 - ii - m_1 + 1) * ...
            Factlist(ii + j_3 - j_2 + m_1 + 1) * Factlist(j_1 + j_2 - j_3 - ii + 1);
        sumres = sumres + (-1) ^ ii / den;
    end


    resultt = res1 * ressqrt * sumres * prefid;



end
