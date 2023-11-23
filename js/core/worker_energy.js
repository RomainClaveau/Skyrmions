importScripts("../util/manipulation_array.js");

var compute_Q = function(s, dx, dy)
{
    Q = 0.0

    for(a = 0; a < dx; a++)
    {
        for(b = 0; b < dy; b++)
        {
            dsdx = [0.0, 0.0, 0.0];
            dsdy = [0.0, 0.0, 0.0];

            if(a == 0) dsdx = matrix.add(s[a+1][b][0], matrix.multiply(s[a][b][0], -1));
            else if(a > 0 && a < dx - 1) dsdx = matrix.multiply(matrix.add(s[a+1][b][0], matrix.multiply(s[a-1][b][0], -1)), 0.5);
            else dsdx = matrix.add(s[a-1][b][0], matrix.multiply(s[a][b][0], -1));

            if(b == 0) dsdy = matrix.add(s[a][b+1][0], matrix.multiply(s[a][b][0], -1));
            else if(b > 0 && b < dy - 1) dsdy = matrix.multiply(matrix.add(s[a][b+1][0], matrix.multiply(s[a][b-1][0], -1)), 0.5);
            else dsdy = matrix.add(s[a][b-1][0], matrix.multiply(s[a][b][0], -1));

            Q += matrix.dot(s[a][b][0], matrix.cross(dsdx, dsdy));
        }
    }

    Q /= (4 * Math.PI);

    return Q;
};

var compute_mean_s = function(s, dx, dy, dz)
{
    sx = 0.0;
    sy = 0.0;
    sz = 0.0;

    for(a = 0; a < dx; a++)
    {
        for(b = 0; b < dy; b++)
        {
            for(c = 0; c < dz; c++)
            {
                sx += s[a][b][c][0] / (dx * dy * dz);
                sy += s[a][b][c][1] / (dx * dy * dz);
                sz += s[a][b][c][2] / (dx * dy * dz);
            }
        }
    }

    return [sx, sy, sz];
};

var compute_mean_norm = function(s, dx, dy, dz)
{
    var norm = 0.0;

    for(a = 0; a < dx; a++)
    {
        for(b = 0; b < dy; b++)
        {
            for(c = 0; c < dz; c++)
            {
                norm += Math.sqrt(s[a][b][c][0] * s[a][b][c][0] + s[a][b][c][1] * s[a][b][c][1] + s[a][b][c][2] * s[a][b][c][2]) / (dx * dy * dz);
            }
        }
    }

    return norm;
};

var compute_mean_W = function(W, dx, dy, dz)
{
    wx = 0.0;
    wy = 0.0;
    wz = 0.0;

    for(a = 0; a < dx; a++)
    {
        for(b = 0; b < dy; b++)
        {
            for(c = 0; c < dz; c++)
            {
                wx += W[a][b][c][0][0] / (dx * dy * dz);
                wy += W[a][b][c][1][1] / (dx * dy * dz);
                wz += W[a][b][c][2][2] / (dx * dy * dz);
            }
        }
    }

    return [wx, wy, wz];
};

var compute_mean_S = function(S, dx, dy, dz)
{
    sx = 0.0;
    sy = 0.0;
    sz = 0.0;

    for(a = 0; a < dx; a++)
    {
        for(b = 0; b < dy; b++)
        {
            for(c = 0; c < dz; c++)
            {
                sx += S[a][b][c][0][0] / (dx * dy * dz);
                sy += S[a][b][c][1][1] / (dx * dy * dz);
                sz += S[a][b][c][2][2] / (dx * dy * dz);
            }
        }
    }

    return [sx, sy, sz];
};

onmessage = function(e)
{
    var data = JSON.parse(e.data);
    
    var hbar = 6.582119514 * 1e-7; // Valeur de ħ en eV.s

    var s = data[0];
    var list_neighbors = data[1];
    var h = data[2];
    var constant_exchange = data[3];
    var constant_dmi = data[4];
    var B = data[5];
    var W = data[6];
    var S = data[7];

    // Calcul des champs du système
    for(a = 0; a < s.length; a++)
    {
        for(b = 0; b < s[a].length; b++)
        {
            for(c = 0; c < s[a][b].length; c++)
            {
                var l = list_neighbors[a][b][c];

                // On parcourt l'ensemble des voisins (pré-calculés) du site
                for(n = 0; n < l.length; n++)
                {
                    var aprime = list_neighbors[a][b][c][n][0];
                    var bprime = list_neighbors[a][b][c][n][1];
                    var cprime = list_neighbors[a][b][c][n][2];

                    var value_exchange = 0.5 * 1e9 * constant_exchange / hbar;

                    var value_dmi = 0.5 * 1e9 * constant_dmi / hbar;
                    
                    // Calcul de la pulsation du champ d'échange
                    h[a][b][c] = matrix.add(matrix.multiply(s[aprime][bprime][cprime], value_exchange), h[a][b][c]);

                    // Calcul de la pulsation du champ DMI
                    var r_ij = [a - aprime, b - bprime, c - cprime];
                    r_ij_norm = Math.sqrt(r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2]);
                    h[a][b][c] = matrix.add(matrix.multiply(matrix.cross(s[aprime][bprime][cprime], r_ij), value_dmi), h[a][b][c]);

                    // Calcul de la pulsation du champ de Zeeman
                    h[a][b][c] = matrix.add(h[a][b][c], [0.0, 0.0, 176e9 * B]);
                }
            }
        }
    }

    // Calcul de la charge topologique
    Q = compute_Q(s, s.length, s[0].length);
    sMean = compute_mean_s(s, s.length, s[0].length, s[0][0].length);
    norm = compute_mean_norm(s, s.length, s[0].length, s[0][0].length);
    WMean = compute_mean_W(W, W.length, W[0].length, W[0][0].length);
    SMean = compute_mean_S(S, S.length, S[0].length, S[0][0].length);

    postMessage(JSON.stringify([h, Q, sMean, norm, WMean, SMean]));
}