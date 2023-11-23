importScripts("../util/manipulation_array.js");

var epsilon = [
    [
        [0, 0, 0],
        [0, 0, 1],
        [0, -1, 0]
    ],
    [
        [0, 0, -1],
        [0, 0, 0],
        [1, 0, 0]
    ],
    [
        [0, 1, 0],
        [-1, 0, 0],
        [0, 0, 0]
    ]
];

F = function(s, h, W, S, tau, coeff, alpha, D)
{
    value = new Array(3);

    for(i = 0; i < 3; i++) value[i] = 0.0;

    for(i = 0; i < 3; i++)
    {
        for(k = 0; k < 3; k++)
        {
            value[i] += coeff * alpha * (h[i] * S[k][k] - h[k] * S[k][i]);

            for(j = 0; j < 3; j++)
            {
                value[i] += coeff * (epsilon[i][j][k] * h[j] * s[k] + epsilon[i][j][k] * W[j][k]);
            }
        }
    }

    return value
},

G = function(s, h, W, S, tau, coeff, alpha, D)
{
    value = Array.from(Array(3), () => new Array(3));

    for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) value[i][j] = 0.0;

    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            value[i][j] -= W[i][j] / tau;

            for(l = 0; l < 3; l++)
            {
                value[i][j] += coeff * D * epsilon[j][i][l] * s[l] / tau - coeff * alpha * (h[l] * (s[l] * W[i][j] + s[j] * W[i][l]) - 2 * h[j] * s[l] * W[i][l]);

                for(k = 0; k < 3; k++)
                {
                    value[i][j] += coeff * (epsilon[j][k][l] * h[k] * W[i][l]);
                }
            }
        }
    }

    return value;                   
},

H = function(s, h, W, S, tau, coeff, alpha, D)
{
    value = Array.from(Array(3), () => new Array(3));

    for(i = 0; i < 3; i++) for(j = 0; j < 3; j++) value[i][j] = 0.0;

    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            for(l = 0; l < 3; l++)
            {
                value[i][j] -= coeff * alpha * (h[l] * (s[i] * S[l][j] + s[l] * S[i][j] + s[j] * S[i][l] - 2 * s[i] * s[j] * s[l]) - h[j] * (s[i] * S[l][l] + 2 * s[l] * S[i][l] - 2 * s[i] * s[l] * s[l]) + h[l] * (s[j] * S[l][i] + s[l] * S[j][i] + s[i] * S[j][l] - 2 * s[j] * s[i] * s[l]) - h[i] * (s[j] * S[l][l] + 2 * s[l] * S[j][l] - 2 * s[j] * s[l] * s[l]));
                for(k = 0; k < 3; k++)
                {
                    value[i][j] += coeff * (epsilon[j][k][l] * h[k] * S[i][l] + epsilon[j][k][l] * (s[i] * W[k][l] + s[l] * W[k][i]) + epsilon[i][k][l] * h[k] * S[j][l] + epsilon[i][k][l] * (s[j] * W[k][l] + s[l] * W[k][j]));
                }
            }
        }
    }

    return value;
},

onmessage = function(e)
{
    var data = JSON.parse(e.data);
    
    if(data != undefined && data[0] != undefined && data[1] != undefined)
    {
        var workerId = data[0];
        var pool = data[1];
        var nbWorkers = data[2];
        var alpha = data[3];
        var tau = data[4];
        var dt = data[5];
        var D = data[6];
        var coeff = 1 / (1 + alpha * alpha);

        var results = new Array(pool.length);
        
        for(site = workerId; site < pool.length; site += nbWorkers)
        {
            var x = pool[site][0];
            var y = pool[site][1];
            var z = pool[site][2];
            var s = pool[site][3];
            var W = pool[site][4];
            var S = pool[site][5];
            var h = data[7][x][y][z];

            // Méthode de Runge-Kutta explicite d'ordre 4
            a1 = matrix.multiply(F(s, h, W, S, tau, coeff, alpha, D), dt);
            b1 = matrix.multiply(G(s, h, W, S, tau, coeff, alpha, D), dt);
            c1 = matrix.multiply(H(s, h, W, S, tau, coeff, alpha, D), dt);

            a2 = matrix.multiply(F(matrix.add(s, matrix.multiply(a1, 0.5)), h, matrix.add(W, matrix.multiply(b1, 0.5)), matrix.add(S, matrix.multiply(c1, 0.5)), tau, coeff, alpha, D), dt);
            b2 = matrix.multiply(G(matrix.add(s, matrix.multiply(a1, 0.5)), h, matrix.add(W, matrix.multiply(b1, 0.5)), matrix.add(S, matrix.multiply(c1, 0.5)), tau, coeff, alpha, D), dt);
            c2 = matrix.multiply(H(matrix.add(s, matrix.multiply(a1, 0.5)), h, matrix.add(W, matrix.multiply(b1, 0.5)), matrix.add(S, matrix.multiply(c1, 0.5)), tau, coeff, alpha, D), dt);

            a3 = matrix.multiply(F(matrix.add(s, matrix.multiply(a2, 0.5)), h, matrix.add(W, matrix.multiply(b2, 0.5)), matrix.add(S, matrix.multiply(c2, 0.5)), tau, coeff, alpha, D), dt);
            b3 = matrix.multiply(G(matrix.add(s, matrix.multiply(a2, 0.5)), h, matrix.add(W, matrix.multiply(b2, 0.5)), matrix.add(S, matrix.multiply(c2, 0.5)), tau, coeff, alpha, D), dt);
            c3 = matrix.multiply(H(matrix.add(s, matrix.multiply(a2, 0.5)), h, matrix.add(W, matrix.multiply(b2, 0.5)), matrix.add(S, matrix.multiply(c2, 0.5)), tau, coeff, alpha, D), dt);

            a4 = matrix.multiply(F(matrix.add(s, a3), h, matrix.add(W, b3), matrix.add(S, c3), tau, coeff, alpha, D), dt);
            b4 = matrix.multiply(G(matrix.add(s, a3), h, matrix.add(W, b3), matrix.add(S, c3), tau, coeff, alpha, D), dt);
            c4 = matrix.multiply(H(matrix.add(s, a3), h, matrix.add(W, b3), matrix.add(S, c3), tau, coeff, alpha, D), dt);

            add_s = matrix.multiply(matrix.add(a1, matrix.add(matrix.multiply(a2, 2),matrix.add(matrix.multiply(a3, 2), a4))), 0.166666666667);
            add_W = matrix.multiply(matrix.add(b1, matrix.add(matrix.multiply(b2, 2),matrix.add(matrix.multiply(b3, 2), b4))), 0.166666666667);
            add_S = matrix.multiply(matrix.add(c1, matrix.add(matrix.multiply(c2, 2),matrix.add(matrix.multiply(c3, 2), c4))), 0.166666666667);

            results[site] = [x, y, z, add_s, add_W, add_S];
        }
    
        postMessage(JSON.stringify([workerId, results]));
    }
}