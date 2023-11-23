importScripts("../util/manipulation_array.js");

onmessage = function(e)
{
    console.log("test");
    
    var data = JSON.parse(e.data);
    var s = data[0];
    var condition = data[1];
    var neighbors = new Array(s.length);

    // Création de la liste des premiers voisins
    for(i = 0; i < s.length; i++)
    {
        neighbors[i] = new Array(s[i].length);

        for(j = 0; j < s[i].length; j++)
        {
            neighbors[i][j] = new Array(s[i][j].length);

            for(k = 0; k < s[i][j].length; k++)
            {
                neighbors[i][j][k] = [];
            }
        }
    }

    if(condition == "strict")
    {
        for(i = 0; i < s.length; i++)
        {
            for(j = 0; j < s[i].length; j++)
            {
                for(k = 0; k < s[i][j].length; k++)
                {
                    var current_neighbors = [];

                    for(a = -1; a < 2; a += 1)
                    {
                        for(b = -1; b < 2; b += 1)
                        {
                            for(c = -1; c < 2; c += 1)
                            {
                                iprime = i + a;
                                jprime = j + b;
                                kprime = k + c;

                                var distance = Math.sqrt((iprime - i) * (iprime - i) + (jprime - j) * (jprime - j) + (kprime - k) * (kprime - k));

                                if(distance < 1.1 && distance > 0 && iprime >= 0 && iprime <= s.length-1 && jprime >= 0 && jprime <= s[i].length-1 && kprime >= 0 && kprime <= s[i][j].length-1)
                                {
                                    current_neighbors.push([iprime, jprime, kprime]);
                                }
                            }
                        }
                    }

                    neighbors[i][j][k] = current_neighbors;
                }
            }
        }
    }

    postMessage(JSON.stringify(neighbors));
}