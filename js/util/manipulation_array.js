var matrix = {
    add: function(a, b)
    {
        var c = new Array(a.length);

        for(i = 0; i < a.length; i++)
        {
            if(Array.isArray(a[i]))
            {
                c[i] = new Array(a[i].length);

                for(j = 0; j < a[i].length; j++)
                {
                    if(Array.isArray(a[i][j]))
                    {
                        c[i][j] = new Array(a[i][j]);

                        for(k = 0; k < a[i][j].length; k++)
                        {
                            if(Array.isArray(a[i][j][k]))
                            {
                                c[i][j][k] = new Array(a[i][j][k].length);

                                for(l = 0; l < a[i][j][k].length; l++)
                                {
                                    if(Array.isArray(a[i][j][k][l]))
                                    {
                                        c[i][j][k][l] = new Array(a[i][j][k][l].length);

                                        for(m = 0; m < a[i][j][k][l].length; m++)
                                        {
                                            if(Array.isArray(a[i][j][k][l][m]))
                                            {
                                                // Do something...
                                            }
                                            else
                                            {
                                                c[i][j][k][l][m] = a[i][j][k][l][m] + b[i][j][k][l][m];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        c[i][j][k][l] = a[i][j][k][l] + b[i][j][k][l];
                                    }
                                }
                            }
                            else
                            {
                                c[i][j][k]= a[i][j][k] + b[i][j][k];
                            }
                        }
                    }
                    else
                    {
                        c[i][j]= a[i][j] + b[i][j];
                    }
                }
            }
            else
            {
                c[i] = a[i] + b[i];
            }
        }

        return c;
    },

    multiply: function(a, lambda)
    {
        var c = new Array(a.length);

        for(i = 0; i < a.length; i++)
        {
            if(Array.isArray(a[i]))
            {
                c[i] = new Array(a[i].length);

                for(j = 0; j < a[i].length; j++)
                {
                    if(Array.isArray(a[i][j]))
                    {
                        c[i][j] = new Array(a[i][j]);

                        for(k = 0; k < a[i][j].length; k++)
                        {
                            if(Array.isArray(a[i][j][k]))
                            {
                                c[i][j][k] = new Array(a[i][j][k].length);

                                for(l = 0; l < a[i][j][k].length; l++)
                                {
                                    if(Array.isArray(a[i][j][k][l]))
                                    {
                                        c[i][j][k][l] = new Array(a[i][j][k][l].length);

                                        for(m = 0; m < a[i][j][k][l].length; m++)
                                        {
                                            if(Array.isArray(a[i][j][k][l][m]))
                                            {
                                                // Do something...
                                            }
                                            else
                                            {
                                                c[i][j][k][l][m] = lambda * a[i][j][k][l][m];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        c[i][j][k][l] = lambda * a[i][j][k][l];
                                    }
                                }
                            }
                            else
                            {
                                c[i][j][k] = lambda * a[i][j][k];
                            }
                        }
                    }
                    else
                    {
                        c[i][j] = lambda * a[i][j];
                    }
                }
            }
            else
            {
                c[i] = lambda * a[i];
            }
        }

        return c;
    },

    cross: function(a, b)
    {
        var c = new Array(a.length);

        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];

        return c;
    },

    dot: function(a, b)
    {
        var c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

        return c;
    }
} || {};