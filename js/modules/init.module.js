atomistic.init =
{
    create_vector: function()
    {
        var arr = new Array(atomistic.constants.system.dx);
        
        for(var i = 0; i < arr.length; i++)
        {
            arr[i] = new Array(atomistic.constants.system.dy);

            for(var j = 0; j < arr[i].length; j++)
            {
                arr[i][j] = new Array(atomistic.constants.system.dz);

                for(var k = 0; k < arr[i][j].length; k++)
                {
                    arr[i][j][k] = new Array(atomistic.constants.system.dim);

                    for(var l = 0; l < 3; l++)
                    {
                        arr[i][j][k][l] = 0.0;
                    }
                }
            }
        }

        return arr;
    },
    
    create_correlation_function: function()
    {
        var arr = new Array(atomistic.constants.system.dx);
        
        for(var i = 0; i < arr.length; i++)
        {
            arr[i] = new Array(atomistic.constants.system.dy);

            for(var j = 0; j < arr[i].length; j++)
            {
                arr[i][j] = new Array(atomistic.constants.system.dz);

                for(var k = 0; k < arr[i][j].length; k++)
                {
                    arr[i][j][k] = new Array(atomistic.constants.system.dim);

                    for(var l = 0; l < 3; l++)
                    {
                        arr[i][j][k][l] = new Array(3);

                        for(var m = 0; m < 3; m++)
                        {
                            arr[i][j][k][l][m] = 0.0;
                        }
                    }
                }
            }
        }

        return arr;
    }
};