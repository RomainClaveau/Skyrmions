atomistic.configuration =
{
    random: function()
    {
        for(var i = 0; i < atomistic.constants.system.dx; i++)
        {
            for(var j = 0; j < atomistic.constants.system.dy; j++)
            {
                for(var k = 0; k < atomistic.constants.system.dz; k++)
                {
                    var rd1 = Math.random() * 2 * Math.PI;
                    var rd2 = Math.random() * 2 * Math.PI;
                    
                    atomistic.constants.observables.s[i][j][k][0] = Math.sin(rd1) * Math.cos(rd2);
                    atomistic.constants.observables.s[i][j][k][1] = Math.sin(rd1) * Math.sin(rd2);
                    atomistic.constants.observables.s[i][j][k][2] = Math.cos(rd1);
                }
            }
        }
    },

    down: function()
    {
        for(var i = 0; i < atomistic.constants.system.dx; i++)
        {
            for(var j = 0; j < atomistic.constants.system.dy; j++)
            {
                for(var k = 0; k < atomistic.constants.system.dz; k++)
                {                        
                    atomistic.constants.observables.s[i][j][k][0] = 0;
                    atomistic.constants.observables.s[i][j][k][1] = 0;
                    atomistic.constants.observables.s[i][j][k][2] = -1;
                }
            }
        }
    },
    
    compute_correlation_S: function()
    {
        for(var i = 0; i < atomistic.constants.system.dx; i++)
        {
            for(var j = 0; j < atomistic.constants.system.dy; j++)
            {
                for(var k = 0; k < atomistic.constants.system.dz; k++)
                {
                    for(var l = 0; l < atomistic.constants.system.dim; l++)
                    {
                        for(var m = 0; m < atomistic.constants.system.dim; m++)
                        {
                            atomistic.constants.observables.S[i][j][k][l][m] = atomistic.constants.observables.s[i][j][k][l] * atomistic.constants.observables.s[i][j][k][m];
                        }
                    }
                }
            }
        }
    }
};