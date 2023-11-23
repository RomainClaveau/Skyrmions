atomistic.multithreading =
{

    createPool: function()
    {
        atomistic.constants.multithreading.currentPool = new Array(atomistic.constants.system.dx * atomistic.constants.system.dy * atomistic.constants.system.dz);
        var counter = 0;
        
        for(var i = 0; i < atomistic.constants.system.dx; i++)
        {
            for(var j = 0; j < atomistic.constants.system.dy; j++)
            {
                for(var k = 0; k < atomistic.constants.system.dz; k++)
                {
                    atomistic.constants.multithreading.currentPool[counter] = [i, j, k, atomistic.constants.observables.s[i][j][k], atomistic.constants.observables.W[i][j][k], atomistic.constants.observables.S[i][j][k]];

                    counter++;
                }
            }
        }
    },

    updatePool: function()
    {
        var counter = 0;
        
        for(var i = 0; i < atomistic.constants.system.dx; i++)
        {
            for(var j = 0; j < atomistic.constants.system.dy; j++)
            {
                for(var k = 0; k < atomistic.constants.system.dz; k++)
                {
                    atomistic.constants.multithreading.currentPool[counter] = [i, j, k, atomistic.constants.observables.s[i][j][k], atomistic.constants.observables.W[i][j][k], atomistic.constants.observables.S[i][j][k]];

                    counter++;
                }
            }
        }
    },

    createWorkers: function()
    {
        atomistic.constants.multithreading.workers = new Array(atomistic.constants.multithreading.nbWorkers);

        for(var i = 0; i < atomistic.constants.multithreading.workers.length; i++)
        {
            atomistic.constants.multithreading.workers[i] = new Worker('js/core/worker_RK4.js');
            atomistic.constants.multithreading.workers[i].onmessage = function(e) { atomistic.multithreading.manage_worker(e); };

            console.log(atomistic.constants.multithreading.workers);
        }
    },
    
    compute_pool_element: function(i)
    {        
        // Passage du contenu au worker
        var content = JSON.stringify([i, atomistic.constants.multithreading.currentPool, atomistic.constants.multithreading.nbWorkers, atomistic.constants.physicalParameters.α, atomistic.constants.physicalParameters.τ, atomistic.constants.physicalParameters.dt, atomistic.constants.physicalParameters.D, atomistic.constants.observables.h]);
        atomistic.constants.multithreading.workers[i].postMessage(content);
    },
    
    manage_worker: function(e)
    {
        // Mise à jour de l'étude
        if(atomistic.constants.simulation.counter % 50 == 0) atomistic.studies.update_study();
        
        var data = JSON.parse(e.data);
        
        atomistic.constants.simulation.counter += 1/atomistic.constants.multithreading.nbWorkers;

        // On met à jour les éléments
        var workerId = data[0];

        for(var site = workerId; site < data[1].length; site += atomistic.constants.multithreading.nbWorkers)
        {
            var x = data[1][site][0];
            var y = data[1][site][1];
            var z = data[1][site][2];

            atomistic.constants.observables.s[x][y][z] = matrix.add(atomistic.constants.observables.s[x][y][z], data[1][site][3]);
            atomistic.constants.observables.W[x][y][z] = matrix.add(atomistic.constants.observables.W[x][y][z], data[1][site][4]);
            atomistic.constants.observables.S[x][y][z] = matrix.add(atomistic.constants.observables.S[x][y][z], data[1][site][5]);
        }

        atomistic.multithreading.updatePool();

        if(atomistic.constants.simulation.debug) console.log("Mise à jour des valeurs : " + window.performance.now());

        // Si tous les processus ont terminé leurs calculs, nous pouvons effectuer un incrément temporel
        if(parseInt(atomistic.constants.simulation.counter) == atomistic.constants.simulation.counter)
        {
            if(atomistic.constants.simulation.debug) console.log("Incrément temporel : " + window.performance.now());

            if(atomistic.constants.simulation.counter % atomistic.constants.simulation.updateEveryStep == 0) atomistic.simulation.saveResult();

            if(atomistic.constants.simulation.counter % 1 == 0) atomistic.simulation.updateValues();

            if(atomistic.constants.renderer.live && atomistic.constants.simulation.counter % atomistic.constants.simulation.updateEveryStep == 0) update_render(atomistic.constants.observables.s, atomistic.constants.renderer.scene, atomistic.constants.renderer.camera, atomistic.constants.system.dx, atomistic.constants.system.dy, atomistic.constants.system.dz);

            // Nouveau calcul de l'énergie
            if(atomistic.constants.simulation.pause === false)
            {
                if(atomistic.constants.multithreading.worker_energy == null) atomistic.constants.multithreading.worker_energy = new Worker("js/core/worker_energy.js");
                var content = JSON.stringify([atomistic.constants.observables.s, atomistic.constants.observables.listNeighbors, atomistic.constants.simulation.struct, atomistic.constants.physicalParameters.J, atomistic.constants.physicalParameters.K, atomistic.constants.physicalParameters.B, atomistic.constants.observables.W, atomistic.constants.observables.S]);
                
                atomistic.constants.multithreading.worker_energy.postMessage(content);

                atomistic.constants.multithreading.worker_energy.onmessage = function(e)
                {
                    var data = JSON.parse(e.data);
                    
                    atomistic.constants.studies.currentStudy["delta"] = Math.abs(atomistic.constants.observables.Q - data[1]);
                    
                    atomistic.constants.observables.h = data[0];
                    atomistic.constants.observables.Q = data[1];
                    atomistic.constants.observables.Q_sum += data[1];
                    atomistic.constants.observables.sMean = data[2];
                    atomistic.constants.observables.sNorm = data[3];
                    atomistic.constants.observables.WMean = data[4];
                    atomistic.constants.observables.SMean = data[5];
                
                    for(var i = 0; i < atomistic.constants.multithreading.nbWorkers; i++) atomistic.multithreading.compute_pool_element(i);
                };
            }
        }
    }
};