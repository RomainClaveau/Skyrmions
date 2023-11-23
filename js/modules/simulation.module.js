atomistic.simulation =
{
    updateValues: function()
    {
        document.querySelectorAll("#values_Q input")[0].value = atomistic.constants.observables.Q;
        document.querySelectorAll("#values_mean input")[0].value = atomistic.constants.observables.sMean[0];
        document.querySelectorAll("#values_mean input")[1].value = atomistic.constants.observables.sMean[1];
        document.querySelectorAll("#values_mean input")[2].value = atomistic.constants.observables.sMean[2];
        document.querySelectorAll("#values_mean2 input")[0].value = atomistic.constants.observables.sNorm;
        document.querySelectorAll("#values_mean3 input")[0].value = atomistic.constants.observables.WMean[0];
        document.querySelectorAll("#values_mean3 input")[1].value = atomistic.constants.observables.WMean[1];
        document.querySelectorAll("#values_mean3 input")[2].value = atomistic.constants.observables.WMean[2];
        document.querySelectorAll("#values_mean4 input")[0].value = atomistic.constants.observables.SMean[0];
        document.querySelectorAll("#values_mean4 input")[1].value = atomistic.constants.observables.SMean[1];
        document.querySelectorAll("#values_mean4 input")[2].value = atomistic.constants.observables.SMean[2];
        document.querySelectorAll("#values_convergence input")[0].value = atomistic.constants.studies.currentStudy.delta;
        

        document.querySelectorAll("#chronology b")[0].innerHTML = atomistic.constants.simulation.results.length;
        document.querySelector("#chronology input").max = atomistic.constants.simulation.results.length - 1;

        document.querySelector(".progress").innerHTML = atomistic.constants.simulation.counter + " pas calculés";
    },

    saveResult: function()
    {
        var current_s = JSON.parse(JSON.stringify(atomistic.constants.observables.s));
        atomistic.constants.simulation.results.push([atomistic.constants.simulation.counter, current_s]);
    },

    printResult: function()
    {
        //document.querySelector("textarea").innerHTML = atomistic.constants.simulation.results;
    },
    
    simulate: function()
    {
        if(atomistic.constants.simulation.debug) console.log("Initialisation : " + window.performance.now());
        
        atomistic.simulation.reinit();
        
        if(atomistic.constants.studies.study) atomistic.studies.set_values_of_study();
        if(atomistic.constants.studies.study && !atomistic.constants.studies.initialized)
        {
            atomistic.studies.initSavedValues();
            atomistic.constants.studies.initialized = true;
        }

        // Création du vecteur de l'aimantation ainsi que les fonctions de corrélation associées au bruit thermique
        atomistic.constants.observables.s = atomistic.init.create_vector();
        atomistic.constants.observables.W = atomistic.init.create_correlation_function();
        atomistic.constants.observables.S = atomistic.init.create_correlation_function();
        atomistic.constants.simulation.struct = atomistic.init.create_vector();

        if(atomistic.constants.simulation.debug) console.log("Création de la configuration : " + window.performance.now());
        
        // Création de la configuration initiale de l'aimantation
        if(atomistic.constants.studies.initialConfig.s == null)
        {
            atomistic.configuration.down();
            atomistic.configuration.compute_correlation_S();
            atomistic.constants.studies.initialConfig.s = JSON.parse(JSON.stringify(atomistic.constants.observables.s));
            atomistic.constants.studies.initialConfig.S = JSON.parse(JSON.stringify(atomistic.constants.observables.S));
        }
        else
        {
            atomistic.constants.observables.s = JSON.parse(JSON.stringify(atomistic.constants.studies.initialConfig.s));
            atomistic.constants.observables.S = JSON.parse(JSON.stringify(atomistic.constants.studies.initialConfig.S));
        }
        
        // Mise à jour des vecteurs sur la grille
        update_render(atomistic.constants.observables.s, atomistic.constants.renderer.scene, atomistic.constants.renderer.camera, atomistic.constants.system.dx, atomistic.constants.system.dy, atomistic.constants.system.dz);

        if(atomistic.constants.simulation.debug) console.log("Création de la liste des tâches : " + window.performance.now());

        // Création du "pool"
        atomistic.multithreading.createPool();

        if(atomistic.constants.simulation.debug) console.log("Création des workers : " + window.performance.now());

        // Création des workers
        atomistic.multithreading.createWorkers();

        if(atomistic.constants.simulation.debug) console.log("Création de la liste des voisins : " + window.performance.now());

        // Création de la liste des voisins
        var worker_neighbors = new Worker("js/core/worker_neighbors.js");
        var content = JSON.stringify([atomistic.constants.observables.s, atomistic.constants.simulation.neighbors])
        worker_neighbors.postMessage(content);

        worker_neighbors.onmessage = function(e)
        {
            var data = JSON.parse(e.data);
            atomistic.constants.observables.listNeighbors = data;
            
            e.target.terminate();

            // La liste des voisins créée, nous pouvons débuter la simulation
            atomistic.simulation.process();
        };
        
    },
    
    reinit: function()
    {
        atomistic.constants.observables.Q = 0;
        atomistic.constants.observables.Q_sum = 0;
        atomistic.constants.observables.S = null;
        atomistic.constants.observables.W = null;
        atomistic.constants.observables.s = null;
        atomistic.constants.observables.h = null;
        atomistic.constants.simulation.counter = 0;
        atomistic.constants.observables.listNeighbors = null;
        atomistic.constants.simulation.results = [];
    },
    
    process: function()
    {
        if(atomistic.constants.simulation.debug) console.log("Calcul de l'énergie du système : " + window.performance.now());

        // Calcul de l'énergie
        var worker_energy = new Worker("js/core/worker_energy.js");
        var content = JSON.stringify([atomistic.constants.observables.s, atomistic.constants.observables.listNeighbors, atomistic.constants.simulation.struct, atomistic.constants.physicalParameters.J, atomistic.constants.physicalParameters.K, atomistic.constants.physicalParameters.B, atomistic.constants.observables.W, atomistic.constants.observables.S]);
        worker_energy.postMessage(content);

        worker_energy.onmessage = function(e)
        {
            var data = JSON.parse(e.data);
            
            atomistic.constants.observables.h = data[0];
            atomistic.constants.observables.Q = data[1];
            atomistic.constants.observables.sMean = data[2];
            atomistic.constants.observables.sNorm = data[3];
            atomistic.constants.observables.WMean = data[4];
            atomistic.constants.observables.SMean = data[5];

            console.log("Début de l'intégration numérique : " + window.performance.now());

            atomistic.constants.simulation.pause = false;
            document.querySelector("#state").innerHTML = `Etude pour D = ${atomistic.constants.physicalParameters.D} et B = ${atomistic.constants.physicalParameters.B}`;
            
            e.target.terminate();
        
            for(var i = 0; i < atomistic.constants.multithreading.nbWorkers; i++) atomistic.multithreading.compute_pool_element(i);
        };
    },

    update_pause: function()
    {
        if(atomistic.constants.simulation.pause === true)
        {
            atomistic.constants.simulation.pause = false;
            atomistic.simulation.process();
            document.querySelector(".stop_simu input").value = "Mettre en pause la simulation";
            document.querySelector("#state").innerHTML = "Simulation en cours...";
        }
        else
        {
            atomistic.constants.simulation.pause = true;
            document.querySelector(".stop_simu input").value = "Continuer la simulation";
            document.querySelector("#state").innerHTML = "Simulation en pause";
        }
    },
    
    update_interface: function()
    {
        document.querySelectorAll("#values_dimension input")[0].value = atomistic.constants.system.dx;
        document.querySelectorAll("#values_dimension input")[1].value = atomistic.constants.system.dy;
        document.querySelectorAll("#values_dimension input")[2].value = atomistic.constants.system.dz;
        
        document.querySelectorAll("#values_exch input")[0].value = atomistic.constants.physicalParameters.J * 1e3;
        document.querySelectorAll("#values_dmi input")[0].value = atomistic.constants.physicalParameters.K * 1e3;
        document.querySelectorAll("#values_dd input")[0].value = atomistic.constants.physicalParameters.L * 1e3;
        document.querySelectorAll("#values_zee input")[0].value = atomistic.constants.physicalParameters.B;
        
        document.querySelectorAll("#values_alpha input")[0].value = atomistic.constants.physicalParameters.α;
        
        document.querySelectorAll("#values_noise input")[0].value = atomistic.constants.physicalParameters.D * 1e-14;

        document.querySelectorAll("#values_dt input")[0].value = atomistic.constants.physicalParameters.dt * 1e15;
        document.querySelectorAll("#values_tau input")[0].value = atomistic.constants.physicalParameters.τ * 1e15;
        document.querySelectorAll("#values_workers input")[0].value = atomistic.constants.multithreading.nbWorkers;

        document.querySelectorAll("#values_live input")[0].checked = atomistic.constants.renderer.live;
        document.querySelectorAll("#values_quality input")[0].value = atomistic.constants.renderer.quality;
        document.querySelectorAll("#values_steps input")[0].value = atomistic.constants.simulation.updateEveryStep;
        document.querySelectorAll("#values_nuances input")[0].value = atomistic.constants.renderer.nuances;

        document.querySelectorAll("#values_Q input")[0].value = 0.0;
        document.querySelectorAll("#values_mean input")[0].value = 0.0;
        document.querySelectorAll("#values_mean input")[1].value = 0.0;
        document.querySelectorAll("#values_mean input")[2].value = 0.0;

        document.querySelectorAll("#values_mean2 input")[0].value = 1.0;

        document.querySelectorAll("#values_mean3 input")[0].value = 0.0;
        document.querySelectorAll("#values_mean3 input")[1].value = 0.0;
        document.querySelectorAll("#values_mean3 input")[2].value = 0.0;

        document.querySelectorAll("#values_mean4 input")[0].value = 0.0;
        document.querySelectorAll("#values_mean4 input")[1].value = 0.0;
        document.querySelectorAll("#values_mean4 input")[2].value = 0.0;
        
        document.querySelectorAll("#values_convergence input")[0].value = 0.0;

        document.querySelector(".progress").innerHTML = atomistic.constants.simulation.counter + " pas calculés";
    },

    viewChronology: function(frame)
    {
        // On ne met plus à jour l'interface
        atomistic.constants.renderer.live = false;
        document.querySelectorAll("#values_live input")[0].checked = atomistic.constants.renderer.live;

        // On affiche la scène selectionnée
        document.querySelectorAll("#chronology b")[1].innerHTML = "Image <i>" + frame + "</i> sélectionnée";
        var frame_s = atomistic.constants.simulation.results[frame][1];
        update_render(frame_s, atomistic.constants.renderer.scene, atomistic.constants.renderer.camera, atomistic.constants.system.dx, atomistic.constants.system.dy, atomistic.constants.system.dz);
    }
};