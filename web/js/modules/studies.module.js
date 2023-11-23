atomistic.studies =
{
    update_study: function()
    {
        // Si ΔQ < 1e-4
        if(atomistic.constants.studies.currentStudy["delta"] < atomistic.constants.studies.currentStudy.maxDelta || atomistic.constants.simulation.counter >= atomistic.constants.studies.currentStudy.maxSteps)
        {
            // On met la simulation en pause
            atomistic.constants.simulation.pause = true;

            // On ajoute le point
            atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.zAxis][atomistic.constants.studies.currentStudy.counterY][atomistic.constants.studies.currentStudy.counterX] = JSON.parse(JSON.stringify(atomistic.constants.observables[atomistic.constants.studies.currentStudy.zAxis])) / JSON.parse(JSON.stringify(atomistic.constants.simulation.counter));

            // Si l'on est sur l'onglet "2D", on met le graphique à jour
            if(atomistic.constants.renderer.currentTab == "2d") atomistic.studies.updateGraphics();

            // On incrémente les valeurs de l'étude
            atomistic.constants.studies.currentStudy.x += atomistic.constants.studies.currentStudy.dx;
            atomistic.constants.studies.currentStudy.counterX += 1;


            // Si x > xmax
            if(atomistic.constants.studies.currentStudy.x > atomistic.constants.studies.currentStudy.xmax)
            {                    
                // Incrémentation des valeurs des paramètres
                atomistic.constants.studies.currentStudy.x = 0;
                atomistic.constants.studies.currentStudy.y += atomistic.constants.studies.currentStudy.dy;

                atomistic.constants.studies.currentStudy.counterX = 0;
                atomistic.constants.studies.currentStudy.counterY += 1;

                // Si y > ymax
                if(atomistic.constants.studies.currentStudy.y > atomistic.constants.studies.currentStudy.ymax)
                {
                    // L'étude est terminée
                    atomistic.studies.stop_study();
                }
                else
                {
                    // On fait une nouvelle simulation avec les nouveaux paramètres
                    atomistic.simulation.simulate();
                }
            }
            else
            {
                // On fait une nouvelle simulation avec les nouveaux paramètres
                atomistic.simulation.simulate();
            }
        }
    },

    stop_study: function()
    {
        for(var i = 0; i < atomistic.constants.multithreading.nbWorkers; i++)
        {
            atomistic.constants.multithreading.workers[i].terminate();
        }
    },

    set_values_of_study: function()
    {
        atomistic.constants.physicalParameters[atomistic.constants.studies.currentStudy.xAxis] = atomistic.constants.studies.currentStudy.x;
        atomistic.constants.physicalParameters[atomistic.constants.studies.currentStudy.yAxis] = atomistic.constants.studies.currentStudy.y;
    },

    initSavedValues: function()
    {
        if(atomistic.constants.studies.dim == "3d")
        {
            for(var i = 0; i < atomistic.constants.studies.currentStudy.ymax; i += atomistic.constants.studies.currentStudy.dy)
            {
                atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.yAxis].push(i);
            }

            for(var i = 0; i < atomistic.constants.studies.currentStudy.xmax; i += atomistic.constants.studies.currentStudy.dx)
            {
                atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.xAxis].push(i);
            }

            for(var i = 0; i < atomistic.constants.studies.currentStudy.ymax; i += atomistic.constants.studies.currentStudy.dy)
            {
                var a = new Array();

                for(var j = 0; j < atomistic.constants.studies.currentStudy.xmax; j += atomistic.constants.studies.currentStudy.dx)
                {
                    a.push(null);
                }

                atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.zAxis].push(a);
            }
        }
    },

    updateGraphics: function()
    {
        var graph = document.querySelector("#graph_Q_B_T");
        
        graph.innerHTML = "";

        var dim = atomistic.constants.studies.dim;
        var xValues = atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.xAxis];
        var yValues = atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.yAxis];
        
        if(atomistic.constants.studies.initialized)
        {
            if(dim == "2d")
            {
                var graph_content =
                {
                    x: xValues,
                    y: yValues,
                    type: 'lines+markers'
                };

                var graph_layout =
                {
                    title: atomistic.constants.studies.currentStudy.graphic.title,
                    xaxis:{title: atomistic.constants.studies.currentStudy.graphic.titleX},
                    yaxis:{title: atomistic.constants.studies.currentStudy.graphic.titleY}
                };
            }
            else
            {
                var zValues = atomistic.constants.studies.savedPoints[atomistic.constants.studies.currentStudy.zAxis]
                var graph_content =
                {
                    x: xValues,
                    y: yValues,
                    z: zValues,
                    type: 'heatmap',
                    zsmooth: 'best',
                    colorscale: atomistic.constants.studies.currentStudy.graphic.palette
                };

                var graph_layout =
                {
                    title: atomistic.constants.studies.currentStudy.graphic.title,
                    xaxis:{title: atomistic.constants.studies.currentStudy.graphic.titleX},
                    yaxis:{title: atomistic.constants.studies.currentStudy.graphic.titleY},
                    colorbar:{title: atomistic.constants.studies.currentStudy.graphic.titleZ}
                };
            }
            
            
            Plotly.newPlot('graph_Q_B_T', [graph_content], graph_layout);
        }
    }
};