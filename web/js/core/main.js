"use strict";

// Méthode globale
var atomistic = 
{
    // Constantes fixes et calculées
    constants:
    {
        // Paramètre de simulation
        simulation:
        {
            counter: 0, // Nombre de pas effectué
            debug: false, // Mode "debug" activé ou non
            struct: null, // Structure vectorielle de l'étude
            results: [], // Résultats de la simulation
            nbSteps: Infinity, // Nombre de pas maximal avant la fin de la simulation (Infinity : pas de limite)
            pause: false, // Simulation en pause ou non
            updateEveryStep: 1, // Taux de rafraîchissement des paramètres calculés
            neighbors: "strict", // Construction de la liste des voisins
        },
        system:
        {
            dx: 15,
            dy: 15,
            dz: 1,
            dim: 3
        },
        physicalParameters:
        {
            α: 1.0,
            τ: 1e-11,
            dt: 1e-14,
            D: 0.0,
            J: 10e-3,
            K: 10e-3,
            L: 0.0,
            B: 0,
        },
        observables:
        {
            s: null,
            W: null,
            S: null,
            h: null,
            listNeighbors: null,
            Q: null,
            Q_sum: 0,
            sMean: null,
            sNorm: null,
            WMean: null,
            SMean: null,
        },
        multithreading:
        {
            currentPool: null,
            nbWorkers: 2,
            workers: null,
            worker_energy: null,
            worker_neighbors: null
        },
        renderer:
        {
            scene: null,
            camera: null,
            renderer: null,
            palette: null,
            live: false,
            quality: 1.5,
            nuances: 1000,
            currentTab: "3d"
        },
        studies:
        {
            study: true,
            initialized: false,
            dim: "3d",
            savedPoints:
            {
                "Q_sum": [],
                "B": [],
                "D": []
            },
            currentStudy:
            {
                x: 0,
                dx: 5e14,
                xmax: 5e15,
                y: 0,
                dy: 0.5,
                ymax: 10,
                delta: 1,
                maxSteps: 2000,
                maxDelta: 1e-5,
                xAxis: "D",
                yAxis: "B",
                zAxis: "Q_sum",
                counterX: 0,
                counterY: 0,
                graphic:
                {
                    title: "Diagramme de stabilité des skyrmions",
                    titleX: "Bruit thermique D (rad.Hz)",
                    titleY: "Intensité du champ magnétique (T)",
                    titleZ: "Charge topologique Q",
                    palette: "Jet"
                }
            },
            initialConfig:
            {
                s: null,
                S: null
            },
        }
    }
} || {};