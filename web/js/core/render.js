var create_palette = function()
{
    atomistic.constants.renderer.palette = [];

    for(i = 0; i < atomistic.constants.renderer.nuances; i++)
    {
        atomistic.constants.renderer.palette.push([Math.round(i * 359 / (atomistic.constants.renderer.nuances-1)), 100, 25]);
    }
}

var create_grid = function(dx, dy, dz)
{
    create_palette();
    
    var scene = new THREE.Scene();
    var renderer = new THREE.WebGLRenderer({antialias: true});
    renderer.setSize( document.querySelector("section").clientWidth, document.querySelector("section").clientHeight );
    document.querySelector("#sectionTab_3d").innerHTML = "";
    document.querySelector("#sectionTab_3d").appendChild( renderer.domElement );
    
    var geometry = new THREE.CylinderGeometry(0.0, 0.4, 1, 10, 0, false);
    var material	= new THREE.MeshPhongMaterial({
		ambient		: 0x444444,
		color		: 0xFFFFFF,
		shininess	: 300, 
		specular	: 0x33AA33,
		shading		: THREE.SmoothShading
	});
    var geometry2 = new THREE.CylinderGeometry(0.2, 0.2, 1.5, 10, 0, false);
    var singleGeometry = new THREE.Geometry();
    var vector = new THREE.Vector3(0, 0, 1);
    var axis = new THREE.Vector3(0, 1, 0);

    for(i = 0; i < dx; i++)
    {
        for(j = 0; j < dy; j++)
        {
            for(k = 0; k < dz; k++)
            {
                var cone = new THREE.Mesh( geometry);
                var cylinder = new THREE.Mesh(geometry2);

                cylinder.position.set(0, -0.5, 0);

                cylinder.updateMatrix();
                singleGeometry.merge(cylinder.geometry, cylinder.matrix);

                cone.position.set(0, 0.75, 0);
                cone.updateMatrix();
                singleGeometry.merge(cone.geometry, cone.matrix);

                var mesh = new THREE.Mesh(singleGeometry, material);


                mesh.position.set(i * 3 - dx * 3 / 2 + 1.5, j * 3 - dy * 3 / 2 + 1.5, k * 3);
                mesh.name = "spin_" + i + "_" + j + "_" + k;
                mesh.quaternion.setFromUnitVectors(axis, vector.clone().normalize());

                scene.add(mesh);
            }
        }
    }

    var ambient	= new THREE.AmbientLight( 0x444444 );
	scene.add( ambient );

    var spotLight = new THREE.SpotLight( 0xFFAA88);
    spotLight.target.position.set( 0, 0, 0 );
    spotLight.position.set( 0, 0, 20 );
    spotLight.shadowCameraNear	= 0.001;

    scene.add( spotLight );

    var axisHelper = new THREE.AxisHelper( 10 );
    scene.add( axisHelper );

    var camera = new THREE.PerspectiveCamera( 45, document.querySelector("section").clientWidth / document.querySelector("section").clientHeight, 0.1, 1000000 );
    camera.name = "camera";
    scene.add(camera);

    document.onkeydown = function(e)
    {
        update_camera(e);
    }

    document.onwheel = function(e)
    {
        update_camera(e);
    }
    
    var camera_pivot = new THREE.Object3D();
    camera_pivot.name = "camera_pivot";
    var AXIS_X = new THREE.Vector3(1, 0, 0);
    scene.add(camera_pivot);
    
    camera_pivot.add(camera);
    camera.position.set(0, 0, 3.2 * dx).length(3.2 * dx);
    camera.lookAt(camera_pivot.position);
    camera_pivot.rotateOnAxis(AXIS_X, 40 * Math.PI / 100);

    //camera.lookAt(new THREE.Vector3(0, 0, 0));

    atomistic.constants.renderer.scene = scene;
    atomistic.constants.renderer.camera = camera;

    renderer.gammaOutput = true;
    renderer.setPixelRatio( window.devicePixelRatio * atomistic.constants.renderer.quality );

    atomistic.constants.renderer.renderer = renderer;

    renderer.render( scene, camera );
}

var update_render = function(s, scene, camera, dx, dy, dz)
{    
    for(i = 0; i < dx; i++)
    {
        for(j = 0; j < dy; j++)
        {
            for(k = 0; k < dz; k++)
            {
                var object = scene.getObjectByName("spin_" + i + "_" + j + "_" + k);

                var vector = new THREE.Vector3(s[i][j][k][0], s[i][j][k][1], s[i][j][k][2]);
                var axis = new THREE.Vector3(0, 1, 0);
                object.quaternion.setFromUnitVectors(axis, vector.clone().normalize());

                var theta = Math.atan2(s[i][j][k][1], s[i][j][k][0]);

                var coeff = Math.round((theta + Math.PI) * (atomistic.constants.renderer.nuances-1) / (2 * Math.PI));
                
                var coeff_2 = Math.round((s[i][j][k][2] + 1.05) * 100 / 2);

                object.material = new THREE.MeshPhongMaterial({
                    ambient		: 0x444444,
                    color		: "hsl("+atomistic.constants.renderer.palette[coeff][0]+", "+atomistic.constants.renderer.palette[coeff][1]+"%, "+coeff_2+"%)",
                    shininess	: 300, 
                    specular	: 0x33AA33,
                    shading		: THREE.SmoothShading
                });
            }
        }
    }

    atomistic.constants.renderer.renderer.render(scene, camera);
}

var update_camera = function(e)
{
    console.log(e);
    var AXIS_X = new THREE.Vector3(1, 0, 0);
    var AXIS_Y = new THREE.Vector3(0, 1, 0);
    var AXIS_Z = new THREE.Vector3(0, 0, 1);

    switch(e.key)
    {
        case "ArrowRight":
            e.preventDefault();
            atomistic.constants.renderer.scene.getObjectByName("camera_pivot").rotation.y += 0.05;

            atomistic.constants.renderer.renderer.render( atomistic.constants.renderer.scene, atomistic.constants.renderer.camera );
            break;

        case "ArrowLeft":
            e.preventDefault();
            atomistic.constants.renderer.scene.getObjectByName("camera_pivot").rotation.y -= 0.05;

            atomistic.constants.renderer.renderer.render( atomistic.constants.renderer.scene, atomistic.constants.renderer.camera );
            break;

        case "ArrowDown":
            e.preventDefault();
            atomistic.constants.renderer.scene.getObjectByName("camera_pivot").rotation.x -= 0.05;

            atomistic.constants.renderer.renderer.render( atomistic.constants.renderer.scene, atomistic.constants.renderer.camera );
            break;

        case "ArrowUp":
            e.preventDefault();
            atomistic.constants.renderer.scene.getObjectByName("camera_pivot").rotation.x += 0.05;

            atomistic.constants.renderer.renderer.render( atomistic.constants.renderer.scene, atomistic.constants.renderer.camera );
            break;
    }

    switch(e.type)
    {
        case "wheel":
            e.preventDefault();

            atomistic.constants.renderer.scene.getObjectByName("camera").position.z += e.deltaY * -1;

            atomistic.constants.renderer.renderer.render( atomistic.constants.renderer.scene, atomistic.constants.renderer.camera );
            break;
    }
}