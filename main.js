import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import Stats from 'three/addons/libs/stats.module.js';
import GUI from 'lil-gui';
import createModule from './em.js';

let scene, camera, renderer;
let geometry, cube, material; 
let coordinateArray;
let xyz_geometry;
let mesh;	// instanced mesh

// Emscripten Exported Function

// 
let stats;
let balls = [];
let sim;	// simulation object;

//const use_instanced_mesh = false;
const use_instanced_mesh = true;

const run_controller = {
	num_particles: 2000,
	initialize_flag: false,
	run_flag: false,
	stop_soon: false,
	n_steps: 1,
	x_scale: 20.0,
	y_scale: 20.0,
	z_scale: 20.0,
	x_shift: 0, 
	y_shift: 0,
	z_shift: 0,
	sphere_scale: 1.0,
	num_forward: 1,
	num_forward_actual: 1,
	start: function() {
		this.initialize_flag = true;
	},
	start_done: function() {
		this.initialize_flag = false;
	},
	step_forward: function () {
		this.num_forward_actual = this.num_forward;
		this.run_flag = true;
	},
	stop: function() {
		this.run_flag = false;
	}
};

const sph_params = {
	//x_limit: 0.2,
	//y_limit: 0.2,
	//z_limit: 0.2,
	//rho0: 1000.0, //[kg/m3]
	//h: 0.026, // m	// stiffness
	//mass: 8.0e-3,	// [kg]
	//stiffness: 10000,	// [Pa]
	//dt: 0.001,	// [s]
	//dx: 0.02,	// [m]
	//visc: 300.0,	// [m]
	//epsilon: 0.01,
	x_limit: 1.0,
	y_limit: 2.0,
	z_limit: 1.0,
	rho0: 15000.0, //[kg/m3]
	h: 0.07, // m	// stiffness
	mass: 1.0,	// [kg]
	stiffness: 20,	// [Pa]
	dt: 0.006,	// [s]
	dx: 0.04,	// [m]
	visc: 100.0,	// [m]
	epsilon: 0.01,
};

createModule().then((Module) => {
	function setup_three() {
		console.log("setup Scene");
		scene = new THREE.Scene();
		camera = new THREE.PerspectiveCamera(
		75,
		window.innerWidth / window.innerHeight,
		0.1,
		1000
		);
		camera.position.x = 5;
		camera.position.y = 5;
		camera.position.z = 5;
		camera.lookAt(15, 0);

		renderer = new THREE.WebGLRenderer();
		renderer.setSize(window.innerWidth, window.innerHeight);
		document.body.appendChild(renderer.domElement);

		geometry = new THREE.BoxGeometry();
		material = new THREE.MeshNormalMaterial();
		cube = new THREE.Mesh(geometry, material);
		//scene.add(cube);

		const controls = new OrbitControls(camera, renderer.domElement);
	}

	function setup_stats() {
		console.log('Stats');
		stats = new Stats();
		stats.showPanel(0);
		document.body.appendChild(stats.dom);
	}

	function setup_gui() {
		const gui = new GUI();
		const folder = gui.addFolder('view settings');
		folder.add(run_controller, 'x_scale').name('x scale factor');
		folder.add(run_controller, 'y_scale').name('y scale factor');
		folder.add(run_controller, 'z_scale').name('z scale factor');
		folder.add(run_controller, 'x_shift').name('x shift after scaling');
		folder.add(run_controller, 'y_shift').name('y shift after scaling');
		folder.add(run_controller, 'z_shift').name('z shift after scaling');

		gui.add(run_controller, 'num_particles').name("Number of Particles");
		gui.add(run_controller, 'start').name('Start');
		gui.add(run_controller, 'step_forward').name('Step Forward');
		gui.add(run_controller, 'stop').name('Stop');
		gui.add(run_controller, 'num_forward').name('Number of steps');
		//gui.add(run_controller, "run_flag").name('Run');
		gui.add(run_controller, 'n_steps').name('n_steps');

		const folder_sph_params = gui.addFolder('SPH Parameters');
		folder_sph_params.add(sph_params, 'x_limit').name('x_limit [m]').onChange(value => {
			console.log('x_limit');
			//sim.x_limit = value;
		});
		folder_sph_params.add(sph_params, 'y_limit').name('y_limit [m]').onChange(value => {
			console.log('y_limit');
			//sim.y_limit = value;
		});
		folder_sph_params.add(sph_params, 'z_limit').name('z_limit [m]');
		folder_sph_params.add(sph_params, 'dx').name('dx [m]');
		folder_sph_params.add(sph_params, 'dt').name('dt [s]');
		folder_sph_params.add(sph_params, 'rho0').name('rho0 [kg/m^3]');
		folder_sph_params.add(sph_params, 'h').name('h [m]');
		folder_sph_params.add(sph_params, 'mass').name('mass [kg]');
		folder_sph_params.add(sph_params, 'stiffness').name('stiffness [Pa]');
		folder_sph_params.add(sph_params, 'visc').name('viscocity constant');

	}

	function setup_helper() {
		const axes= new THREE.AxesHelper();
		const gridhelper = new THREE.GridHelper(10, 10);
		scene.add(axes);
		scene.add(gridhelper);
	}

	function setup_balls() {
		if (0 < balls.length) {
			balls.splice(0);
		}
		for(let i = 0; i < n; i++) {
			let ball_geometry = new THREE.SphereGeometry(0.05, 8, 8);
			let ball_material = new THREE.MeshBasicMaterial({color:0x0000ff});
			let ball= new THREE.Mesh(ball_geometry, ball_material)
			balls.push(ball);
			scene.add(ball);
		}
	}
	function update_ball_positions() {
		let n = get_num_particles();
		for(let i = 0; i < n; i++) {
			let x = coordinateArray[i*3+0]; //300;
			let y = coordinateArray[i*3+1]; //300;
			let z = coordinateArray[i*3+2]; //300;
			//console.log(x)
			x *= run_controller.x_scale;
			y *= run_controller.y_scale;
			z *= run_controller.z_scale;
			x += run_controller.x_shift;
			y += run_controller.y_shift;
			z += run_controller.z_shift;
			//z -= 5;
			balls[i].position.set(x, y, z);
		}
	}
	function setup_instanced_mesh() {
		if (mesh instanceof THREE.InstancedMesh) {
			console.log('dispose');
			mesh.dispose();
		}
		let n = sim.get_num_particles();
		let ball_geometry = new THREE.SphereGeometry(0.05, 8, 8);
		let ball_material = new THREE.MeshBasicMaterial({color:0x0000ff});
		mesh = new THREE.InstancedMesh(ball_geometry, ball_material, n);
		console.log(sim.x_limit, sim.y_limit, sim.z_limit);
		scene.add(mesh);
	}
	function update_instanced_mesh_positions() {
		const dummy = new THREE.Object3D();
		let n = sim.get_num_particles();
		console.log(n);
		for(let i = 0; i < n; i++) {
			let x = coordinateArray[i*3+0];
			let y = coordinateArray[i*3+1];
			let z = coordinateArray[i*3+2];
			x *= run_controller.x_scale;
			y *= run_controller.y_scale;
			z *= run_controller.z_scale;
			x += run_controller.x_shift;
			y += run_controller.y_shift;
			z += run_controller.z_shift;
			//z -= 5;
			dummy.position.set(x, y, z)
			dummy.updateMatrix();
			mesh.setMatrixAt(i, dummy.matrix)
		}
		mesh.instanceMatrix.needsUpdate = true;
	}

	function setup_simulation() {
		if (sim instanceof Module.SPHSimulator) { sim.delete(); }
		sim = new Module.SPHSimulator(run_controller.num_particles);

		sim.x_limit = sph_params.x_limit;
		sim.y_limit = sph_params.y_limit;
		sim.z_limit = sph_params.z_limit;
		sim.rho0 = sph_params.rho0;
		sim.h = sph_params.h;
		sim.mass = sph_params.mass;
		sim.stiffness = sph_params.stiffness;
		sim.dt = sph_params.dt;
		sim.dx = sph_params.dx;
		sim.epsilon = sph_params.epsilon;
		sim.visc = sph_params.visc;

		sim.init_simulation(true);
		console.log(sim.get_num_particles());
		coordinateArray = new Float32Array(
			Module.HEAPF32.buffer, 
			sim.get_serialized_coordinate_buffer_ptr(), 
			sim.get_serialized_coordinate_buffer_size()
		);
	}

	setup_three();
	setup_stats();
	setup_helper();
	setup_gui();
	setup_simulation();
	//const box_geometry = new THREE.BoxGeometry(sim.x_limit * run_controller.x_scale, sim.y_limit * run_controller.y_scale, sim.z_limit * run_controller.z_scale);
	//const box_geometry = new THREE.BoxGeometry(sim.x_limit, sim.y_limit, sim.z_limit);
	const box_geometry = new THREE.BoxGeometry(1.0, 1.0, 1.0);
	//console.log(sim.x_limit, sim.y_limit, sim.z_limit);
	const box_edges = new THREE.EdgesGeometry(box_geometry);
	// 線として描画（白線など）
	const box_material = new THREE.LineBasicMaterial({ color: 0xffffff });
	const wireframe = new THREE.LineSegments(box_edges, box_material);
	scene.add(wireframe);
	if (use_instanced_mesh == true) {
		setup_instanced_mesh();
		update_instanced_mesh_positions();
	} else {
		setup_balls();
		update_ball_positions();
	}
	//balls = []

	animate();

	function animate() {
		wireframe.scale.set(
			sim.x_limit * run_controller.x_scale, 
			sim.y_limit * run_controller.y_scale, 
			sim.z_limit * run_controller.z_scale);
		wireframe.position.set(
			sim.x_limit * run_controller.x_scale / 2.0, 
			sim.y_limit * run_controller.y_scale / 2.0, 
			sim.z_limit * run_controller.z_scale / 2.0 );
		if (run_controller.initialize_flag == true) {
			setup_simulation();
			run_controller.start_done();
			if (use_instanced_mesh == true) {
				setup_instanced_mesh();
			} else {
				setup_balls();
			}
		}
		if (run_controller.run_flag == true) {
			if (0 == run_controller.num_forward_actual) {
				run_controller.run_flag = false;
				console.log('stop');
			} else if (0 < run_controller.num_forward_actual) {
				sim.step(run_controller.n_steps, true);
				run_controller.num_forward_actual -= 1;
				console.log('Rest');
				console.log(run_controller.num_forward_actual);
			} else if (-1 == run_controller.num_forward_actual) {
				sim.step(run_controller.n_steps, true);
				console.log('inf');
			}
		}

		if (use_instanced_mesh == true) {
			update_instanced_mesh_positions();
		} else {
			update_ball_positions();
		}
		if (run_controller.stop_soon == true) {
			run_controller.run_flag = false;
		}
	requestAnimationFrame(animate);
	cube.rotation.x += 0.01;
	cube.rotation.y += 0.01;
	renderer.render(scene, camera);
	stats.update();
	}
});

//animate();

