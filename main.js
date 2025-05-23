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
let InitSPH, Integrate, ComputeForces, ComputeDensityPressure;
let init_serialized_coordinate, get_serialized_coordinate_buffer_ptr, get_serialized_coordinate_buffer_size;
let get_num_particles;
let set_num_particles;

// 
let stats;
let balls = [];

//const use_instanced_mesh = false;
const use_instanced_mesh = true;

const run_controller = {
	num_particles: 2000,
	initialize_flag: false,
	start: function() {
		this.initialize_flag = true;
	},
	start_done: function() {
		this.initialize_flag = false;
	}
};

createModule().then((Module) => {
	InitSPH = Module.InitSPH;
	Integrate = Module.Integrate;
	ComputeForces = Module.ComputeForces;
	ComputeDensityPressure = Module.ComputeDensityPressure;
	init_serialized_coordinate = Module.init_serialized_coordinate;
	get_serialized_coordinate_buffer_ptr = Module.get_serialized_coordinate_buffer_ptr;
	get_serialized_coordinate_buffer_size = Module.get_serialized_coordinate_buffer_size;
	get_num_particles = Module.get_num_particles;
	set_num_particles = Module.set_num_particles;

	set_num_particles(run_controller.num_particles);
	InitSPH();
	console.log("InitSPH done!");
	init_serialized_coordinate();
	coordinateArray = new Float32Array(Module.HEAPF32.buffer, get_serialized_coordinate_buffer_ptr(), get_serialized_coordinate_buffer_size());
	console.log(coordinateArray);
	console.log("size: " + get_serialized_coordinate_buffer_size());
	console.log(coordinateArray.length);
	console.log(get_serialized_coordinate_buffer_ptr());

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
		gui.add(run_controller, 'num_particles').name("Number of Particles");
		gui.add(run_controller, 'start').name('Start');
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
			x /= 100;
			y /= 100;
			z /= 100;
			x -= 5;
			//z -= 5;
			balls[i].position.set(x, y, z);
		}
	}
	function setup_instanced_mesh() {
		if (mesh instanceof THREE.InstancedMesh) {
			console.log('dispose');
			mesh.dispose();
		}
		let n = get_num_particles();
		let ball_geometry = new THREE.SphereGeometry(0.05, 8, 8);
		let ball_material = new THREE.MeshBasicMaterial({color:0x0000ff});
		mesh = new THREE.InstancedMesh(ball_geometry, ball_material, n);
		scene.add(mesh);
	}
	function update_instanced_mesh_positions() {
		const dummy = new THREE.Object3D();
		console.log('instanced mesh');
		let n = get_num_particles();
		for(let i = 0; i < n; i++) {
			let x = coordinateArray[i*3+0];
			let y = coordinateArray[i*3+1];
			let z = coordinateArray[i*3+2];
			x /= 100;
			y /= 100;
			z /= 100;
			x -= 5;
			z -= 5;
			dummy.position.set(x, y, z)
			dummy.updateMatrix();
			mesh.setMatrixAt(i, dummy.matrix)
		}
		mesh.instanceMatrix.needsUpdate = true;
	}

	setup_three();
	setup_stats();
	setup_helper();
	if (use_instanced_mesh == true) {
		setup_instanced_mesh();
		update_instanced_mesh_positions();
	} else {
		setup_balls();
		update_ball_positions();
	}
	setup_gui();
	//balls = []

	animate();

	function animate() {
		if (run_controller.initialize_flag == true) {
			set_num_particles(run_controller.num_particles);
			InitSPH();
			run_controller.start_done();
			if (use_instanced_mesh == true) {
				setup_instanced_mesh();
			} else {
				setup_balls();
			}
		}
		for(let i = 0; i < 10; i++) {
			ComputeDensityPressure();
			ComputeForces();
			Integrate();
		}
		init_serialized_coordinate();
		console.log('calculation done');
		if (use_instanced_mesh == true) {
			update_instanced_mesh_positions();
	 	} else {
			update_ball_positions();
		}
	requestAnimationFrame(animate);
	cube.rotation.x += 0.01;
	cube.rotation.y += 0.01;
	renderer.render(scene, camera);
	stats.update();
	}
});

//animate();

