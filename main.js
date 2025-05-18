import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import createModule from './em.js';

let scene, camera, renderer;
let geometry, cube, material; 
let coordinateArray;
let xyz_geometry;

let InitSPH, Integrate, ComputeForces, ComputeDensityPressure;
let init_serialized_coordinate, get_serialized_coordinate_buffer_ptr, get_serialized_coordinate_buffer_size;
let num_particles;
let set_num_particles;
const n = 2000;
let balls = [];
createModule().then((Module) => {
	InitSPH = Module.InitSPH;
	Integrate = Module.Integrate;
	ComputeForces = Module.ComputeForces;
	ComputeDensityPressure = Module.ComputeDensityPressure;
	init_serialized_coordinate = Module.init_serialized_coordinate;
	get_serialized_coordinate_buffer_ptr = Module.get_serialized_coordinate_buffer_ptr;
	get_serialized_coordinate_buffer_size = Module.get_serialized_coordinate_buffer_size;
	num_particles = Module.num_particles;
	set_num_particles = Module.set_num_particles;

	set_num_particles(n);
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

	function setup_helper() {
		const axes= new THREE.AxesHelper();
		const gridhelper = new THREE.GridHelper(10, 10);
		scene.add(axes);
		scene.add(gridhelper);
	}


	setup_three();
	setup_helper();
	//balls = []
	for(let i = 0; i < n; i++) {
		let ball_geometry = new THREE.SphereGeometry(0.05, 8, 8);
		let ball_material = new THREE.MeshBasicMaterial({color:0x0000ff});
		let ball= new THREE.Mesh(ball_geometry, ball_material)
		let x = coordinateArray[i*3+0]; //300;
		let y = coordinateArray[i*3+1]; //300;
		let z = coordinateArray[i*3+2]; //300;
		x /= 100;
		y /= 100;
		z /= 100;
		x -= 5;
		ball.position.set(x, y, z);
		balls.push(ball);
		scene.add(ball);
	}



	animate();

	function animate() {
		for(let i = 0; i < 10; i++) {
			ComputeDensityPressure();
			ComputeForces();
			Integrate();
		}
		init_serialized_coordinate();
		console.log('calculation done');
		for(let i = 0; i < n; i++) {
			let x = coordinateArray[i*3+0]; //300;
			let y = coordinateArray[i*3+1]; //300;
			let z = coordinateArray[i*3+2]; //300;
			//console.log(x)
			x /= 100;
			y /= 100;
			z /= 100;
			x -= 5;
			balls[i].position.set(x, y, z);
			//balls[i].position.set(x, y, z);
			//ball.position.set(0 + i, 0, 0);
		}
	requestAnimationFrame(animate);
	cube.rotation.x += 0.01;
	cube.rotation.y += 0.01;
	renderer.render(scene, camera);

	}
});

//animate();

