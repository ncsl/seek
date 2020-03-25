
import { OrbitControls } from './OrbitControls.js';
import { GLTFLoader } from './GLTFLoader.js'
import * as THREE from './three.module.js'

let scene;
let gyri;
let wm;
let subStructures;
let gyriScale = 1;
let wmScale = 1;
let subStructuresScale = 1;
window.onload = async () => {
	scene = await load3DBrain_gltf();
	wm = scene.getObjectByName("WhiteMatter");
	gyri = scene.getObjectByName("Gyri");
	electrodes = scene.getObjectByName("Electrodes");
	subStructures = scene.getObjectByName("Brain");
}


let load3DBrain_gltf = () => {
	return new Promise((resolve, reject) => {
		let brainContainer = document.getElementById('fm-brain-3D');
		let scene = new THREE.Scene();
		scene.background = new THREE.Color(0xffffff);
		let camera = new THREE.PerspectiveCamera(45, 640 / 480, 0.1, 5000);
		let renderer = new THREE.WebGLRenderer({
			antialias: true
		});
		let controls = new OrbitControls(camera, renderer.domElement);
		let light = new THREE.HemisphereLight(0xffffff, 0x444444);
		camera.position.set(-100, 0, -100);
		renderer.setSize(640, 480);
		light.position.set(0, 0, 10);
		controls.target.set(10, 20, 0);
		controls.update();
		scene.add(light);
		Cache.enabled = true;
		let loader = new GLTFLoader()
		loader.load('/brain', object3d => {
			scene.add(object3d.scene)
			let mainScene = scene.getObjectByName("Scene");
			mainScene.rotation.set(-Math.PI / 2, 0, 0)
			resolve(scene)
		})
		const animate = () => {
			requestAnimationFrame(animate);
			renderer.render(scene, camera);
		};
		animate();
		brainContainer.appendChild(renderer.domElement);
	})
}
document.getElementById('transparencyToggle_g').onclick = () => {
	gyri.visible = true;
	wm.visible = false;
	if (gyriScale < .1) {
		gyriScale = 1;
	}
	else {
		gyriScale = gyriScale - .1
	}

	gyri.traverse(child => {
		if (child instanceof THREE.Mesh) {
			if (child.material.opacity < .0) {
				child.material.transparent = false;
				child.material.opacity = gyriScale;
			}
			else {
				child.material.transparent = true;
				child.material.opacity = gyriScale;
			}
		}
	});
}
document.getElementById('transparencyToggle_WM').onclick = () => {
	wm.visible = true;
	gyri.visible = false;
	if (wmScale < .1) {
		wmScale = 1;
	}
	else {
		wmScale = wmScale - .1
	}
	wm.traverse(child => {
		if (child instanceof THREE.Mesh) {
			if (child.material.opacity < .0) {
				child.material.transparent = false;
				child.material.opacity = wmScale;
			}
			else {
				child.material.transparent = true;
				child.material.opacity = wmScale;
			}
		}
	});
}
document.getElementById('transparencyToggle_sub').onclick = () => {
	subStructures.visible = true;
	gyri.visible = false;
	wm.visible = false;
	if (subStructuresScale < .1) {
		subStructuresScale = 1;
	}
	else {
		subStructuresScale = subStructuresScale - .1
	}
	subStructures.traverse(child => {
		if (child instanceof THREE.Mesh) {
			if (child.material.opacity < .0) {
				child.material.transparent = false;
				child.material.opacity = subStructuresScale;
			}
			else {
				child.material.transparent = true;
				child.material.opacity = subStructuresScale;
			}
		}
	});
}

// document.getElementById('electrodeScaler').onclick = () => {
// 	let geometry, line;
// 	let stimPosition, stimElectrode;
// 	let empty;
// 	let material
// 	electrodes.traverse(child => {
// 		electrodeNames.forEach(elec => {
// 			if (elec == child.name) {
// 				stimElectrode = scene.getObjectByName("LAA1")
// 				geometry = new THREE.Geometry();
// 				geometry.vertices.push(new THREE.Vector3(0, 0, 0));
// 				let x = child.getWorldPosition().x - stimElectrode.getWorldPosition().x
// 				let y = child.getWorldPosition().y - stimElectrode.getWorldPosition().y
// 				let z = child.getWorldPosition().z - stimElectrode.getWorldPosition().z
// 				geometry.vertices.push(new THREE.Vector3(x, y, z));
// 				material = new THREE.LineBasicMaterial({ color: 0x0000ff });
// 				material.color.setHSL(resInfo[elec] / Math.max(...elecValues), 1, .5)
// 				line = new THREE.Line(geometry, material);
// 				scene.add(line);
// 				stimElectrode.getWorldPosition(line.position)
// 			}
// 		})

// 		// geometry = new THREE.BoxBufferGeometry(3, 3, 3);
// 		// material = new THREE.MeshBasicMaterial({ color: 0xff0000 });
// 		// mesh = new THREE.Mesh(geometry, material);
// 		// mesh.position.copy(child.position);
// 		if (child instanceof THREE.Mesh) {
// 			if (child.scale.x > 2) {
// 				child.scale.set(1, 1, 1)
// 			}
// 			else {
// 				child.scale.set(child.scale.x + .1, child.scale.y + .1, child.scale.z + .1)
// 			}
// 		}
// 	})
// }
