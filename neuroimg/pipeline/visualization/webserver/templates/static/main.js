
const StackHelper = AMI.stackHelperFactory(THREE);
const CamerasOrthographic = AMI.orthographicCameraFactory(THREE);
const ControlsOrthographic = AMI.trackballOrthoControlFactory(THREE);
const ControlsTrackball = AMI.trackballControlFactory(THREE);
const HelpersLut = AMI.lutHelperFactory(THREE);

// import { GLTFLoader } from './GLTFLoader.js'


let lut;
let ready = false;
let elecs;
let brainScene, wm, gyri, substructures;

const r0 = {
	domId: 'r0',
	domElement: null,
	renderer: null,
	color: 0x212121,
	targetID: 0,
	camera: null,
	controls: null,
	scene: null,
	light: null,
};

// 2d axial renderer
const r1 = {
	domId: 'r1',
	domElement: null,
	renderer: null,
	color: 0x121212,
	sliceOrientation: 'axial',
	sliceColor: 0xff1744,
	targetID: 1,
	camera: null,
	controls: null,
	scene: null,
	light: null,
	stackHelper: null,
};

// 2d sagittal renderer
const r2 = {
	domId: 'r2',
	domElement: null,
	renderer: null,
	color: 0x121212,
	sliceOrientation: 'sagittal',
	sliceColor: 0xffea00,
	targetID: 2,
	camera: null,
	controls: null,
	scene: null,
	light: null,
	stackHelper: null,
};

// 2d coronal renderer
const r3 = {
	domId: 'r3',
	domElement: null,
	renderer: null,
	color: 0x121212,
	sliceOrientation: 'coronal',
	sliceColor: 0x76ff03,
	targetID: 3,
	camera: null,
	controls: null,
	scene: null,
	light: null,
	stackHelper: null,
};


let data = []

let sceneClip = new THREE.Scene();
let clipPlane1 = new THREE.Plane(new THREE.Vector3(0, 0, 0), 0);
let clipPlane2 = new THREE.Plane(new THREE.Vector3(0, 0, 0), 0);
let clipPlane3 = new THREE.Plane(new THREE.Vector3(0, 0, 0), 0);

function initRenderer3D(renderObj) {
	// renderer
	renderObj.domElement = document.getElementById(renderObj.domId);
	renderObj.renderer = new THREE.WebGLRenderer({
		antialias: true,
	});

	renderObj.renderer.setSize(renderObj.domElement.clientWidth, renderObj.domElement.clientHeight);
	renderObj.renderer.setClearColor(renderObj.color, 1);
	renderObj.renderer.domElement.id = renderObj.targetID;
	renderObj.domElement.appendChild(renderObj.renderer.domElement);

	// camera
	renderObj.camera = new THREE.PerspectiveCamera(
		45,
		renderObj.domElement.clientWidth / renderObj.domElement.clientHeight,
		0.1,
		100000
	);
	renderObj.camera.position.x = 250;
	renderObj.camera.position.y = 250;
	renderObj.camera.position.z = 250;

	// controls
	renderObj.controls = new ControlsTrackball(renderObj.camera, renderObj.domElement);
	renderObj.controls.rotateSpeed = 5.5;
	renderObj.controls.zoomSpeed = 1.2;
	renderObj.controls.panSpeed = 0.8;
	renderObj.controls.staticMoving = true;
	renderObj.controls.dynamicDampingFactor = 0.3;

	renderObj.scene = new THREE.Scene();

	renderObj.light = new THREE.DirectionalLight(0xffffff, 1);
	renderObj.light.position.copy(renderObj.camera.position);
	renderObj.scene.add(renderObj.light);
}

function initRenderer2D(rendererObj) {
	// renderer
	rendererObj.domElement = document.getElementById(rendererObj.domId);
	rendererObj.renderer = new THREE.WebGLRenderer({
		antialias: true,
	});
	rendererObj.renderer.autoClear = false;
	rendererObj.renderer.localClippingEnabled = true;
	rendererObj.renderer.setSize(
		rendererObj.domElement.clientWidth,
		rendererObj.domElement.clientHeight
	);
	rendererObj.renderer.setClearColor(0x121212, 1);
	rendererObj.renderer.domElement.id = rendererObj.targetID;
	rendererObj.domElement.appendChild(rendererObj.renderer.domElement);

	// camera
	rendererObj.camera = new CamerasOrthographic(
		rendererObj.domElement.clientWidth / -2,
		rendererObj.domElement.clientWidth / 2,
		rendererObj.domElement.clientHeight / 2,
		rendererObj.domElement.clientHeight / -2,
		1,
		1000
	);

	// controls
	rendererObj.controls = new ControlsOrthographic(rendererObj.camera, rendererObj.domElement);
	rendererObj.controls.staticMoving = true;
	rendererObj.controls.noRotate = true;
	rendererObj.camera.controls = rendererObj.controls;

	// scene
	rendererObj.scene = new THREE.Scene();
}

function initHelpersStack(rendererObj, stack) {
	rendererObj.stackHelper = new StackHelper(stack);
	rendererObj.lutLayer0 = new HelpersLut(
		'main',
		'default',
		'linear',
		[[0, 0, 0, 0], [1, 1, 1, 1]],
		[[0, 1], [1, 1]]
	);
	rendererObj.lutLayer0.luts = HelpersLut.presetLuts();
	rendererObj.stackHelper.bbox.visible = false;
	rendererObj.stackHelper.borderColor = rendererObj.sliceColor;
	rendererObj.stackHelper._slice.borderColor = 0x000000;
	rendererObj.stackHelper.slice.canvasWidth = rendererObj.domElement.clientWidth;
	rendererObj.stackHelper.slice.canvasHeight = rendererObj.domElement.clientHeight;
	// console.log(rendererObj.lutLayer0)
	rendererObj.stackHelper.slice.lutTexture = rendererObj.lutLayer0.texture

	// set camera
	let worldbb = stack.worldBoundingBox();
	let lpsDims = new THREE.Vector3(
		(worldbb[1] - worldbb[0]) / 2,
		(worldbb[3] - worldbb[2]) / 2,
		(worldbb[5] - worldbb[4]) / 2
	);

	// box: {halfDimensions, center}
	let box = {
		center: stack.worldCenter().clone(),
		halfDimensions: new THREE.Vector3(lpsDims.x + 50, lpsDims.y + 50, lpsDims.z + 50),
	};

	// init and zoom
	let canvas = {
		width: rendererObj.domElement.clientWidth,
		height: rendererObj.domElement.clientHeight,
	};

	rendererObj.camera.directions = [stack.xCosine, stack.yCosine, stack.zCosine];
	rendererObj.camera.box = box;
	rendererObj.camera.canvas = canvas;
	rendererObj.camera.orientation = rendererObj.sliceOrientation;
	rendererObj.camera.update();
	rendererObj.camera.fitBox(2, 1);

	rendererObj.stackHelper.orientation = rendererObj.camera.stackOrientation;
	rendererObj.stackHelper.index = Math.floor(rendererObj.stackHelper.orientationMaxIndex / 2);
	rendererObj.scene.add(rendererObj.stackHelper);

}
function render() {
	// we are ready when both meshes have been loaded
	if (ready) {
		// render
		r0.controls.update();
		r1.controls.update();
		r2.controls.update();
		r3.controls.update();

		r0.light.position.copy(r0.camera.position);
		r0.renderer.render(r0.scene, r0.camera);

		r1.renderer.clear();
		r1.renderer.render(r1.scene, r1.camera);
		r1.renderer.clearDepth();
		data.forEach(object => {
			object.materialFront.clippingPlanes = [clipPlane1];
			object.materialBack.clippingPlanes = [clipPlane1];
		});
		r1.renderer.render(sceneClip, r1.camera);
		r1.renderer.clearDepth();

		r2.renderer.clear();
		r2.renderer.render(r2.scene, r2.camera);
		r2.renderer.clearDepth();
		data.forEach(object => {
			object.materialFront.clippingPlanes = [clipPlane2];
			object.materialBack.clippingPlanes = [clipPlane2];
		});
		r2.renderer.render(sceneClip, r2.camera);
		r2.renderer.clearDepth();

		r3.renderer.clear();
		r3.renderer.render(r3.scene, r3.camera);
		r3.renderer.clearDepth();
		data.forEach(object => {
			object.materialFront.clippingPlanes = [clipPlane3];
			object.materialBack.clippingPlanes = [clipPlane3];
		});
		r3.renderer.render(sceneClip, r3.camera);
		r3.renderer.clearDepth();
	}

}



function init() {
	function animate() {
		render();
		requestAnimationFrame(() => animate());
	}

	initRenderer3D(r0);
	initRenderer2D(r1);
	initRenderer2D(r2);
	initRenderer2D(r3);

	animate();
}



window.onload = async () => {

	init();
	let brainReq = await fetch(`/nifti`)
	let brain = await brainReq.blob()
	let objectURL = URL.createObjectURL(brain);
	console.log(brain)
	let loader = new AMI.VolumeLoader();
	loader.load('./static/brain.mgz')
		.then(() => {
			let series = loader.data[0].mergeSeries(loader.data)[0];
			loader.free();
			loader = null;
			let stack = series.stack[0];
			stack._frame.forEach(frame => {
				frame._imagePosition = [-128, 128, 128]
			})

			stack.prepare();
			// center 3d camera/control on the stack
			let centerLPS = stack.worldCenter();

			r0.camera.lookAt(centerLPS.x, centerLPS.y, centerLPS.z);
			r0.camera.updateProjectionMatrix();
			r0.controls.target.set(centerLPS.x, centerLPS.y, centerLPS.z);


			// red slice
			initHelpersStack(r1, stack);
			r0.scene.add(r1.scene);

			// yellow slice
			initHelpersStack(r2, stack);
			r0.scene.add(r2.scene);

			// green slice
			initHelpersStack(r3, stack);
			r0.scene.add(r3.scene);

			let gui = new dat.GUI({
				autoPlace: false,
			});

			let customContainer = document.getElementById('my-gui-container');
			customContainer.appendChild(gui.domElement);


			lut = new HelpersLut(
				'slices',
				'default',
				'linear',
				[[0, 0, 0, 0], [1, 1, 1, 1]],
				[[0, 1], [1, 1]]
			);
			lut.luts = HelpersLut.presetLuts();

			let LutFolder = gui.addFolder('LUT');
			let lutUpdate1 = LutFolder.add(r1.stackHelper.slice, 'lut', lut.lutsAvailable());
			let lutUpdate2 = LutFolder.add(r2.stackHelper.slice, 'lut', lut.lutsAvailable());
			let lutUpdate3 = LutFolder.add(r3.stackHelper.slice, 'lut', lut.lutsAvailable());
			lutUpdate1.onChange(value => {
				lut.lut = value;
				r1.stackHelper.slice.lutTexture = lut.texture;
			});
			lutUpdate2.onChange(value => {
				lut.lut = value;
				r2.stackHelper.slice.lutTexture = lut.texture;
			});
			lutUpdate3.onChange(value => {
				lut.lut = value;
				r3.stackHelper.slice.lutTexture = lut.texture;
			});

			// Red
			let stackFolder = gui.addFolder('Slicer');
			let redChanged = stackFolder
				.add(r1.stackHelper, 'index', 0, r1.stackHelper.orientationMaxIndex)
				.step(1)
				.listen();
			let yellowChanged = stackFolder
				.add(r2.stackHelper, 'index', 0, r2.stackHelper.orientationMaxIndex)
				.step(1)
				.listen();
			let greenChanged = stackFolder
				.add(r3.stackHelper, 'index', 0, r3.stackHelper.orientationMaxIndex)
				.step(1)
				.listen();

			var meshController = function () {
				this.Mesh = true
				this.Cortex = true
				this.WM = true
				this.Substructures = true
				this.xTransform = 0;
				this.yTransform = 0;
				this.zTransform = 0;
				this.screenshot = () => {
					html2canvas(document.querySelector("#main")).then(canvas => {
						document.body.appendChild(canvas)
					});
				}
				this.Transparency = true;
			}

			let electrodeController = function () {
				this.Display = true;
			}
			let elecCtrl = new electrodeController();
			let electrodeMenu = gui.addFolder('Electrodes');
			let electrodeToggler = electrodeMenu.add(elecCtrl, "Display", true).listen()
			electrodeToggler.onChange(val => {
				elecs.visible = val;

			})

			let text = new meshController();
			let meshToggler = gui.addFolder('Mesh');
			let transparencyToggler = gui.addFolder('Transparency');
			let transToggler = transparencyToggler.add(text,'Transparency', true).listen();
			// let screenshot = meshToggler.add(text, "screenshot").listen();
			let fullMeshToggler = meshToggler.add(text, "Mesh", false).listen()
			let cortexMeshToggler = meshToggler.add(text, "Cortex", false).listen()
			let wmMeshToggler = meshToggler.add(text, "WM", false).listen()
			let subMeshToggler = meshToggler.add(text, "Substructures", false).listen()

			transToggler.onChange(val=> {
				if(val == true){
					brainScene.traverse(child => {
						if (child instanceof THREE.Mesh && child.parent.name != "Electrodes") {
							child.material.transparent = true;
							child.material.opacity = .5
						}
					})	
				}
				else{
					brainScene.traverse(child => {
						if (child instanceof THREE.Mesh && child.parent.name != "Electrodes") {
							child.material.transparent = false;
						}
					})	
				}
			})


			fullMeshToggler.onChange(val => {
				if (val == false) {
					wm.visible = false
					gyri.visible = false
					substructures.visible = false
					text.Cortex = false
					text.WM = false
					text.Substructures = false
				}
				else {
					wm.visible = true
					gyri.visible = true
					substructures.visible = true
					text.Cortex = true
					text.WM = true
					text.Substructures = true
				}
			})
			cortexMeshToggler.onChange(val => {
				if (val == false) {
					gyri.visible = false
				}
				else {
					gyri.visible = true
				}
			})
			wmMeshToggler.onChange(val => {
				if (val == false) {
					wm.visible = false
				}
				else {
					wm.visible = true
				}
			})
			subMeshToggler.onChange(val => {
				if (val == false) {
					substructures.visible = false
				}
				else {
					substructures.visible = true
				}
			})



			function updateClipPlane(refObj, clipPlane) {
				const stackHelper = refObj.stackHelper;
				const camera = refObj.camera;
				let vertices = stackHelper.slice.geometry.vertices;
				let p1 = new THREE.Vector3(vertices[0].x, vertices[0].y, vertices[0].z).applyMatrix4(
					stackHelper._stack.ijk2LPS
				);
				let p2 = new THREE.Vector3(vertices[1].x, vertices[1].y, vertices[1].z).applyMatrix4(
					stackHelper._stack.ijk2LPS
				);
				let p3 = new THREE.Vector3(vertices[2].x, vertices[2].y, vertices[2].z).applyMatrix4(
					stackHelper._stack.ijk2LPS
				);

				clipPlane.setFromCoplanarPoints(p1, p2, p3);

				let cameraDirection = new THREE.Vector3(1, 1, 1);
				cameraDirection.applyQuaternion(camera.quaternion);

				if (cameraDirection.dot(clipPlane.normal) > 0) {
					clipPlane.negate();
				}
			}

			function onYellowChanged() {
				updateClipPlane(r2, clipPlane2);
			}

			yellowChanged.onChange(onYellowChanged);

			function onRedChanged() {
				updateClipPlane(r1, clipPlane1);

				// if (redContourHelper) {
				// 	redContourHelper.geometry = r1.stackHelper.slice.geometry;
				// }
			}

			redChanged.onChange(onRedChanged);

			function onGreenChanged() {
				updateClipPlane(r3, clipPlane3);
			}

			greenChanged.onChange(onGreenChanged);

			function onDoubleClick(event) {
				const canvas = event.target.parentElement;
				const id = event.target.id;
				const mouse = {
					x: ((event.clientX - canvas.offsetLeft) / canvas.clientWidth) * 2 - 1,
					y: -((event.clientY - canvas.offsetTop) / canvas.clientHeight) * 2 + 1,
				};
				//
				let camera = null;
				let stackHelper = null;
				let scene = null;
				switch (id) {
					case '0':
						camera = r0.camera;
						stackHelper = r1.stackHelper;
						scene = r0.scene;
						break;
					case '1':
						camera = r1.camera;
						stackHelper = r1.stackHelper;
						scene = r1.scene;
						break;
					case '2':
						camera = r2.camera;
						stackHelper = r2.stackHelper;
						scene = r2.scene;
						break;
					case '3':
						camera = r3.camera;
						stackHelper = r3.stackHelper;
						scene = r3.scene;
						break;
				}

				const raycaster = new THREE.Raycaster();
				raycaster.setFromCamera(mouse, camera);

				const intersects = raycaster.intersectObjects(scene.children, true);
				if (intersects.length > 0) {
					let ijk = UtilsCore.worldToData(stackHelper.stack.lps2IJK, intersects[0].point);

					r1.stackHelper.index = ijk.getComponent((r1.stackHelper.orientation + 2) % 3);
					r2.stackHelper.index = ijk.getComponent((r2.stackHelper.orientation + 2) % 3);
					r3.stackHelper.index = ijk.getComponent((r3.stackHelper.orientation + 2) % 3);

					onGreenChanged();
					onRedChanged();
					onYellowChanged();
				}
			}

			// event listeners
			r0.domElement.addEventListener('dblclick', onDoubleClick);
			r1.domElement.addEventListener('dblclick', onDoubleClick);
			r2.domElement.addEventListener('dblclick', onDoubleClick);
			r3.domElement.addEventListener('dblclick', onDoubleClick);


			function onScroll(event) {
				const id = event.target.domElement.id;
				let stackHelper = null;
				switch (id) {
					case 'r1':
						stackHelper = r1.stackHelper;
						break;
					case 'r2':
						stackHelper = r2.stackHelper;
						break;
					case 'r3':
						stackHelper = r3.stackHelper;
						break;
				}

				if (event.delta > 0) {
					if (stackHelper.index >= stackHelper.orientationMaxIndex - 1) {
						return false;
					}
					stackHelper.index += 1;
				} else {
					if (stackHelper.index <= 0) {
						return false;
					}
					stackHelper.index -= 1;
				}

				onGreenChanged();
				onRedChanged();
				onYellowChanged();
			}

			// event listeners
			r1.controls.addEventListener('OnScroll', onScroll);
			r2.controls.addEventListener('OnScroll', onScroll);
			r3.controls.addEventListener('OnScroll', onScroll);

			function windowResize2D(rendererObj) {
				rendererObj.camera.canvas = {
					width: rendererObj.domElement.clientWidth,
					height: rendererObj.domElement.clientHeight,
				};
				rendererObj.camera.fitBox(2, 1);
				rendererObj.renderer.setSize(
					rendererObj.domElement.clientWidth,
					rendererObj.domElement.clientHeight
				);

				// update info to draw borders properly
				rendererObj.stackHelper.slice.canvasWidth = rendererObj.domElement.clientWidth;
				rendererObj.stackHelper.slice.canvasHeight = rendererObj.domElement.clientHeight;

			}

			function onWindowResize() {
				// update 3D
				r0.camera.aspect = r0.domElement.clientWidth / r0.domElement.clientHeight;
				r0.camera.updateProjectionMatrix();
				r0.renderer.setSize(r0.domElement.clientWidth, r0.domElement.clientHeight);

				// update 2d
				windowResize2D(r1);
				windowResize2D(r2);
				windowResize2D(r3);
			}

			window.addEventListener('resize', onWindowResize, false);


			let RASToLPS = new THREE.Matrix4();
			const worldCenter = r1.stackHelper.stack.worldCenter();
			RASToLPS.set(
				-1,
				0,
				0,
				worldCenter.x,
				0,
				-1,
				0,
				worldCenter.y,
				0,
				0,
				1,
				worldCenter.z,
				0,
				0,
				0,
				1
			);


			function loadElectrodes(elec, i) {

				let meshOpacity = .8;
				let transparency = false;
				if (elec.parent.name == "Electrodes") {
					transparency = false;
				}
				else {
					meshOpacity = .5
					transparency = true;
				}

				data[i].scene = new THREE.Scene();
				data[i].materialFront = new THREE.MeshBasicMaterial({
					color: elec.material.color,
					side: THREE.FrontSide,
					depthWrite: true,
					opacity: 0,
					transparent: true,
					clippingPlanes: [],
				});
				data[i].materialBack = new THREE.MeshBasicMaterial({
					color: elec.material.color,
					side: THREE.BackSide,
					depthWrite: true,
					opacity: meshOpacity,
					transparent: true,
					clippingPlanes: [],
				});
				data[i].meshFront = new THREE.Mesh(elec.geometry, data[i].materialFront);
				data[i].meshBack = new THREE.Mesh(elec.geometry, data[i].materialBack);
				data[i].meshFront.position.set(elec.position.x, elec.position.y, elec.position.z)
				data[i].meshBack.position.set(elec.position.x, elec.position.y, elec.position.z)
				if (elec.parent.name == "Electrodes") {
					data[i].meshFront.scale.set(2, 2, 2)
					data[i].meshBack.scale.set(2, 2, 2)
				}
				data[i].scene.add(data[i].meshFront);
				data[i].scene.add(data[i].meshBack);
				sceneClip.rotation.set(-Math.PI / 2, 0, Math.PI)
				sceneClip.add(data[i].scene);
			}

			onGreenChanged();
			onRedChanged();
			onYellowChanged();
			let load3DBrain_gltf = () => {
				return new Promise((resolve, reject) => {
					let loader = new THREE.GLTFLoader()
					loader.load(`./static/reconstruction.glb`, object3d => {
						object3d.scene.traverse(child => {
							if (child instanceof THREE.Mesh && child.parent.name != "Electrodes") {
								child.material.transparent = true;
								child.material.opacity = .5
							}
							if (child instanceof THREE.Mesh && child.parent.name == "Electrodes") {
								let colorVals = {
									"LPI": [0.000 ,0.000 ,1.000],
									"LMI": [1.000 ,0.000 ,0.000],
									"LAI": [0.471 ,0.000 ,1.000],
									"LFC": [1.000 ,0.569 ,0.482],
									"LBT": [0.000 ,0.745 ,1.000],
									"LPH": [1.000 ,0.000 ,1.000],
									"LTA": [0.000 ,0.337 ,0.000],
									"LAH": [1.000 ,0.431 ,0.000],
									"LA": [0.388 ,0.263 ,0.137],
									"LMH": [0.000 ,1.000 ,0.000],
									"LTP": [1.000 ,0.851 ,0.000]
																	}
								let [elecName, elecNumber] = child.name.replace(/\'/g, '').split(/(\d+)/).filter(Boolean);
								if (colorVals[elecName] != undefined) {
									child.material = new THREE.MeshLambertMaterial({
										emissive: new THREE.Color(colorVals[elecName][0], colorVals[elecName][1], colorVals[elecName][2])
									});
								} else {
									child.material = new THREE.MeshLambertMaterial({
										emissive: new THREE.Color(1, 1, 1)
									});
								}
								child.scale.set(.75, .75, .75)
							}
						})
						object3d.scene.rotation.set(0, 0, Math.PI)
						r0.scene.add(object3d.scene);
						brainScene = object3d.scene
						wm = brainScene.children[3]
						gyri = brainScene.children[2]
						substructures = brainScene.children[1]
						elecs = object3d.scene.children[0]
						elecs.rotation.set(0, 0, Math.PI)
						elecs.position.set(128, 128, 128)
						resolve(brainScene);
						console.log(object3d)
					})

				})
			}
			load3DBrain_gltf().then((scene) => {
				let index = 0;
				scene.traverse(child => {
					if (child instanceof THREE.Mesh) {
						if (child.parent.parent.name != "WhiteMatter") {
							data[index] = {}
							loadElectrodes(child, index);
							index++;
						}
					}
				})
				ready = true;
				render();
			})
		})
		.catch(function (error) {
			window.console.log('oops... something went wrong...');
			window.console.log(error);
		});
		// $('#alertModal').modal({
		// 	show:true
		// })
	
};
