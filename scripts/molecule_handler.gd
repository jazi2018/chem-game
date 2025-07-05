extends Node2D

#temporary test string
var smiles: String = "CN(CC#C)C(=O)CCNC1=NC2=C(C(=C(C=C2C(=N1)N3CCN(CC3)C(=O)C=C)Cl)C4=C(C=CC=C4F)O)F"
@onready var http = $HTTPRequest
@onready var molecule_root = $MoleculeRoot

#scale for visualization spacing
@export var view_scale = 40.0
@export var atom_scene: PackedScene
@export var bond_scene: PackedScene

func _ready() -> void:
	#make request only temporarily here - eventually game controller will call make_request and update_smiles
	make_request()

func make_request() -> void:
	var body = JSON.stringify({"smiles": self.smiles})
	var headers = ["Content-Type: application/json"]
	http.request("http://127.0.0.1:8000/parse", headers, HTTPClient.METHOD_POST, body)

func update_smiles(new_smiles: String) -> void:
	self.smiles = new_smiles

func _on_httprequest_request_completed(_result: int, response_code: int,
_headers: PackedStringArray, body: PackedByteArray) -> void:
	if response_code != 200:
		print("Request failed: ", response_code)
		return

	var data: Dictionary = JSON.parse_string(body.get_string_from_utf8())
	if not data:
		print("Failed to parse JSON")
		return

	#draw_molecule(data)
	render_molecule(data["atoms"], data["bonds"])

func render_molecule(atoms_data: Array, bonds_data: Array, v_scale: float = 40.0) -> void:
	# Clear any existing visuals
	for child in molecule_root.get_children():
		remove_child(child)
		child.queue_free()

	# Map index → position
	var atom_list: Dictionary = {}

	for atom_info in atoms_data:
		var idx: int = atom_info["index"]
		var symbol: String = atom_info["symbol"]
		var pos: Vector2 = Vector2(atom_info["x"], -atom_info["y"]) * v_scale
		
		var atom = atom_scene.instantiate() as Node2D
		atom.position = pos
		atom.setup(symbol, idx)
		molecule_root.add_child(atom)
		atom_list[idx] = atom

	#draw bonds
	for bond_info in bonds_data:
		var i1: int = bond_info["begin"]
		var i2: int = bond_info["end"]
		var order: float = bond_info.get("order", 1.0)
		var stereo: int = bond_info["stereo"] #TODO: APPLY BOND STEREO TO BONDS INGAME THEN DONE WITH RENDERING!!!
		#ensure correct type of bond based on order
		for i in range(int(order)):
			var bond = bond_scene.instantiate() as Node2D
			bond.setup(atom_list[i1], atom_list[i2], i, order, stereo)
			molecule_root.add_child(bond)
