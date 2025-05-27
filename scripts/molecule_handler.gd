extends Node2D

@onready var http = $HTTPRequest
@onready var molecule_root = $MoleculeRoot

#scale for visualization spacing
@export var view_scale = 40.0

func _ready() -> void:
	#TEMP SMILES TEST
	var body = JSON.stringify({"smiles": "CC(C(=O)O)N"})
	var headers = ["Content-Type: application/json"]
	http.request("http://127.0.0.1:8000/parse", headers, HTTPClient.METHOD_POST, body) 

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

func render_molecule(atoms_data: Array, bonds_data: Array, scale: float = 40.0) -> void:
	# Clear any existing visuals
	for child in molecule_root.get_children():
		remove_child(child)
		child.queue_free()

	# Map index → position
	var atom_positions: Dictionary = {}

	for atom_info in atoms_data:
		var idx: int = atom_info["index"]
		var symbol: String = atom_info["symbol"]
		var pos: Vector2 = Vector2(atom_info["x"], -atom_info["y"]) * scale
		atom_positions[idx] = pos

		# Only show labels for non-carbon atoms
		if symbol != "C":
			var label := Label.new()
			label.text = symbol
			#label.theme_override_fonts/font = load("res://Fonts/ChemFont.tres") # optional
			label.add_theme_color_override("font_color", Color.BLACK)
			label.anchor_left = 0.5
			label.anchor_top = 0.5
			label.position = pos - Vector2(10, 10)
			add_child(label)

	# Draw bonds
	for bond_info in bonds_data:
		var i1: int = bond_info["begin"]
		var i2: int = bond_info["end"]
		var order: float = bond_info.get("order", 1.0)

		var p1: Vector2 = atom_positions.get(i1, Vector2.ZERO)
		var p2: Vector2 = atom_positions.get(i2, Vector2.ZERO)

		var bond := Line2D.new()
		bond.width = 2
		bond.default_color = Color(0, 0, 0)
		bond.add_point(p1)
		bond.add_point(p2)
		add_child(bond)


#func draw_molecule(data: Dictionary) -> void:
	## Clear previous visuals
	#for child: Node in molecule_root.get_children():
		#child.queue_free()
#
	#var atom_nodes = {}  # Map from index -> position
#
	## Draw atoms
	#for atom in data["atoms"]:
		#var pos = Vector2(atom["x"], atom["y"]) * view_scale
		#var atom_node = create_atom(atom["symbol"], pos)
		#molecule_root.add_child(atom_node)
		#atom_nodes[atom["index"]] = pos
#
	## Draw bonds
	#for bond in data["bonds"]:
		#var start = atom_nodes[bond["begin"]]
		#var end = atom_nodes[bond["end"]]
		#var bond_node = create_bond(start, end)
		#molecule_root.add_child(bond_node)
#
#func create_atom(symbol: String, pos: Vector2) -> Node2D:
	#var atom_node = Node2D.new()
	#atom_node.position = pos
#
	## Draw a circle
	##atom_node.connect("draw", Callable(atom_node, "_on_draw_atom"))
	#atom_node.set_script(_make_atom_draw_script())
	#
	#var label = Label.new()
	#label.text = symbol
	#label.align = HORIZONTAL_ALIGNMENT_CENTER
	#label.vertical_alignment = VERTICAL_ALIGNMENT_CENTER
	#label.position = Vector2(-10, -10)
	#atom_node.add_child(label)
#
	#return atom_node
#
#func _make_atom_draw_script():
	#var script := GDScript.new()
	#script.source_code = '''
#extends Node2D
#
#func _draw():
	#draw_circle(Vector2.ZERO, 12, Color(0.8, 0.8, 1))
#'''
	#script.reload()
	#return script
#
#func create_bond(start: Vector2, end: Vector2) -> Line2D:
	#var line = Line2D.new()
	#line.width = 2.0
	#line.default_color = Color(0.2, 0.2, 0.2)
	#line.points = [start, end]
	#return line
