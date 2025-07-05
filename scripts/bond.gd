extends Node2D

var start_atom: Atom
var end_atom: Atom
var offset_idx: int = 0
var bond_order: float = 1.0
var stereo: LineType = LineType.NONE
enum LineType {NONE, WEDGE, DASH}

@onready var start: Line2D = $Start
@onready var end: Line2D = $End

@export var spacing: float = 6.0
@export var dash_texture: Texture

func setup(from: Atom, to: Atom, offset: int, order: float = 1.0, ster: int = 0) -> void:
	start_atom = from
	end_atom = to
	offset_idx = offset
	bond_order = order
	match ster:
		1:
			self.stereo = LineType.WEDGE
		2:
			self.stereo = LineType.DASH
		_:
			self.stereo = LineType.NONE

func _ready() -> void:
	start.default_color = start_atom.get_color()
	end.default_color = end_atom.get_color()
	#create wedge shape
	if self.stereo == LineType.WEDGE or self.stereo == LineType.DASH:
		#setup curves
		var start_curve := Curve.new()
		start_curve.add_point(Vector2(0, 0.3))
		start_curve.add_point(Vector2(1, 0.7))
		var end_curve := Curve.new()
		end_curve.add_point(Vector2(0, 0.7))
		end_curve.add_point(Vector2(1, 1))
		#assign curves
		start.width = 10.0
		start.width_curve = start_curve
		end.width = 10.0
		end.width_curve = end_curve
	#assign texture if necessary
	if self.stereo == LineType.DASH:
		start.texture = dash_texture
		end.texture = dash_texture
	setup_line()

func setup_line() -> void:
	var start_radius = 0.0 if start_atom.is_carbon else start_atom.get_label_radius()
	var end_radius = 0.0 if end_atom.is_carbon else end_atom.get_label_radius()
	var uv := (end_atom.position - start_atom.position).normalized()
	var normal := Vector2(-uv.y, uv.x)
	
	var offset := 0.0
	#calculate offset
	if bond_order == 2.0:
		offset = spacing * (offset_idx * 2 - 1) * 0.5 #[-spacing/2, spacing/2]
	elif bond_order == 3.0:
		offset = spacing * (offset_idx - 1) #[-spacing, 0, spacing]
	
	var midpoint := (start_atom.position + end_atom.position) * 0.5
	var perp_offset = normal * offset
	#first line
	start.clear_points()
	start.add_point(start_atom.position + uv * start_radius + perp_offset)
	start.add_point(midpoint + perp_offset)
	
	#second line
	end.clear_points()
	end.add_point(midpoint + perp_offset)
	end.add_point(end_atom.position - uv * end_radius + perp_offset)

func _process(_delta: float) -> void:
	#self.points = [
		#start_atom.position,
		#end_atom.position
	#]
	setup_line() #NOTE: CHANGE THIS AT SOME POINT - this is a TEMPORARY FIX for testing
	#should only be run WHEN THE ATOM IS MOVED LOL
	#NOTE: Maybe add a setting at some point that can disable moving atoms for computers
	#that are not very performative
	pass
