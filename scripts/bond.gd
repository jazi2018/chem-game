extends Node2D

var start_atom: Atom
var end_atom: Atom
@onready var start: Line2D = $Start
@onready var end: Line2D = $End
#var bond_order: int = 1
#var gradient := Gradient.new()
#set gradient interpolation_mode to Gradient.INTERP_CONSTANT
#then add_point at 0.5 (halfway thru the line) and set color to relevant color on each side

func setup(from: Atom, to: Atom) -> void: #order: int = 1
	start_atom = from
	end_atom = to
	#bond_order = order

func _ready() -> void:
	start.default_color = start_atom.get_color()
	end.default_color = end_atom.get_color()
	setup_line()

func setup_line() -> void:
	var midpoint := (start_atom.position + end_atom.position) * 0.5
	var start_radius = 0.0 if start_atom.is_carbon else start_atom.get_label_radius()
	var end_radius = 0.0 if end_atom.is_carbon else end_atom.get_label_radius()
	var uv := (end_atom.position - start_atom.position).normalized()
	
	#first line
	start.clear_points()
	start.add_point(start_atom.position + uv * start_radius)
	start.add_point(midpoint)
	
	#second line
	end.clear_points()
	end.add_point(midpoint)
	end.add_point(end_atom.position - uv * end_radius)

func _process(_delta: float) -> void:
	#self.points = [
		#start_atom.position,
		#end_atom.position
	#]
	setup_line() #NOTE: CHANGE THIS AT SOME POINT - this is a TEMPORARY FIX for testing
	#should only be run WHEN THE ATOM IS MOVED LOL
	pass
