extends Line2D

var start_atom: Atom
var end_atom: Atom
#var bond_order: int = 1
#can use var gradient := Gradient.new()
#set gradient interpolation_mode to Gradient.INTERP_CONSTANT
#then add_point at 0.5 (halfway thru the line) and set color to relevant color on each side

func _init(from: Atom, to: Atom) -> void: #order: int = 1
	start_atom = from
	end_atom = to
	#bond_order = order

func _ready() -> void:
	setup_line()

func setup_line() -> void:
	clear_points()
	add_point(start_atom.position)
	add_point(end_atom.position)
