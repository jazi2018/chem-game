extends Node2D
class_name Atom

var symbol: String = "C" #carbon by default
var index: int = 0
#rendering variables
var color: Color = Color.BLACK
var is_carbon: bool = true #carbon flag

@onready var label: Label = $Label

func setup(symb: String, idx: int) -> void:
	symbol = symb
	if symbol != "C":
		is_carbon = false
	index = idx
	set_color()

func _ready() -> void:
	#label atom if not carbon
	if not is_carbon:
		label.visible = true
		label.set_text(symbol)
		label.add_theme_color_override("font_color", color)

func set_color() -> void:
	#assess only first 2 characters of symbol
	var temp = symbol
	if len(symbol) > 1:
		temp = symbol.left(2)
	
	match temp:
		'C' : color = Color.BLACK
		'N' : color = Color.BLUE
		'O' : color = Color.RED
		'S' : color = Color.YELLOW
		'Cl' : color = Color.GREEN
		'Br' : color = Color.BROWN
		'I' : color = Color.REBECCA_PURPLE
		_ : color = Color.BLACK

func get_color() -> Color:
	return color

func get_label_radius() -> float:
	var size = label.get_minimum_size()
	#centered, so * 0.6 for slightly larger than text cutoff
	return max(size.x, size.y) * 0.5
