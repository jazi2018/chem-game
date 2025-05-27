extends Node2D
class_name Atom

var symbol: String = "C" #carbon by default
var index: int = 0
var color: Color = Color.BLACK

func setup(symb: String, idx: int) -> void:
	symbol = symb
	index = idx

func _ready() -> void:
	set_color()

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
