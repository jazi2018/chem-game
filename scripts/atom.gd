extends Node2D
class_name Atom

var symbol: String = "C" #carbon by default
var index: int = 0
#rendering variables
var color: Color = Color.BLACK
var is_carbon: bool = true #carbon flag

@onready var TEMP_IDX_LABEL : Label = $DEBUG_idx_label
@onready var label : RichTextLabel = $Label

func setup(symb: String, idx: int) -> void:
	symbol = symb
	if symbol != "C":
		is_carbon = false
	index = idx
	set_color()

func _ready() -> void:
	#label atom if not carbon and apply subscript affects
	label.install_effect(preload("res://scripts/subscript.gd").new())
	if not is_carbon:
		label.set_text(symbol)
		label.add_theme_color_override("default_color", color)
	### DEBUG ###
	TEMP_IDX_LABEL.set_text(str(index))
	### DEBUG ###

func set_color() -> void:
	#assess only first 2 characters of symbol
	if len(symbol) > 1:
		var temp = symbol[0]
		match temp:
			'C' : color = Color.GREEN #Cl
			'B' : color = Color.BROWN #Br
			'N' : color = Color.BLUE #N / NH / ...
			'O' : color = Color.RED #O / OH
			'S' : color = Color.YELLOW # S / SH / ...
			'I' : color = Color.REBECCA_PURPLE # I / IH? / ...
			_ : color = Color.BLACK
	else:
		match symbol:
			'C' : color = Color.BLACK
			'N' : color = Color.BLUE
			'O' : color = Color.RED
			'S' : color = Color.YELLOW
			'I' : color = Color.REBECCA_PURPLE
			'F' : color = Color.MAGENTA
			_ : color = Color.BLACK

func get_color() -> Color:
	return color

func get_label_radius() -> float:
	var size = label.get_minimum_size()
	#centered, so * 0.6 for slightly larger than text cutoff
	return max(size.x, size.y) * 0.5
