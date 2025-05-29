@tool
extends RichTextEffect
class_name SubscriptEffect

const bbcode = "sub"

func _process_custom_fx(char_fx: CharFXTransform) -> bool:
	#adjust subscript offset
	char_fx.offset.y += 13
	char_fx.offset.x += 7

	#transform scaled
	var scale_factor = 0.75
	char_fx.transform = char_fx.transform.scaled(Vector2(scale_factor, scale_factor))

	return true
