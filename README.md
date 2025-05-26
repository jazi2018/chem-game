#Some kinda chemistry video game
##Idea
A roguelite video game in which the player attempts to use / develop organic chemistry
knowledge to last as long as they can in a lab which is inevitably doomed to fail. (think Balatro)

Players will be given access to various reagents, which will enable them to manipulate
a compound to meet criterion provided to them.

##Implementation
This game uses various python libraries to enable certain functionalities.
	- [RDKit](https://www.rdkit.org/)
	- [FastAPI](https://fastapi.tiangolo.com/)
	- [uvicorn](https://www.uvicorn.org/)
The most central among these is RDKit. It enables the utilization of SMILES to both visualize
and manipulate molecules in a way that is mostly accurate to organic chemistry taught in college.

FastAPI and uvicorn allow the RDKit script to be run on a local server, which enables Godot (through the
HTTPNode) to interface with it and thus utilize its functionality.

This local server and RDKit exist only for the computational logic that is present in the game, in order
to maintain a high degree of chemical accuracy (after all - the game is meant to be educational)

All other functionality is / will be done locally within the Godot engine.

##TODO
	- Implement Atom and Bond node for better visuals and more extensible code
		- Atom node functionality:
			- Extends Node2D - relative position is important
			- Visual clarity
			- Animations
			- Varying colors for different elements
		- Bond node functionality:
			- Extends Line2D
			- Connected to neigboring atoms
			- Remains connected when atom moves (for animations)
			- Eventually need to implement stereochemistry (can just adjust thickness over length for wedge)
		- Children of Molecule, which should act as a container for atoms and bonds and store SMILE for whole molecule
	- Implement chemical reactions
		- Utilize reaction SMILES in RDKit
		- Cover reactions from Ochem 1 - potentially implement more later
		- Note all reagents that will be necessary
	- Actually develop game
		- Has to be fun or why play it?
		- Potentially randomly generated SMILE strings for when the run gets out of hand
		- Some kind of special powers that maintain educational value but enable the roguelite feeling
			- Thinking maybe a scoring system like balatro? Like "+$500 every time a stereocenter is flipped" and you get to use X reagents Y times
		- Up in the air rn - just want to get fundamentals up and running
