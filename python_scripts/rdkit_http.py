#SMILES parsing
from rdkit import Chem
from rdkit.Chem import AllChem

#http requests
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

app = FastAPI()

#data model for request payloads - json with string field "smiles"
#for ACTUAL game, probably posting reaction type and letting smiles be stored locally
#or can pass both smiles and reaction type - either way payload will be different
class SmilesInput(BaseModel):
    smiles: str

def parse_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    AllChem.Compute2DCoords(mol)

    mol_with_h = Chem.AddHs(mol) #add hydrogens to different molecules
    h_counts = {}
    for atom in mol_with_h.GetAtoms(): #get hydrogen counts on each non hydrogen atom - store in dict
        if atom.GetSymbol() != 'H':
            idx = atom.GetIdx()
            h_count = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
            h_counts[idx] = h_count

    #set symbol to include hydrogens
    #TODO: Still need to implement fully - possibly need to include stereochem but we'll get there

    atoms = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        symbol = atom.GetSymbol()
        h_count = h_counts.get(idx, 0) #get h count from dict - lets geometry ignore H count but still be recalled for symbol
        if h_count and symbol != 'C':
            symbol += f"H[sub]{h_count}[/sub]" if h_count > 1 else "H" #format with bbcode for godot
        pos = mol.GetConformer().GetAtomPosition(idx)
        atoms.append({
            'symbol': symbol,
            'x': pos.x,
            'y': pos.y,
            'index': idx
        })
    
    bonds = []
    for bond in mol.GetBonds():
        bgn = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        bonds.append({
            'begin': bgn,
            'end': end,
            'order': bond.GetBondTypeAsDouble()
        })
    
    #DEBUG:
    # print(atoms)
    # print(f'\n{bonds}')
    
    return {'atoms': atoms, 'bonds': bonds}

#define post endpoint
@app.post("/parse")
def parse(input: SmilesInput):
    try:
        result = parse_smiles(input.smiles)
        return result
    except Exception as e:
        print(e)
        raise HTTPException(status_code=400, detail=str(e))

if __name__ == '__main__':
    import uvicorn
    uvicorn.run("rdkit_http:app", host="127.0.0.1", port=8000)