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
    mol = Chem.AddHs(mol) #adds ALL hydrogens - will filter them out later
    AllChem.Compute2DCoords(mol)

    include = set()
    #determine which atoms to include
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            include.add(atom.GetIdx())
        else:
            #only add if not attached to a carbon
            neighbors = atom.GetNeighbors()
            if neighbors and neighbors[0].GetSymbol() != 'C':
                include.add(atom.GetIdx())

    atoms = []
    for atom in mol.GetAtoms():
        #filter out non-included atoms
        idx = atom.GetIdx()
        if idx not in include:
            continue
        
        #set symbol to include hydrogens
        #TODO: Still need to implement fully - possibly need to include stereochem but we'll get there
        # symbol = atom.GetSymbol()
        # hydrogens = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')
        # if symbol != 'H' and hydrogens:
        #     symbol = f"{symbol}(H{hydrogens if hydrogens > 1 else ''})"

        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        atoms.append({
            'symbol': atom.GetSymbol(),
            'x': pos.x,
            'y': pos.y,
            'index': atom.GetIdx()
        })
    
    bonds = []
    for bond in mol.GetBonds():
        bgn = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        #filter out omitted atom bonds
        if bgn in include and end in include:
            bonds.append({
                'begin': bgn,
                'end': end,
                'order': bond.GetBondTypeAsDouble()
            })
    
    #DEBUG:
    print(atoms)
    print(f'\n{bonds}')
    
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