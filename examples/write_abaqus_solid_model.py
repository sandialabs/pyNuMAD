import pynumad as pynu
import numpy as np
from os.path import join

from pynumad.mesh_gen.mesh_gen import get_solid_mesh

#print("Notice: example write_abacus_solid_model.py temporarily down for updates.  To be restored shortly.")

abqFileName = "shellBlade.inp"
adhesiveMat = "Adhesive"

## Read blade data from yaml file
blade = pynu.Blade()
file_name = join("example_data","blade.yaml")
blade.read_yaml(file_name)

## Set the airfoil point resolution
for stat in blade.definition.stations:
    stat.airfoil.resample(n_samples=300)
    
blade.update_blade()
nStations = blade.geometry.coordinates.shape[2]
minTELengths = 0.001*np.ones(nStations)
blade.expand_blade_geometry_te(minTELengths)


## Set the target element size for the mesh
elementSize = 0.2

## Specify the elements per primary layer and generate mesh
layNumEls = [1,1,1]
bladeMesh = get_solid_mesh(blade,layNumEls,elementSize)

## Check for negative jacobians

failedEls = pynu.mesh_gen.mesh_tools.check_all_jacobians(bladeMesh['nodes'],bladeMesh['elements'])

totEls = len(bladeMesh['elements'])
newELabel = -np.ones(totEls,dtype=int)
lab = 0
for i in range(0,totEls):
    if(i not in failedEls):
        newELabel[i] = lab
        lab = lab + 1

## Write Abaqus input

outFile = open('solidBlade.inp','w')

outFile.write('*Part, name=Blade\n')
outFile.write('*Node\n')

i = 1
for nd in bladeMesh['nodes']:
    lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
    ln = ', '.join(lst)
    outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D8\n')
i = 0
for el in bladeMesh['elements']:
    if(newELabel[i] > -1 and el[7] != -1):
        el = el + 1
        lst = [str(newELabel[i]+1),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),str(el[6]),str(el[7]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D6\n')
i = 0
for el in bladeMesh['elements']:
    if(newELabel[i] > -1 and el[7] == -1):
        el[0:6] = el[0:6] + 1
        lst = [str(newELabel[i]+1),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
for es in bladeMesh['sets']['element']:
    ln = '*Elset, elset=set_' + es['name'] + '\n'
    outFile.write(ln)
    for el in es['labels']:
        if(newELabel[el] > -1):
            ln = '  ' + str(newELabel[el]+1) + '\n'
            outFile.write(ln)
        
for sec in bladeMesh['sections']:
    ln = '*Orientation, name=ori_' + sec['elementSet'] + '\n'
    outFile.write(ln)
    dataLn = list()
    for d in sec['xDir']:
        dataLn.append(str(d))
    for d in sec['xyDir']:
        dataLn.append(str(d))
    dataStr = ', '.join(dataLn) + '\n'
    outFile.write(dataStr)
    outFile.write('1, 0.\n')
    
for sec in bladeMesh['sections']:
    snm = sec['elementSet']
    ln = '*Solid Section, elset=set_' + snm + ', material=' + sec['material'] + ', orientation=ori_' + snm + '\n'
    outFile.write(ln)
    outFile.write(', \n')

outFile.write('*End Part\n')

outFile.write('*Part, name=Adhesive\n')
outFile.write('*Node\n')

i = 1
for nd in bladeMesh['adhesiveNds']:
    lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
    ln = ', '.join(lst)
    outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D8\n')
i = 1
for el in bladeMesh['adhesiveEls']:
    if(el[7] != -1):
        el = el + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),str(el[6]),str(el[7]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D6\n')
i = 1
for el in bladeMesh['adhesiveEls']:
    if(el[7] == -1):
        el[0:6] = el[0:6] + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1

es = bladeMesh['adhesiveElSet']
ln = '*Elset, elset=set_' + es['name'] + '\n'
outFile.write(ln)
for el in es['labels']:
    ln = '  ' + str(el+1) + '\n'
    outFile.write(ln)

outFile.write('*Orientation, name=global\n')
outFile.write('1., 0., 0., 0., 1., 0.\n')
outFile.write('1, 0.\n')

ln = '*Solid Section, elset=set_' + es['name'] + ', material=' + adhesiveMat + ', orientation=global\n'
outFile.write(ln)
outFile.write(', \n')

outFile.write('*End Part\n')

outFile.write('*Assembly, name=Assembly\n')
outFile.write('**\n')
outFile.write('*Instance, name=BladeInst, part=Blade\n')
outFile.write('*End Instance\n')
outFile.write('*Instance, name=AdhesiveInst, part=Adhesive\n')
outFile.write('*End Instance\n')
outFile.write('**\n')

for ns in bladeMesh['sets']['node']:
    ln = '*Nset, nset=set_' + ns['name'] + ', instance=BladeInst\n'
    outFile.write(ln)
    for nd in ns['labels']:
        ln = str(nd + 1) + ',\n'
        outFile.write(ln)
        
tiedAdNds = set()
tgtBlNds = set()

constraints = bladeMesh['constraints']
for c in constraints:
    terms = c['terms']
    for t in terms:
        lab = t['node'] + 1
        if(t['nodeSet'] == 'tiedMesh'):
            tiedAdNds.add(lab)
        else:
            tgtBlNds.add(lab)
            
for nd in tiedAdNds:
    ndstr = str(nd)
    ln = '*Nset, nset=a' + ndstr + ', instance=AdhesiveInst\n'
    outFile.write(ln)
    ln = ndstr + ',\n'
    outFile.write(ln)
    
for nd in tgtBlNds:
    ndstr = str(nd)
    ln = '*Nset, nset=b' + ndstr + ', instance=BladeInst\n'
    outFile.write(ln)
    ln = ndstr + ',\n'
    outFile.write(ln)
    
for c in constraints:
    numTerms = str(len(c['terms'])) + '\n'
    for i in range(1,4):
        outFile.write('*Equation\n')
        outFile.write(numTerms)
        for t in c['terms']:
            lab = t['node'] + 1
            if(t['nodeSet'] == 'tiedMesh'):
                ns = 'a' + str(lab)
            else:
                ns = 'b' + str(lab)
            ln = ', '.join([ns,str(i),str(t['coef'])]) + '\n'
            outFile.write(ln)
            
outFile.write('*End Assembly\n')

for mn in blade.definition.materials:
    mat = blade.definition.materials[mn]
    ln = '*Material, name=' + mat.name + '\n'
    outFile.write(ln)
    outFile.write('*Density\n')
    ln = str(mat.density) + ',\n'
    outFile.write(ln)
    outFile.write('*Elastic, type=ENGINEERING CONSTANTS\n')
    eProps = [str(mat.ex),str(mat.ey),str(mat.ez)]
    if(mat.type == "isotropic"):
        nu = str(mat.prxy)
        pr = [nu,nu,nu]
    else:
        pr = [str(mat.prxy),str(mat.prxz),str(mat.pryz)]
    eProps.extend(pr)
    eProps.extend([str(mat.gxy),str(mat.gxz),str(mat.gyz)])
    ln = ', '.join(eProps[0:8]) + '\n'
    outFile.write(ln)
    ln = eProps[8] + ',\n'
    outFile.write(ln)

outFile.write('*Boundary\n')
for i in range(1,4):
    sti = str(i)
    ln = 'set_RootNodes, ' + sti + ', ' + sti + '\n'
    outFile.write(ln)

outFile.close()    

## End write Abaqus