import pynumad as pynu
import numpy as np
import os
from os.path import join

from pynumad.shell.shell import get_shell_mesh

abqFileName = "shellBlade.inp"
adhesiveMat = 'Adhesive'

## Read blade data from yaml file
blade = pynu.Blade()
fileName = join("example_data","blade.yaml")
blade.read_yaml(fileName)

## Set the airfoil point resolution
for stat in blade.definition.stations:
    stat.airfoil.resample(n_samples=300)
    
#blade.generate_geometry()
blade.update_blade()
nStations = blade.geometry.coordinates.shape[2]
minTELengths = 0.001*np.ones(nStations)
blade.expand_blade_geometry_te(minTELengths)

## Set the target element size for the mesh
elementSize = 0.2

## Generate mesh
adhes = 1
bladeMesh = get_shell_mesh(blade, adhes, elementSize)

## Write Abaqus input

outFile = open(abqFileName,'w')

outFile.write('*Part, name=Blade\n')
outFile.write('*Node\n')

i = 1
for nd in bladeMesh['nodes']:
    lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
    ln = ', '.join(lst)
    outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=S4\n')
i = 1
for el in bladeMesh['elements']:
    if(el[3] != -1):
        el = el + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=S3\n')
i = 1
for el in bladeMesh['elements']:
    if(el[3] == -1):
        el[0:3] = el[0:3] + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
for es in bladeMesh['sets']['element']:
    ln = '*Elset, elset=set_' + es['name'] + '\n'
    outFile.write(ln)
    for el in es['labels']:
        ln = '  ' + str(el+1) + '\n'
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
    ln = '*Shell General Section, elset=set_' + snm + ', composite, orientation=ori_' + snm + ', offset=0.5, layup=lyp_' + snm + '\n'
    outFile.write(ln)
    for lay in sec['layup']:
        laylst = [str(lay[1])," ",lay[0],str(lay[2])]
        layStr = ', '.join(laylst) + '\n'
        outFile.write(layStr)

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
ln = '*Elset, elset=' + es['name'] + '\n'
outFile.write(ln)
for el in es['labels']:
    ln = '  ' + str(el+1) + '\n'
    outFile.write(ln)
    
outFile.write('*Orientation, name=global\n')
outFile.write('1., 0., 0., 0., 1., 0.\n')
outFile.write('1, 0.\n')

ln = '*Solid Section, elset=' + es['name'] + ', material=' + adhesiveMat + ', orientation=global\n'
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
    ln = '*Nset, nset=' + ns['name'] + ', instance=BladeInst\n'
    outFile.write(ln)
    for nd in ns['labels']:
        ln = str(nd + 1) + ',\n'
        outFile.write(ln)
        
outFile.write('*Nset, nset=allShellNds, instance=BladeInst, generate\n')
ln = '1, ' + str(len(bladeMesh['nodes'])) + ', 1\n'
outFile.write(ln)

outFile.write('*Elset, elset=allShellEls, instance=BladeInst, generate\n')
ln = '1, ' + str(len(bladeMesh['elements'])) + ', 1\n'
outFile.write(ln)

outFile.write('*Nset, nset=allAdhesiveNds, instance=AdhesiveInst, generate\n')
ln = '1, ' + str(len(bladeMesh['adhesiveNds'])) + ', 1\n'
outFile.write(ln)

outFile.write('*Elset, elset=allAdhesiveEls, instance=AdhesiveInst, generate\n')
ln = '1, ' + str(len(bladeMesh['adhesiveEls'])) + ', 1\n'
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
        elif(t['nodeSet'] == 'targetMesh'):
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
for i in range(1,7):
    sti = str(i)
    ln = 'RootNodes, ' + sti + ', ' + sti + '\n'
    outFile.write(ln)

outFile.close()    

## End write Abaqus
