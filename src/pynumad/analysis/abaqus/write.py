# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:53:59 2024

@author: evaande
"""

from pynumad.objects.blade import *

def write_shell_general(fileName,bladeMesh):
    outFile = open(fileName,'w')

    outFile.write('*Part, name=part1\n')
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
        outFile.write('3, 0.\n')
        
    for sec in bladeMesh['sections']:
        snm = sec['elementSet']
        ln = '*Shell General Section, elset=set_' + snm + ', composite, orientation=ori_' + snm + ', offset=0.0, layup=lyp_' + snm + '\n'
        outFile.write(ln)
        for lay in sec['layup']:
            laylst = [str(lay[1])," ",lay[0],str(lay[2])]
            layStr = ', '.join(laylst) + '\n'
            outFile.write(layStr)

    outFile.write('*End Part\n')

    outFile.write('*Assembly, name=Assembly\n')
    outFile.write('**\n')
    outFile.write('*Instance, name=part1Inst, part=part1\n')
    outFile.write('*End Instance\n')
    outFile.write('**\n')

    for ns in bladeMesh['sets']['node']:
        ln = '*Nset, nset=' + ns['name'] + ', instance=part1Inst\n'
        outFile.write(ln)
        for nd in ns['labels']:
            ln = str(nd + 1) + ',\n'
            outFile.write(ln)
                
    outFile.write('*End Assembly\n')

    for mat in bladeMesh['materials']:
        ln = '*Material, name=' + mat['name'] + '\n'
        outFile.write(ln)
        outFile.write('*Density\n')
        ln = str(mat['density']) + ',\n'
        outFile.write(ln)
        outFile.write('*Elastic, type=ENGINEERING CONSTANTS\n')
        E = mat['elastic']['E']
        nu = mat['elastic']['nu']
        G = mat['elastic']['G']
        eProps = [str(E[0]),str(E[1]),str(E[2])]
        eProps.extend([str(nu[0]),str(nu[1]),str(nu[2])])
        eProps.extend([str(G[0]),str(G[1]),str(G[2])])
        ln = ', '.join(eProps[0:8]) + '\n'
        outFile.write(ln)
        ln = eProps[8] + ',\n'
        outFile.write(ln)

    outFile.close()
    return
    
def write_solid_general(fileName,bladeMesh):
    outFile = open(fileName,'w')

    outFile.write('*Part, name=part1\n')
    outFile.write('*Node\n')

    i = 1
    for nd in bladeMesh['nodes']:
        lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
        i = i + 1
        
    outFile.write('*Element, type=C3D8I\n')
    i = 1
    for el in bladeMesh['elements']:
        if(el[6] != -1):
            el = el + 1
            lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),str(el[6]),str(el[7]),'\n']
            ln = ', '.join(lst)
            outFile.write(ln)
        i = i + 1
        
    outFile.write('*Element, type=C3D6\n')
    i = 1
    for el in bladeMesh['elements']:
        if(el[6] == -1):
            el[0:6] = el[0:6] + 1
            lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),'\n']
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
        outFile.write('3, 0.\n')
        
    for sec in bladeMesh['sections']:
        snm = sec['elementSet']
        ln = '*Solid Section, elset=set_' + snm + ', orientation=ori_' + snm + ', material=' + sec['material'] + '\n'
        outFile.write(ln)
        outFile.write(', \n')

    outFile.write('*End Part\n')

    outFile.write('*Assembly, name=Assembly\n')
    outFile.write('**\n')
    outFile.write('*Instance, name=part1Inst, part=part1\n')
    outFile.write('*End Instance\n')
    outFile.write('**\n')

    try:
        for ns in bladeMesh['sets']['node']:
            ln = '*Nset, nset=' + ns['name'] + ', instance=part1Inst\n'
            outFile.write(ln)
            for nd in ns['labels']:
                ln = str(nd + 1) + ',\n'
                outFile.write(ln)
    except:
        pass
        
    for ni in range(0, len(bladeMesh['nodes'])):
        sni = str(ni+1)
        ln = '*Nset, nset=n_' + sni + ', instance=part1Inst\n'
        outFile.write(ln)
        ln = sni + ',\n'
        outFile.write(ln)
                
    outFile.write('*End Assembly\n')

    for mat in bladeMesh['materials']:
        ln = '*Material, name=' + mat['name'] + '\n'
        outFile.write(ln)
        outFile.write('*Density\n')
        ln = str(mat['density']) + ',\n'
        outFile.write(ln)
        outFile.write('*Elastic, type=ENGINEERING CONSTANTS\n')
        E = mat['elastic']['E']
        nu = mat['elastic']['nu']
        G = mat['elastic']['G']
        eProps = [str(E[0]),str(E[1]),str(E[2])]
        eProps.extend([str(nu[0]),str(nu[1]),str(nu[2])])
        eProps.extend([str(G[0]),str(G[1]),str(G[2])])
        ln = ', '.join(eProps[0:8]) + '\n'
        outFile.write(ln)
        ln = eProps[8] + ',\n'
        outFile.write(ln)

    outFile.close()
    return

def write_shell_modal_input(fileName,blade,bladeMesh,adhesiveMat,numModes):
    outFile = open(fileName,'w')

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
        outFile.write('3, 0.\n')
        
    for sec in bladeMesh['sections']:
        snm = sec['elementSet']
        ln = '*Shell General Section, elset=set_' + snm + ', composite, orientation=ori_' + snm + ', offset=0.5, layup=lyp_' + snm + '\n'
        outFile.write(ln)
        for lay in sec['layup']:
            laylst = [str(lay[1])," ",lay[0],str(lay[2])]
            layStr = ', '.join(laylst) + '\n'
            outFile.write(layStr)

    outFile.write('*End Part\n')

    if('adhesiveNds' in bladeMesh):
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
    if('adhesiveNds' in bladeMesh):
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

    if('adhesiveNds' in bladeMesh):
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
        
    outFile.write('*Step, name=modalAnalysis, nlgeom=NO, perturbation\n')
    outFile.write('*Frequency, eigensolver=Lanczos, sim, acoustic coupling=on, normalization=mass\n')
    ln = str(numModes) + ', , , , , \n'
    outFile.write(ln)
    outFile.write('*Restart, write, frequency=0\n')
    outFile.write('*Output, field\n')
    outFile.write('*Node Output\n')
    outFile.write('U,\n')
    outFile.write('*Element Output, position=CENTROIDAL, directions=YES\n')
    outFile.write('E, S\n')
    outFile.write('*End Step\n')

    outFile.close()
    return

def write_damping_script(scriptFile,inputFile,dampFile):
    outFile = open(scriptFile,'w')

    outFile.write("from abaqus import *\n")
    outFile.write("from abaqusConstants import *\n")
    outFile.write("\n")

    outFile.write("print('creating job')\n")
    ln = "mdb.JobFromInputFile(name='modalAnalysis', inputFileName='" + inputFile + "', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)\n"
    outFile.write(ln)
    outFile.write("print('submitting job')\n")
    outFile.write("mdb.jobs['modalAnalysis'].submit(consistencyChecking=OFF)\n")
    outFile.write("mdb.jobs['modalAnalysis'].waitForCompletion()\n")
    outFile.write("\n")
    outFile.write("print('extracting results')\n")
    outFile.write("jobDB = session.openOdb(name='modalAnalysis.odb')\n")
    outFile.write("\n")
    ln = "outFile = open('" + dampFile + "','w')\n"
    outFile.write(ln)
    outFile.write("outFile.write('modes:\\n')\n")
    outFile.write("for fi, f in enumerate(jobDB.steps['modalAnalysis'].frames):\n")
    outFile.write("    if(fi > 0):\n")
    outFile.write("        mdStr = '  mode_' + str(fi) + ':\\n'\n")
    outFile.write("        outFile.write(mdStr)\n")
    outFile.write("        label = list()\n")
    outFile.write("        part = list()\n")
    outFile.write("        layer = list()\n")
    outFile.write("        stress = list()\n")
    outFile.write("        strain = list()\n")
    outFile.write("        stressVals = f.fieldOutputs['S'].values\n")
    outFile.write("        strainVals = f.fieldOutputs['E'].values\n")
    outFile.write("        for i, v in enumerate(stressVals):\n")
    outFile.write("            label.append(str(v.elementLabel))\n")
    outFile.write("            part.append(v.instance.name)\n")
    outFile.write("            try:\n")
    outFile.write("                desc = v.sectionPoint.description\n")
    outFile.write("                lst = desc.split('Layer =')\n")
    outFile.write("                layer.append(lst[1])\n")
    outFile.write("            except:\n")
    outFile.write("                layer.append('-1')\n")
    outFile.write("            if(len(v.data) == 4):\n")
    outFile.write("                s = v.data\n")
    outFile.write("                stress.append(str([s[0],s[1],s[2],s[3],0.0,0.0]))\n")
    outFile.write("                e = strainVals[i].data\n")
    outFile.write("                strain.append(str([e[0],e[1],e[2],e[3],0.0,0.0]))\n")
    outFile.write("            else:\n")
    outFile.write("                s = v.data\n")
    outFile.write("                stress.append(str([s[0],s[1],s[2],s[3],s[4],s[5]]))\n")
    outFile.write("                e = strainVals[i].data\n")
    outFile.write("                strain.append(str([e[0],e[1],e[2],e[3],e[4],e[5]]))\n")
    outFile.write("        outFile.write('    label: [')\n")
    outFile.write("        outFile.write(','.join(label))\n")
    outFile.write("        outFile.write(']\\n')\n")
    outFile.write("        outFile.write('    part: [')\n")
    outFile.write("        outFile.write(','.join(part))\n")
    outFile.write("        outFile.write(']\\n')\n")
    outFile.write("        outFile.write('    layer: [')\n")
    outFile.write("        outFile.write(','.join(layer))\n")
    outFile.write("        outFile.write(']\\n')\n")
    outFile.write("        outFile.write('    stress: [')\n")
    outFile.write("        outFile.write(','.join(stress))\n")
    outFile.write("        outFile.write(']\\n')\n")
    outFile.write("        outFile.write('    strain: [')\n")
    outFile.write("        outFile.write(','.join(strain))\n")
    outFile.write("        outFile.write(']\\n')\n")
    outFile.write("outFile.write('\\n')\n")
    outFile.write("outFile.close()\n")
    outFile.write("\n")

    outFile.close()
    return