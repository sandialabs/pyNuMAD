########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import yaml
import numpy as np
from scipy.stats import mode

from pynumad.utils.misc_utils import LARCetaT, LARCetaL, _parse_data,fullKeysFromSubStrings
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.Component import Component
from pynumad.objects.Airfoil import Airfoil
from pynumad.objects.Material import Material


def yaml_to_blade(blade, filename: str, write_airfoils: bool = False):
    """
    This method writes blade information from a .yaml file to a Blade object.
    The yaml file is expected to be formatted according to the WindIO ontology.
    See https://windio.readthedocs.io/en/stable/source/turbine.html.

    Parameters
    ----------
    blade : Blade
    filename : string 
        path to .yaml file
    write_airfoils : bool
        Set true to write airfoil files while reading in data. Defaults to false.

    Returns
    -------
    blade : Blade
        input blade object populated with yaml data
    """

    # Read in yaml file as a nested dictionary
    with open(filename) as blade_yaml:
        # data = yaml.load(blade_yaml,Loader=yaml.FullLoader)
        data = yaml.load(blade_yaml,Loader=yaml.Loader)

    # Name some key subdata
    blade_outer_shape_bem = data['components']['blade']['outer_shape_bem']
    
    # older versions of wind ontology do not have 'outer_shape_bem' subsection for hub data
    try:
        hub_outer_shape_bem = data['components']['hub']['outer_shape_bem']
    except KeyError:
        hub_outer_shape_bem = data['components']['hub']
    
    blade_internal_structure = data['components']['blade']['internal_structure_2d_fem']
    af_data = data['airfoils']
    mat_data = data['materials']

    
    ### STATIONS / AIRFOILS
    _add_stations(blade, blade_outer_shape_bem, hub_outer_shape_bem, 
                    af_data, filename, write_airfoils)
    
    ### MATERIALS
    _add_materials(blade, mat_data)

    ## Blade Components
    N_layer_comp = len(blade_internal_structure['layers'])
    
    # Spar Cap Width and Offset
    # Obtain component name and index for hp and lp sides of sparcap


    for i in range(N_layer_comp):
        print(blade_internal_structure['layers'][i]['name'])
        if 'spar' in blade_internal_structure['layers'][i]['name'].lower():
            name = blade_internal_structure['layers'][i]['name']
            
            if 'suc' in blade_internal_structure['layers'][i]['side'].lower():
                spar_lp_index = i
                spar_lp_name = name
            if 'pres' in blade_internal_structure['layers'][i]['side'].lower():
                spar_hp_index = i
                spar_hp_name = name

    bladeParts=['layers','webs']
    # Make sure each blade.ispan has layer thicknesses and widths
    fullSpanGrid=np.array(blade_outer_shape_bem['chord']['grid'])
    nStations=len(fullSpanGrid)
    keysToModify={'offset_y_pa','thickness','fiber_orientation','width'}
    for partName in bladeParts:
        N_layer_comp = len(blade_internal_structure[partName])
        for currentLayer in range(N_layer_comp):
            layerKeys=set(blade_internal_structure[partName][currentLayer].keys())
            
            for currentKey in keysToModify.intersection(layerKeys):
                grid=blade_internal_structure[partName][currentLayer][currentKey]['grid']
                values=blade_internal_structure[partName][currentLayer][currentKey]['values']
                startStationLoc=grid[0]
                endStationLoc=grid[-1]

                subSpanGridIndex=np.where((fullSpanGrid>=startStationLoc) & (fullSpanGrid <=endStationLoc))[0]

                #iterpolate fullSpanGrid locations onto layer grid defined in the yamle file for the layer
                subSpanValues=interpolator_wrap(grid,values,fullSpanGrid[subSpanGridIndex],'pchip')
                fullSpanValues=np.zeros(nStations)

                fullSpanValues[subSpanGridIndex]=subSpanValues
            
                #Reset
                blade_internal_structure[partName][currentLayer][currentKey]['grid']=fullSpanGrid
                blade_internal_structure[partName][currentLayer][currentKey]['values']=fullSpanValues

    # Change layers from integer indexing to name indexing
    bladeStructureDict={}
    for i in range(len(blade_internal_structure['layers'])):
        bladeStructureDict[blade_internal_structure['layers'][i]['name'].lower()]=blade_internal_structure['layers'][i]

    #Spar caps
    sparCapKeys=fullKeysFromSubStrings(bladeStructureDict.keys(),['spar'])
    if len(sparCapKeys)==2:
        for iSparCap in range(2):
            if 'suc' in bladeStructureDict[sparCapKeys[iSparCap]]['side'].lower():
                lpSideIndex=iSparCap
            if 'pres' in bladeStructureDict[sparCapKeys[iSparCap]]['side'].lower():
                hpSideIndex=iSparCap
    else:
        raise ValueError('Incorrect number of spar cap components')
    

    blade.sparcapwidth_lp  = bladeStructureDict[sparCapKeys[lpSideIndex]]['width']['values']*1000
    blade.sparcapoffset_lp = bladeStructureDict[sparCapKeys[lpSideIndex]]['offset_y_pa']['values']*1000

    blade.sparcapwidth_hp  = bladeStructureDict[sparCapKeys[hpSideIndex]]['width']['values']*1000
    blade.sparcapoffset_hp = bladeStructureDict[sparCapKeys[hpSideIndex]]['offset_y_pa']['values']*1000

    # print(blade.sparcapwidth_lp-np.array(blade_internal_structure['layers'][spar_lp_index]['width']['values'])*1000)
    # print(blade.sparcapoffset_lp-np.array(blade_internal_structure['layers'][spar_lp_index]['offset_y_pa']['values'])*1000)
    # print(blade.sparcapwidth_hp-np.array(blade_internal_structure['layers'][spar_hp_index]['width']['values'])*1000)
    # print(blade.sparcapoffset_hp-np.array(blade_internal_structure['layers'][spar_hp_index]['offset_y_pa']['values'])*1000)

    # blade.sparcapwidth_hp = np.array(blade_internal_structure['layers'][spar_hp_index]['width']['values'])*1000
    # blade.sparcapwidth_lp  = np.array(blade_internal_structure['layers'][spar_lp_index]['width']['values'])*1000

    # blade.sparcapoffset_hp = np.array(blade_internal_structure['layers'][spar_hp_index]['offset_y_pa']['values'])*1000
    # blade.sparcapoffset_lp = np.array(blade_internal_structure['layers'][spar_lp_index]['offset_y_pa']['values'])*1000
    
    # TE Bands
    teReinfKeys=fullKeysFromSubStrings(bladeStructureDict.keys(),['te','reinf'])
    if len(teReinfKeys)==1:
        blade.teband = bladeStructureDict[teReinfKeys[0]]['width']['values']*1000 / 2
    elif len(teReinfKeys)==2:
        blade.teband = (bladeStructureDict[teReinfKeys[0]]['width']['values']+bladeStructureDict[teReinfKeys[1]]['width']['values'])*1000 / 2
    else:
        raise ValueError('Unknown number of TE reinforcements')

    # LE Bands
    leReinfKeys=fullKeysFromSubStrings(bladeStructureDict.keys(),['le','reinf'])
    if len(leReinfKeys)==1:
        blade.leband = bladeStructureDict[leReinfKeys[0]]['width']['values']*1000 / 2
    elif len(leReinfKeys)==2:
        blade.leband = (bladeStructureDict[leReinfKeys[0]]['width']['values']+bladeStructureDict[leReinfKeys[1]]['width']['values'])*1000 / 2
    else:
        raise ValueError('Invalid number of LE reinforcements')
   
    
    ### COMPONENTS
    _add_components(blade, blade_internal_structure, bladeStructureDict)
    
    blade.updateBlade()
    # save(blade_name)
    # BladeDef_to_NuMADfile(obj,numad_name,matdb_name,numad_af_folder)
    return blade


def _add_stations(blade,blade_outer_shape_bem, hub_outer_shape_bem,
                    af_data, file: str, write_airfoils):

    # Obtaining some parameters not explicitly given in YAML file
    L = np.ceil(blade_outer_shape_bem['reference_axis']['z']['values'][-1])
    R = L + hub_outer_shape_bem['diameter'] / 2
    L = R - hub_outer_shape_bem['diameter'] / 2
    blade.ispan = np.multiply(np.transpose(blade_outer_shape_bem['chord']['grid']),L)
    
    
    #Aerodynamic properties
    # using interp because yaml can have different r/R for twist and chord
    temp_x = np.transpose(blade_outer_shape_bem['twist']['grid'])
    temp_y = blade_outer_shape_bem['twist']['values']
    blade.degreestwist = interpolator_wrap(np.multiply(temp_x,L),np.transpose(temp_y),blade.ispan) * 180.0 / np.pi
    blade.chord = interpolator_wrap(
        np.multiply(np.transpose(blade_outer_shape_bem['chord']['grid']),L),
        np.transpose(blade_outer_shape_bem['chord']['values']),blade.ispan)
    af_dir_names = []
    for i in range(len(af_data)):
        af_dir_names.append(af_data[i]['name'])
    numstations = len(blade_outer_shape_bem['airfoil_position']['labels'])
    tc = [None]*numstations
    aero_cent = [None]*numstations

    for i in range(numstations):
        _,_,iaf_temp = np.intersect1d(blade_outer_shape_bem['airfoil_position']['labels'][i],af_dir_names,'stable',return_indices=True)
        IAF = iaf_temp[0] # Expect only one index of intersection
        tc[i] = af_data[IAF]['relative_thickness']
        tc_xL = blade_outer_shape_bem['airfoil_position']['grid'][i]
        aero_cent[i] = af_data[IAF]['aerodynamic_center']
        x = np.array(af_data[IAF]['coordinates']['x'], dtype=float)
        y = np.array(af_data[IAF]['coordinates']['y'], dtype=float)
        xf_coords = np.stack((x,y),1)

        # find coordinate direction (clockwise or counter-clockwise) Winding
        # Number. clockwise starting at (1,0) is correct
        with np.errstate(divide='ignore', invalid='ignore'):
            if np.nanmean(np.gradient(np.arctan(xf_coords[:,1] / xf_coords[:,0]))) > 0:
                xf_coords = np.flipud(xf_coords)

        if write_airfoils:
            import os
            out_folder = 'yaml2BladeDef_' + file.replace('.yaml','')
            # blade_name = out_folder + '/' + file.replace('.yaml','') + '_blade.mat'
            # matdb_name =...
            # numade_name =...

            # Creating folders
            os.makedirs(out_folder+'/af_coords/', exist_ok = True)
            # os.makedirs(out_folder+'/af_polars/', exist_ok = True)
            os.makedirs(out_folder+'/airfoil/', exist_ok = True)
            writeNuMADAirfoil(xf_coords, 
                blade_outer_shape_bem['airfoil_position']['labels'][i],
                out_folder + '/af_coords/' +
                blade_outer_shape_bem['airfoil_position']['labels'][i]+'.txt')

        ref = blade_outer_shape_bem['airfoil_position']['labels'][i]
        af = Airfoil(coords = xf_coords, ref = ref)
        af.resample(spacing='half-cosine')
        blade.addStation(af,tc_xL*L)
    # Obtain some key blade attributes
    blade.span = blade.ispan
    blade.percentthick = np.multiply(interpolator_wrap(np.multiply(blade_outer_shape_bem['airfoil_position']['grid'],L),tc,blade.ispan),100)
    blade.aerocenter = interpolator_wrap(np.multiply(blade_outer_shape_bem['airfoil_position']['grid'],L),aero_cent,blade.span)
    blade.chordoffset = interpolator_wrap(np.multiply(np.transpose(blade_outer_shape_bem['pitch_axis']['grid']),L),
        np.transpose(blade_outer_shape_bem['pitch_axis']['values']),blade.span)
    blade.naturaloffset = 0
    blade.prebend = interpolator_wrap(np.multiply(np.transpose(blade_outer_shape_bem['reference_axis']['x']['grid']),L),
        np.transpose(blade_outer_shape_bem['reference_axis']['x']['values']),blade.span)
    blade.sweep = interpolator_wrap(np.multiply(np.transpose(blade_outer_shape_bem['reference_axis']['y']['grid']),L),
        np.transpose(blade_outer_shape_bem['reference_axis']['y']['values']),blade.span)

    # for i in range(len(tc)):
    #     afc = AirfoilDef(out_folder + 
    #     '/af_coords/' + 
    #     blade_outer_shape_bem['airfoil_position']['labels'][i] +
    #     '.txt')
    #     blade.addStation(afc,np.multiply(tc_xL[i],L))

    #NOTE nothing happens to afc? Tentatively ignoring...
    # If i return to this make sure to listify the afcs
    ### AIRFOILS
    # for i in range(len(tc)):
    #     afc = AirfoilDef(out_folder + '/af_coords/' + 
    #         blade_outer_shape_bem['airfoil_position']['labels'][i] + 
    #         '.txt')
    #     blade.addStation(afc,np.multiply(tc_xL[i],L))
    # afc.resample #NOTE afc isn't used after this... why resample?
    return


def _add_materials(blade, material_data):
    materials_dict =dict()
    for i in range(len(material_data)):
        cur_mat = Material()
        cur_mat.name = material_data[i]['name']
        if material_data[i]['orth'] == 1:
            cur_mat.type = 'orthotropic'
        else:
            cur_mat.type = 'isotropic'
        # Add ply thickness option if ply thickness exists in yaml
        try:
            cur_mat.layerthickness = material_data[i]['ply_t'] * 1000
        except KeyError:
                print('Warning! material ply thickness ' + 
                        material_data[i]['name'] + 
                        ' not defined, assuming 1 mm thickness')
                cur_mat.layerthickness = 1
            
        finally:
            pass

        # first
        cur_mat.uts = _parse_data(material_data[i]['Xt'])
        cur_mat.ucs = -_parse_data(material_data[i]['Xc'])
        cur_mat.uss = _parse_data(material_data[i]['S'])
        cur_mat.xzit = 0.3
        cur_mat.xzic = 0.25
        cur_mat.yzit = 0.3
        cur_mat.yzic = 0.25
        try: 
            cur_mat.g1g2 = material_data[i].get('GIc',0) / material_data[i].get('GIIc',0)
        except ZeroDivisionError:
            cur_mat.g1g2 = np.nan
        if 'alp0' in material_data[i]:
            cur_mat.alp0 = _parse_data(material_data[i]['alp0'])
            cur_mat.etat = LARCetaT(cur_mat.alp0)
        else: 
            cur_mat.alp0 = None
            cur_mat.etat = None
        try:
            #test if property is a list
            material_data[i]['E']+[]
        except TypeError:
            cur_mat.ex = _parse_data(material_data[i]['E'])
            cur_mat.ey = _parse_data(material_data[i]['E'])
            cur_mat.ez = _parse_data(material_data[i]['E'])
            cur_mat.gxy = _parse_data(material_data[i]['G'])
            cur_mat.gxz = _parse_data(material_data[i]['G'])
            cur_mat.gyz = _parse_data(material_data[i]['G'])
            cur_mat.prxy = _parse_data(material_data[i]['nu'])
            cur_mat.prxz = _parse_data(material_data[i]['nu'])
            cur_mat.pryz = _parse_data(material_data[i]['nu'])
            cur_mat.etal = LARCetaL(cur_mat.uss,cur_mat.ucs,cur_mat.alp0)
        else:
            cur_mat.ex = _parse_data(material_data[i]['E'][0])
            cur_mat.ey = _parse_data(material_data[i]['E'][1])
            cur_mat.ez = _parse_data(material_data[i]['E'][2])
            cur_mat.gxy = _parse_data(material_data[i]['G'][0])
            cur_mat.gxz = _parse_data(material_data[i]['G'][1])
            cur_mat.gyz = _parse_data(material_data[i]['G'][2])
            cur_mat.prxy = _parse_data(material_data[i]['nu'][0])
            cur_mat.prxz = _parse_data(material_data[i]['nu'][1])
            cur_mat.pryz = _parse_data(material_data[i]['nu'][2])
            cur_mat.etal = LARCetaL(cur_mat.uss[0],cur_mat.ucs[1],cur_mat.alp0)
        try:
            cur_mat.m = material_data[i]['m']
        except KeyError:
            print(f"No fatigue exponent found for material: {material_data[i]['name']}")
        cur_mat.density = material_data[i]['rho']
        # cur_mat.dens = mat_data[i]['rho']
        cur_mat.drydensity = material_data[i]['rho']
        if 'description' in material_data[i].keys() and 'source' in material_data[i].keys():
            desc_sourc = [material_data[i]['description'],', ',material_data[i]['source']]
            cur_mat.reference = ''.join(desc_sourc)
        else:
            cur_mat.reference = []

        materials_dict[cur_mat.name] = cur_mat
    blade.materials = materials_dict
    return


def _add_components(blade, blade_internal_structure, bladeStructureDict):
    N_layer_comp = len(blade_internal_structure['layers'])
    component_list = list()
    for i in range(N_layer_comp):
        i_component_data = blade_internal_structure['layers'][i]
        cur_comp = Component()
        cur_comp.group = 0
        cur_comp.name = i_component_data['name']
        #   comp['material'] = blade_internal_structure['layers']{i}['material'];
        # mat_names = [mat.name for mat in blade.materials]
        # C,IA,IB = np.intersect1d(mat_names,i_component_data['material'],return_indices=True)
        cur_comp.materialid = i_component_data['material']
        try:
            cur_comp.fabricangle = np.mean(i_component_data['fiber_orientation']['values'])
        finally:
            pass
        if 'spar' in i_component_data['name'].lower():
            cur_comp.imethod = 'pchip'
        else:
            cur_comp.imethod = 'linear'
        # cur_comp.cp[:,0] = np.transpose(i_component_data['thickness']['grid'])
        cptemp1 = np.transpose(i_component_data['thickness']['grid'])
        temp_n_layer = np.multiply(np.transpose(i_component_data['thickness']['values']),1000.0) / blade.materials[cur_comp.materialid].layerthickness
        I_round_up = np.flatnonzero((temp_n_layer > 0.05) & (temp_n_layer < 0.5))
        cptemp2 = np.round(np.multiply(np.transpose(i_component_data['thickness']['values']),1000.0) / blade.materials[cur_comp.materialid].layerthickness)
        cur_comp.cp = np.stack((cptemp1,cptemp2),axis=1)
        # if I_round_up.size > 0:
        #     cur_comp.cp[I_round_up,1] = 1 # increase n_layers from 0 to 1 for 0.05<n_layers<0.5
        #     comp['cp'](:,2) = cell2mat(blade_internal_structure['layers']{i}['thickness']['values'])'.*1000;  # use when each material ply is 1 mm
        cur_comp.pinnedends = 0
        component_list.append(cur_comp)

    component_dict = dict()
    for comp in component_list:
        component_dict[comp.name] = comp

    # Spar Caps (pressure and suction)
    keyList=fullKeysFromSubStrings(component_dict.keys(),['spar','ps'])
    component_dict[keyList[0]].hpextents = ['b','c']

    keyList=fullKeysFromSubStrings(component_dict.keys(),['spar','ss'])
    component_dict[keyList[0]].lpextents = ['b','c']


        # uv coating
    keyList=fullKeysFromSubStrings(component_dict.keys(),['uv'])  #Try 1
    if len(keyList)==0:
        keyList=fullKeysFromSubStrings(component_dict.keys(),['gel']) #Try 2

    if len(keyList)==1:
        component_dict[keyList[0]].hpextents = ['le','te']
        component_dict[keyList[0]].lpextents = ['le','te']
    elif len(keyList)==0:
        raise ValueError('No UV or gelcoat found')
    else:
        raise ValueError('Too many uv or gelcoat components')

    # Shell skin
    keyList=fullKeysFromSubStrings(component_dict.keys(),['shell']) 
    if len(keyList)==2:
        component_dict[keyList[0]].hpextents = ['le','te']
        component_dict[keyList[1]].lpextents = ['le','te']
    else:
        raise ValueError('Incorrect number of shell components')

    # TE Band(s)
    keyList=fullKeysFromSubStrings(component_dict.keys(),['te','reinf'])
    if len(keyList)==1:
        component_dict[keyList[0]].hpextents = ['d','te']
        component_dict[keyList[0]].lpextents = ['d','te'] 
    elif len(keyList)==2:
        tempKeyList=fullKeysFromSubStrings(keyList,['ss'])
        if len(tempKeyList)==1:
            component_dict[tempKeyList[0]].lpextents = ['d','te'] 
        else:
            ValueError('Incorrect number of te reinf ss components')

        tempKeyList=fullKeysFromSubStrings(keyList,['ps'])
        if len(tempKeyList)==1:
            component_dict[tempKeyList[0]].hpextents = ['d','te'] 
        else:
            ValueError('Incorrect number of te reinf ps components')      
    else:
        raise ValueError('Invalid number of LE reinforcements')


    # LE Band(s)
    keyList=fullKeysFromSubStrings(component_dict.keys(),['le','reinf'])
    if len(keyList)==1:
        component_dict[keyList[0]].hpextents = ['le','a']
        component_dict[keyList[0]].lpextents = ['le','a'] 
    elif len(keyList)==2:
        tempKeyList=fullKeysFromSubStrings(keyList,['ss'])
        if len(tempKeyList)==1:
            component_dict[tempKeyList[0]].lpextents = ['le','a'] 
        else:
            ValueError('Incorrect number of te reinf ss components')

        tempKeyList=fullKeysFromSubStrings(keyList,['ps'])
        if len(tempKeyList)==1:
            component_dict[tempKeyList[0]].hpextents = ['le','a'] 
        else:
            ValueError('Incorrect number of te reinf ps components')      
    else:
        raise ValueError('Invalid number of LE reinforcements')
    
    # Trailing edge suction-side panel
    keyList=fullKeysFromSubStrings(component_dict.keys(),['te_ss'])
    if len(keyList)==1:
        component_dict[keyList[0]].lpextents = ['c','d']
    else:
        raise ValueError('Invalid number of trailing edge suction-side panels')

    # Leading edge suction-side panel
    keyList=fullKeysFromSubStrings(component_dict.keys(),['le_ss'])
    if len(keyList)==1:
        component_dict[keyList[0]].lpextents = ['a','b']
    else:
        raise ValueError('Invalid number of leading edge suction-side panels')  


    # Trailing edge suction-side panel
    keyList=fullKeysFromSubStrings(component_dict.keys(),['le_ps'])
    if len(keyList)==1:
        component_dict[keyList[0]].hpextents = ['a','b']
    else:
        raise ValueError('Invalid number of leading edge pressure-side panels')

    # Leading edge suction-side panel
    keyList=fullKeysFromSubStrings(component_dict.keys(),['te_ps'])
    if len(keyList)==1:
        component_dict[keyList[0]].hpextents = ['c','d']
    else:
        raise ValueError('Invalid number of trailing edge pressure-side panels')  
    
    #Web

    for comp in component_dict:
        print(comp)
    keyList=fullKeysFromSubStrings(component_dict.keys(),['web','fore']) #Try 1
    if len(keyList)==0:
        keyList=fullKeysFromSubStrings(component_dict.keys(),['web','1']) #Try 2

    if len(keyList)>0:
        for key in keyList:
            component_dict[key].hpextents = ['b']
            component_dict[key].lpextents = ['b']
            component_dict[key].group = 1
    elif len(keyList)==0:
        raise ValueError('No fore web layers found found') 


    keyList=fullKeysFromSubStrings(component_dict.keys(),['web','aft']) #Try 1
    if len(keyList)==0:
        keyList=fullKeysFromSubStrings(component_dict.keys(),['web','0']) #Try 2
    if len(keyList)==0:
        keyList=fullKeysFromSubStrings(component_dict.keys(),['web','rear']) #Try 3

    if len(keyList)>0:
        for key in keyList:
            component_dict[key].hpextents = ['c']
            component_dict[key].lpextents = ['c']
            component_dict[key].group = 2
    elif len(keyList)==0:
        raise ValueError('No rear web layers found found')    

    

    ### add components to blade
    blade.components = component_dict
    return


def writeNuMADAirfoil(coords, reftext, fname): 
    """
    WriteNuMADAirfoil  Write NuMAD airfoil files
    **********************************************************************
    *                   Part of the SNL NuMAD Toolbox                    *
    * Developed by Sandia National Laboratories Wind Energy Technologies *
    *             See license.txt for disclaimer information             *
    **********************************************************************
      WriteNuMADAirfoil(coords,reftext,fname)
        
            fname - full filename, incl extension, of NuMAD airfoil file to write
        coords - Nx2 array of airfoil coordinate data.  First column contains
        x-values, second column contains y-values.  Airfoil coordinates are in
        order as specified by NuMAD (i.e. trailing edge = (1,0) and leading
        edge = (0,0)
        reftext = string representing reference text
    """
    with open(fname,'wt') as fid:
        fid.write('<reference>\n%s</reference>\n' % (reftext))
        fid.write('<coords>\n' % ())
        for i in range(coords.shape[0]):
            fid.write('%8.12f\t%8.12f\n' % tuple(coords[i,:]))
        
        fid.write('</coords>' % ())