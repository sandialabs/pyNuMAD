# -*- coding: utf-8 -*-
"""
@author: evaande
"""

## Note:
## This script is meant to be run after the accompanying example, 'write_abaqus_shell_model.py',
## and after running the resulting modal analysis execution script in Abaqus FEA software.
## The blade model/mesh yaml file 'BAR0.yaml' and the Abaqus results yaml 'BAR0DampingInp.yaml'
## produced by those steps are needed inputs for the stuctural damping loss factor calculations
## performed below.

import pynumad as pynu
import yaml
from yaml import CLoader as Loader
from pynumad.analysis.abaqus.read import *
from pynumad.utils.damping import *
from pynumad.io.mesh_to_yaml import *

modelYaml = 'BAR0.yaml'
abqResYaml = 'BAR0DampingInp.yaml'

## Define material loss factors for every relevant material in the blade
mat_lf = {'Gelcoat': [1.,1.,1.,1.,1.,1.], 
          'Adhesive': [0.3,0.3,0.3,0.3,0.3,0.3],
          'glass_biax': [0.8,0.8,0.5,0.3,0.3,0.3],
          'glass_triax': [0.85,0.85,0.5,0.4,0.4,0.4],
          'glass_uni': [0.9,0.5,0.5,0.3,0.3,0.3],
          'medium_demsity_foam': [0.3,0.3,0.3,0.3,0.3,0.3]}

print('reading model yaml...')
## Read in model/mesh data

model_data = yaml_to_mesh(modelYaml)

print('constructing damping input...')

## Extract results from abaqus output
damp_inp = build_damping_input(abqResYaml, mat_lf, model_data)

print('getting structural damping factors...')

## Get the structural loss factors

struct_lf = get_modal_loss_factors(damp_inp)

print('structural damping factors:')
print(struct_lf)