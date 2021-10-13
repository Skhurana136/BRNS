# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

%reset

from ogs5py import OGS
from ogs5py.reader import readtec_point
from ogs5py.reader import readtec_polyline
from matplotlib import pyplot as plt
import os

os.chdir('C:/Users/sanew/Desktop/uni/MScEnviroFood/4Sem/MA/OGSBRNS/Implements')

model = OGS(task_root='TRACER_V2', task_id='1d_tr',output_dir='out')

# generate a radial mesh
#model.msh.generate("rectangular", dim=2)
# generate a radial outer boundary
#model.gli.generate("rectangular", dim=2, x=0.5, y=0.5)
#model.gli.add_points([0., 0., 0.], "pwell")
#model.gli.add_points([1., 0., 0.], "owell")

# Work around for mesh generation
model.gli.read_file('C:/Users/sanew/Desktop/uni/MScEnviroFood/4Sem/MA/OGSBRNS/Implements/TRACER_V2/1d_tr.gli')

model.msh.read_file('C:/Users/sanew/Desktop/uni/MScEnviroFood/4Sem/MA/OGSBRNS/Implements/TRACER_V2/1d_tr.msh')

# Define boundary conditions: ONLY REQUIRED FOR MOBILE SPECIES!
model.bc.add_block(  # boundary condition
    PCS_TYPE='RICHARDS_FLOW',
    PRIMARY_VARIABLE='PRESSURE1',
    GEO_TYPE=['POINT', "POINT0"],
    DIS_TYPE=['CONSTANT', -9000.0],
)
model.bc.add_block(  # boundary condition
    PCS_TYPE='MASS_TRANSPORT',
    PRIMARY_VARIABLE='tracer',
    GEO_TYPE=['POINT', "POINT1"],
    DIS_TYPE=['CONSTANT', 2],
)

# Define Source and Sink Terms
model.st.add_block(  # source term
    PCS_TYPE='RICHARDS_FLOW',
    PRIMARY_VARIABLE='PRESSURE1',
    GEO_TYPE=['POINT', "POINT1"],
    DIS_TYPE=['CONSTANT', 0.015],
)

# Define Initial Conditions

#model.ic.read_file('C:/Users/sanew/Desktop/uni/MScEnviroFood/4Sem/MA/OGSBRNS/Implements/1D_FULL_DORM_V2/1d_fd.ic')

model.ic.add_block(  # initial condition
    PCS_TYPE='RICHARDS_FLOW',
    PRIMARY_VARIABLE='PRESSURE1',
    GEO_TYPE='DOMAIN',
    DIS_TYPE=['CONSTANT', -15000.0],
)
model.ic.add_block(  # initial condition
    PCS_TYPE='MASS_TRANSPORT',
    PRIMARY_VARIABLE='tracer',
    GEO_TYPE='DOMAIN',
    DIS_TYPE=['CONSTANT', 1.0e-20],
)

# Processes
model.pcs.add_block(  # set the process type
    PCS_TYPE='RICHARDS_FLOW',
    NUM_TYPE='NEW',
)
model.pcs.add_block(  # set the process type
    PCS_TYPE='MASS_TRANSPORT',
    NUM_TYPE='NEW',
)

# Component Properties
model.mcp.add_block(
    NAME='tracer',
    MOBILE=1,
    DIFFUSION=[1, 1.0e-09],
)

# rfd (CURVES) file
model.rfd.add_block(
    CURVE=[[0, 1],
           [10, 1],
           [11, 0],
           [100, 0],
            ],
)
model.rfd.add_block(
    CURVE=[[0, 1],
           [10, 1],
           [11, 0],
           [100, 0],
            ],
)

# Define Medium Properties
model.mmp.add_block(  # medium properties
    GEOMETRY_DIMENSION=1,
    GEOMETRY_AREA=1.0,
    # STORAGE=[1, 1.0e-04],
    PERMEABILITY_TENSOR=['ISOTROPIC', 1.20e-3],
    PERMEABILITY_SATURATION=[4, 0.1548, 0.5226, 0.7814],
    CAPILLARY_PRESSURE=[4, 0.02076],
    POROSITY=[1, 0.39],
    TORTUOSITY=[1, 1],
    MASS_DISPERSION=[1, 0.001, 0.1],
    DENSITY=[1, 1.2e03],
)

# Define Solid Properties
model.msp.add_block(
    DENSITY=[1, 1.2e03],        
)

# Define Fluid Properties
model.mfp.add_block(
    FLUID_TYPE='WATER',
    PCS_TYPE='PRESSURE1',
    DENSITY=[1, 1000],
    VISCOSITY=[1, 1.0e-03],
)

# Numerical Solver
model.num.add_block(  # numerical solver
    OVERALL_COUPLING=[1, 12],
    PCS_TYPE='RICHARDS_FLOW',
    LINEAR_SOLVER=[3, 6, 1.0e-14, 1000, 1, 100, 4],
    NON_LINEAR_ITERATION=['PICARD', 'ERNORM', 1000, 0, 1e-07],
    COUPLING_CONTROL=['ERNORM', 5e-03],
)
model.num.add_block(  # numerical solver
    PCS_TYPE='MASS_TRANSPORT',
    LINEAR_SOLVER=[2, 6, 1.0e-14, 5000, 1, 100, 4],
    COUPLING_CONTROL=['ERNORM', 5e-03],
)

# Output
model.out.add_block(  # point observation
    PCS_TYPE=['RICHARDS_FLOW'],
    NOD_VALUES=[['PRESSURE1'],
                ['SATURATION1']],
    GEO_TYPE='DOMAIN',
    DAT_TYPE='TECPLOT',
    TIM_TYPE=['STEPS', 10],
)
model.out.add_block(  # point observation
    PCS_TYPE=['MASS_TRANSPORT'],
    NOD_VALUES=['tracer'],
    GEO_TYPE='DOMAIN',
    DAT_TYPE='TECPLOT',
    TIM_TYPE=['STEPS', 10],
)

# Time settings
model.tim.add_block(  # set the timesteps
    PCS_TYPE='RICHARDS_FLOW',
    TIME_START=0.0,
    TIME_END=100,
    TIME_CONTROL=[['SELF_ADAPTIVE'],
                  [3, 1.3],
                  [7, 0.8],
                  ['MIN_TIME_STEP'],
                  [1.0e-03],
                  ['MAX_TIME_STEP'],
                  [1],
                  ['INITIAL_STEP_SIZE'],
                  [1],
                  ['ITERATIVE_TYPE'],
                  ['COUPLED'],
                  ['STAY'],
                  [3],
            ]
)
model.tim.add_block(  # set the timesteps
    PCS_TYPE='MASS_TRANSPORT',
    TIME_START=0.0,
    TIME_END=100,
    TIME_CONTROL=[['SELF_ADAPTIVE'],
                  [3, 1.3],
                  [7, 0.8],
                  ['MIN_TIME_STEP'],
                  [1.0e-03],
                  ['MAX_TIME_STEP'],
                  [1],
                  ['INITIAL_STEP_SIZE'],
                  [1],
                  ['ITERATIVE_TYPE'],
                  ['COUPLED'],
                  ['STAY'],
                  [3],
            ]
)

# Build Model
model.write_input()

# OGS Stand alone
#success = model.run_model(ogs_root='C:/Users/sanew/Desktop/uni/MScEnviroFood/4Sem/MA/OGS/ogs5/build/bin/Release')

# OGS-BRNS --- update folder
success = model.run_model(ogs_root='C:/Users/sanew/Desktop/uni/MScEnviroFood/4Sem/MA/OGSBRNS/Implements/1D_FULL_DORM_V2/ogs.exe')

point = readtec_polyline(
    task_root="1D_FULL_DORM_V2",
    task_id="1d_fd",
    pcs=None,
)

time = point['']["X"]
head = point['owell']["HEAD"]

#plt.plot(time, head)
#plt.show()