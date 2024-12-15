beamlines = {
    'SNAP': {
        'Name': 'SNAP',
        'InstrumentName': 'SNAP',
        'Facility': 'SNS',
        'Wavelength': [1, 4],
        'Grouping': '4x4',
        'BankPixels': [256, 256],
        'Goniometers': ['BL3:Mot:omega,0,1,0,1'],
        'Goniometer': {
            'Goniometer': {
                'BL3:Mot:omega': [0, 1, 0, 1, 0, 360]
                }
            },
        'Motor': {
            'det_lin1': 0,
            'det_lin2': 0,
            'det_arc1': -75,
            'det_arc2': 105,
        },
        'RawFile': 'nexus/SNAP_{}.nxs.h5',
    },
    'CORELLI': {
        'Name': 'CORELLI',
        'InstrumentName': 'CORELLI',
        'Facility': 'SNS',
        'Wavelength': [0.6, 2.5],
        'Grouping': '1x4',
        'BankPixels': [16, 256],
        'Goniometers': [
            'BL9:Mot:Sample:Axis1,0,1,0,1',
            'BL9:Mot:Sample:Axis2,0,1,0,1',
            'BL9:Mot:Sample:Axis3,0,1,0,1',
        ],
        'Goniometer': {
            'Goniometer': {
                'BL9:Mot:Sample:Axis1': [0, 1, 0, 1, 0, 0],
                'BL9:Mot:Sample:Axis2': [0, 1, 0, 1, 0, 0],
                'BL9:Mot:Sample:Axis3': [0, 1, 0, 1, 0, 360]
             }
        },
        'RawFile': 'nexus/CORELLI_{}.nxs.h5',
    },
    'TOPAZ': {
        'Name': 'TOPAZ',
        'InstrumentName': 'TOPAZ',
        'Facility': 'SNS',
        'Wavelength': [0.4, 3.5],
        'Grouping': '4x4',
        'BankPixels': [256, 256],
        'Goniometers': [
            'BL12:Mot:omega,0,1,0,1',
            'BL12:Mot:chi,0,0,1,1',
            'BL12:Mot:phi,0,1,0,1',
        ],
        'Goniometer': {
            'Ambient': {
                'BL12:Mot:goniokm:omega': [0, 1, 0, 1, -180, 180],
                'BL12:Mot:goniokm:chi':  [0, 0, 1, 1, 135, 135],
                'BL12:Mot:goniokm:phi':  [0, 1, 0, 1, -180, 180]
                },
            'Cryogenic': {
                'BL12:Mot:Gonioc:Omega': [0, 1, 0, 1, -180, 180],
                'BL12:Mot:Gonioc:Chi':  [0, 0, 1, 1, 0, 0],
                'BL12:Mot:Gonioc:Phi':  [0, 1, 0, 1, 0, 0]
                }
             },
        'RawFile': 'nexus/TOPAZ_{}.nxs.h5',
    },
    'MANDI': {
        'Name': 'MANDI',
        'InstrumentName': 'MANDI',
        'Facility': 'SNS',
        'Wavelength': [2, 4],
        'Grouping': '4x4',
        'BankPixels': [256, 256],
        'Goniometers': [
            'BL11B:Mot:omega,0,1,0,1',
            'BL11B:Mot:chi,0,0,1,1',
            'BL11B:Mot:phi,0,1,0,1',
        ],
        'Goniometer': {
            'Goniometer': {
                'BL11B:Mot:omega': [0, 1, 0, 1, 0, 90],
                'BL11B:Mot:chi': [0, 0, 1, 1, 130, 130],
                'BL11B:Mot:phi': [0, 1, 0, 1, -180, 180]
            }
        },
        'RawFile': 'nexus/MANDI_{}.nxs.h5',
    },
    'WAND²': {
        'Name': 'WAND',
        'InstrumentName': 'HB2C',
        'Facility': 'HFIR',
        'Wavelength': 1.486,
        'Grouping': '4x4',
        'BankPixels': [480, 512],
        'Goniometers': ['s1,0,1,0,1'],
        'Goniometer': {
            'Goniometer': {
                'HB2C:Mot:sgl': [1, 0, 0, -1, 0, 0],
                'HB2C:Mot:sgu': [0, 0, 1, -1, 0, 0],
                'HB2C:Mot:s1': [0, 1, 0, 1, -180, 180]
            }
        },
        'Motor': {
            'HB2C:Mot:s2.RBV': 30,
            'HB2C:Mot:detz.RBV': 0,
        },
        'RawFile': 'nexus/HB2C_{}.nxs.h5',
    },
    'DEMAND': {
        'Name': 'HB3A',
        'InstrumentName': 'HB3A',
        'Facility': 'HFIR',
        'Wavelength': 1.546,
        'Grouping': '4x4',
        'BankPixels': [512, 512],
        'Goniometers': [
            'omega,0,1,0,-1',
            'chi,0,0,1,-1',
            'phi,0,1,0,-1',
        ],
        'Goniometer': {
            'Goniometer': {
                'omega': [0, 1, 0, -1, -13, 45],
                'chi': [0, 0, 1, -1, -91, 91],
                'phi': [0, 1, 0, -1, -180, 180]
            }
        },
        'Motor': {
            '2theta': 30,
            'det_trans': 410.38595,
        },
        'RawFile': 'shared/autoreduce/HB3A_exp{:04}_scan{:04}.nxs',
    },
}
