from mantid.simpleapi import (CreateSampleWorkspace,
                              SetSample,
                              SetGoniometer,
                              mtd)

class SampleModel():

    def __init__(self):

        CreateSampleWorkspace(OutputWorkspace='sample')

    def get_shape_dict(self, shape, params, alpha=0, beta=0, gamma=0):

        if shape == 'Sphere':
            radius = params[0]/200
            shape = ' \
            <sphere id="sphere"> \
              <radius val="{}" /> \
              <centre x="0.0" y="0.0" z="0.0" /> \
              <rotate x="{}" y="{}" z="{}" /> \
            </sphere> \
            '.format(radius,alpha,beta,gamma)
        elif shape == 'Cylinder':
            radius, height = params[0]/200, params[1]/100
            shape = ' \
            <cylinder id="cylinder"> \
              <centre z="0.0" y="0.0" p="0.0" /> \
              <axis x="0.0" y="1.0" z="0" /> \
              <radius val="{}" /> \
              <height val="{}" /> \
              <rotate x="{}" y="{}" z="{}" /> \
            </cylinder> \
            '.format(radius,height,alpha,beta,gamma)
        else:
            width, height, depth = params[0]/100, params[1]/100, params[1]/100
            shape = ' \
            <cuboid id="cuboid"> \
              <width val="{}" /> \
              <height val="{}" /> \
              <depth val="{}" /> \
              <centre x="0.0" y="0.0" z="0.0" /> \
              <rotate x="{}" y="{}" z="{}" /> \
            </cuboid> \
            '.format(width,height,depth,alpha,beta,gamma)

            return {'Shape': 'CSG', 'Value': shape}

    def get_material_dict(self, chemical_formula, z_parameter, volume):

        mat_dict = {'ChemicalFormula': chemical_formula,
                    'ZParameter': z_parameter,
                    'UnitCellVolume': volume}

        return mat_dict

    def set_sample(self, shape_dict, mat_dict, omega=0, chi=0, phi=0):

        SetGoniometer(Workspace='sample',
                      Axis0='{},0,1,0,1'.format(omega),
                      Axis1='{},0,0,1,1'.format(chi),
                      Axis2='{},0,1,0,1'.format(phi))

        SetSample(InputWorkspace='sample',
                  Geometry=shape_dict,
                  Material=mat_dict)

    def get_absorption_dict(self):

        mat = mtd['sample'].sample().getMaterial()

        sigma_a = mat.absorbXSection()
        sigma_s = mat.totalScatterXSection()

        M = mat.relativeMolecularMass()
        n = mat.numberDensityEffective 
        N = mat.totalAtoms

        V = abs(mtd['sample'].sample().getShape().volume()*100**3)

        rho = (n/N)/0.6022*M
        m = rho*V

        mu_s = n*sigma_s
        mu_a = n*sigma_a

        abs_dict = {'sigma_a': sigma_a, # barn
                    'sigma_s': sigma_s, # barn
                    'mu_a': mu_a, # 1/cm
                    'mu_s': mu_s, # 1/cm
                    'N': N, # atoms
                    'M': M, # g/mol
                    'n': n, # 1/A^3
                    'rho': rho, # g/cm^3
                    'V': V, # cm^3 
                    'm': m} # g

        return abs_dict