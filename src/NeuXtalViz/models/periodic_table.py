from mantid.kernel import V3D

from mantid.geometry import (CrystalStructure,
                             ReflectionGenerator,
                             ReflectionConditionFilter,
                             PointGroup,
                             PointGroupFactory,
                             SpaceGroupFactory)

from mantid.simpleapi import (CreateSampleWorkspace,
                              LoadCIF,
                              mtd)

import numpy as np
import scipy.linalg

from NeuXtalViz.models.base_model import NeuXtalVizModel

class CrystalStructureModel(NeuXtalVizModel):

    def __init__(self, ref_ws=None):

        CreateSampleWorkspace(OutputWorkspace='crystal')

    def generate_space_groups_from_crystal_system(self, system):

        pg_system = getattr(PointGroup.CrystalSystem, system)
        pgs = list(PointGroupFactory.getPointGroupSymbols(pg_system))
        pgs = [PointGroupFactory.createPointGroup(pg) for pg in pgs]
        sgs = [SpaceGroupFactory.getSpaceGroupsForPointGroup(pg) for pg in pgs]
        sgs = [sg for sg_list in sgs for sg in sg_list]
        sgs = [SpaceGroupFactory.createSpaceGroup(sg) for sg in sgs]

        nos = np.unique([sg.getNumber() for sg in sgs]).tolist()

        space_group = []
        for no in nos:
            symbol = SpaceGroupFactory.subscribedSpaceGroupSymbols(no)[0]
            space_group.append('{}: {}'.format(no,symbol))

        return space_group

    def generate_settings_from_space_group(self, sg):

        no, symbol = sg.split(': ')

        return list(SpaceGroupFactory.subscribedSpaceGroupSymbols(int(no)))

    def load_CIF(self, filename):

        LoadCIF(Workspace='crystal', InputFile=filename)

        self.calculate_UB()

    def set_crystal_structure(self, params, space_group, scatterers):

        line = ' '.join(['{}']*6)

        constants = line.format(*params)

        atom_info = ';'.join([line.format(*s) for s in scatterers])

        cs = CrystalStructure(constants, space_group, atom_info)

        mtd['crystal'].sample().setCrystalStructure(cs)

        self.calculate_UB()

    def update_lattice_parameters(self, a, b, c, alpha, beta, gamma):

        scatterers = self.get_scatterers()

        sg = self.get_setting()

        params = [a, b, c, alpha, beta, gamma]

        self.set_crystal_structure(params, sg, scatterers)

    def calculate_UB(self):

        cs = mtd['crystal'].sample().getCrystalStructure()

        uc = cs.getUnitCell()

        G = uc.getG()
        G_star = np.linalg.inv(G)

        A = scipy.linalg.cholesky(G, lower=False)
        B = scipy.linalg.cholesky(G_star, lower=False)

        U = np.linalg.inv(A).T @ np.linalg.inv(B)

        self.set_UB(np.dot(U, B))

    def generate_F2(self, d_min=0.7):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        generator = ReflectionGenerator(cryst_struct)

        sf_filt = ReflectionConditionFilter.StructureFactor

        unit_cell = cryst_struct.getUnitCell()

        d_max = np.max([unit_cell.a(), unit_cell.b(), unit_cell.c()])

        hkls = generator.getUniqueHKLsUsingFilter(d_min, d_max, sf_filt)

        ds = generator.getDValues(hkls)

        F2s = generator.getFsSquared(hkls)

        return hkls, ds, F2s

    def calculate_F2(self, h, k, l):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        generator = ReflectionGenerator(cryst_struct)

        hkl = V3D(h, k, l)

        d = generator.getDValues([hkl])[0]

        F2 = generator.getFsSquared([hkl])[0]

        pg = cryst_struct.getSpaceGroup().getPointGroup()

        equivalents = pg.getEquivalents(hkl)

        return equivalents, d, F2

    def get_crystal_system(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        pg = cryst_struct.getSpaceGroup().getPointGroup()

        return pg.getCrystalSystem().name

    def get_lattice_system(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        pg = cryst_struct.getSpaceGroup().getPointGroup()

        return pg.getLatticeSystem().name

    def get_point_group_name(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        pg = cryst_struct.getSpaceGroup().getPointGroup()

        return pg.getName()

    def get_space_group(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        sg = cryst_struct.getSpaceGroup()

        no = sg.getNumber()
        symbol = SpaceGroupFactory.subscribedSpaceGroupSymbols(no)[0]

        return '{}: {}'.format(no,symbol)

    def get_setting(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        return cryst_struct.getSpaceGroup().getHMSymbol()

    def get_lattice_constants(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        uc = cryst_struct.getUnitCell()

        params = uc.a(), uc.b(), uc.c(), uc.alpha(), uc.beta(), uc.gamma()

        return params

    def get_unit_cell_volume(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()        

        return cryst_struct.getUnitCell().volume()

    def get_scatterers(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        scatterers = cryst_struct.getScatterers()

        scatterers = [atm.split(' ') for atm in list(scatterers)]

        scatterers = [[val if val.isalpha() else float(val) \
                       for val in scatterer] for scatterer in scatterers]

        return scatterers

    def get_chemical_formula_z_parameter(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        sg = cryst_struct.getSpaceGroup()

        scatterers = self.get_scatterers() 

        atom_dict = {}

        for scatterer in scatterers:
            atom, x, y, z, occ, Uiso = scatterer
            n = len(sg.getEquivalentPositions([x,y,z]))
            if atom_dict.get(atom) is None:
                atom_dict[atom] = [n], [occ]
            else:
                ns, occs = atom_dict[atom]
                ns.append(n)
                occs.append(occ)
                atom_dict[atom] = ns, occs

        chemical_formula = []

        n_atm = []
        n_wgt = []

        for key in atom_dict.keys():
            ns, occs = atom_dict[key]
            n_atm.append(np.sum(ns))
            n_wgt.append(np.sum(np.multiply(ns, occs)))
            if key.isalpha():
                chemical_formula.append(key+'{}')
            else:
                chemical_formula.append('('+key+')'+'{}')

        Z = np.gcd.reduce(n_atm)
        n = np.divide(n_wgt, Z)

        chemical_formula = '-'.join(chemical_formula).format(*n)

        return chemical_formula, Z

    def generate_atom_positions(self):

        scatterers = self.get_scatterers()

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        sg = cryst_struct.getSpaceGroup()

        A = self.get_unit_cell_transform()

        atom_dict = {}

        corners = [0,0,0], [1,0,0], [0,1,0], [0,0,1], \
                  [1,1,1], [0,1,1], [1,0,1], [1,1,0]

        for ind, scatterer in enumerate(scatterers):

            atom, x, y, z, occ, U = scatterer

            xyz = np.array(sg.getEquivalentPositions([x,y,z]))

            xyz = np.mod(xyz, 1)

            xyz = np.row_stack([xyz+corner for corner in corners])

            xyz = xyz[np.all(xyz <= 1, axis=1)]

            r_xyz = np.einsum('ij,kj->ki', A, xyz).tolist()
            r_occ = np.full(len(xyz), float(occ)).tolist()
            r_ind = np.full(len(xyz), ind).tolist()

            if atom_dict.get(atom) is None:
                atom_dict[atom] = r_xyz, r_occ, r_ind
            else:
                R_xyz, R_occ, R_ind = atom_dict[atom]
                R_xyz += r_xyz
                R_occ += r_occ
                R_ind += r_ind
                atom_dict[atom] = R_xyz, R_occ, R_ind

        return atom_dict

    def get_unit_cell_transform(self):

        cryst_struct = mtd['crystal'].sample().getCrystalStructure()

        uc = cryst_struct.getUnitCell()

        G = uc.getG()

        A = scipy.linalg.cholesky(G, lower=False)

        return A

    def constrain_parameters(self):

        params = np.array([False]*6)

        lattice_system = self.get_lattice_system()

        if lattice_system == 'Cubic':
            params[1:6] = True
        elif lattice_system == 'Rhombohedral':
            params[1:3] = True
            params[4:6] = True
        elif lattice_system == 'Hexagonal' or lattice_system == 'Tetragonal':
            params[2] = True
            params[3:6] = True
        elif lattice_system == 'Orthorhombic':
            params[3:6] = True
        elif lattice_system == 'Monoclinic':
            if 'unique axis b' in self.get_point_group_name():
                params[3] = True
                params[5] = True
            else:
                params[3:4] = True

        return params.tolist()

    def update_parameters(self, params):

        params = np.array(params)

        lattice_system = self.get_lattice_system()

        if lattice_system == 'Cubic':
            params[1:3] = params[0]
        elif lattice_system == 'Rhombohedral':
            params[1:3] = params[0]
            params[4:6] = params[3]
        elif lattice_system == 'Hexagonal' or lattice_system == 'Tetragonal':
            params[1] = params[0]

        return params.tolist()