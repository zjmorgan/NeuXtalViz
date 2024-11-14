from mantid.simpleapi import mtd

import numpy as np

def SaveMDToAscii(workspace, filename, exclude_integrated=True, format='%.6e'):
    """
    Save an MDHistoToWorkspace to an ASCII file (column format).

    workspace : str
      Name of workspace as string.
    filename : str
      Path to output file.
    exclude_integrated : bool, optional
      Exclude integrated dimensions with bin size of one. Default is `True`.
    format : str
      Column format.

    """

    ws = mtd[workspace]

    if ws.id() != 'MDHistoWorkspace':
        raise ValueError('The workspace is not an MDHistoToWorkspace')

    if exclude_integrated:
        dims = ws.getNonIntegratedDimensions()
    else:
        dims = [ws.getDimension(i) for i in range(ws.getNumDims())]

    dimarrays = [np.linspace(d.getMinimum()+(d.getX(1)-d.getX(0))*0.5, \
                             d.getMaximum()-(d.getX(1)-d.getX(0))*0.5, \
                             d.getNBins()) for d in dims]

    if len(dimarrays) > 1:
        newdimarrays = np.meshgrid(*dimarrays, indexing='ij')
    else:
        newdimarrays = dimarrays

    data = ws.getSignalArray()*1.
    err = np.sqrt(ws.getErrorSquaredArray())

    header = 'Intensity Error '+' '.join([d.getName() for d in dims])
    header += '\n shape: '+'x'.join([str(d.getNBins()) for d in dims])

    to_save = np.column_stack([data.flatten(), err.flatten()])

    for d in newdimarrays:
        to_save = np.c_[to_save, d.flatten()]

    np.savetxt(filename, to_save, fmt=format, header=header)