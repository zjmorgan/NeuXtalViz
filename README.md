# NeuXtalViz
Neutron Scattering Single Crystal Visualization

*new-crystal-vis*

NeXtalViz is an application for three-dimensional (3d) visualization of neutron scattering data.
It brings together two main libraries; PyVista is the main tool for displaying the 3d data while Mantid serves as the main library for working with reduced single crystal neutron diffraction data.

- [PyVista](https://pyvista.org/)
- [Mantid](https://github.com/mantidproject/mantid/)

### Getting started

Create conda environment
`conda env create -f environment.yml`

Activate garnet environment
`conda activate nxv`

Install in editable mode for developlment
`python -m pip install -e .`

Run the GUI
`python src/nxv.py`

### Additional information

The application is designed with a model-view-presenter pattern. This makes it possible to replace the view with a different one for other applications platforms. 
This could include [trame](https://kitware.github.io/trame/) for web-based applications.

[![](https://mermaid.ink/img/pako:eNptUM0KgzAMfpWSk4Ied5Gx044T3AY79RJsnIW2So0TEd999WcXWQ4h-f4ImaBsFEEGlWmGskbP4vaQToRqPXXkmHwUFb8xjg-cOKfpRdiQYaJ86btgRXYSHWv11_fRNESv0HbXsm9MMd75dATZoyVIwJK3qFW4elokErgmSxKyMCqqsDcsQbo5SLHn5jm6EjL2PSXQtwqZrhrfIQuyCk0XUFKaG59vn1gfMn8BnIFa3w?type=png)](https://mermaid.live/edit#pako:eNptUM0KgzAMfpWSk4Ied5Gx044T3AY79RJsnIW2So0TEd999WcXWQ4h-f4ImaBsFEEGlWmGskbP4vaQToRqPXXkmHwUFb8xjg-cOKfpRdiQYaJ86btgRXYSHWv11_fRNESv0HbXsm9MMd75dATZoyVIwJK3qFW4elokErgmSxKyMCqqsDcsQbo5SLHn5jm6EjL2PSXQtwqZrhrfIQuyCk0XUFKaG59vn1gfMn8BnIFa3w)
