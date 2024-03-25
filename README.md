# NeuXtalViz
Neutron Scattering Single Crystal Visualization

*new-crystal-vis*

NeuXtalViz is an application for three-dimensional (3d) visualization of neutron scattering data.
It brings together two main libraries; PyVista is the main tool for displaying the 3d data while Mantid serves as the main library for working with reduced single crystal neutron diffraction data.

- [PyVista](https://pyvista.org/)
- [Mantid](https://github.com/mantidproject/mantid/)

The application also relises on several other packages.

- [Matplotlib](https://matplotlib.org/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [scikit-learn](https://scikit-learn.org/stable/)

## Description

The current features of NeuXtalViz including the following tools

- Crystal Structure: load a CIF, display and modify crystal information.
- Sample Tools: load a UB, define a sample shape, and calculate absorption information.
- Modulation Information: load an indexed peaks workspace and use DBSCAN to determine modulation offsets.
- Volume Slicer: load a NXS histogram file and display slices from volume data.

[![](https://mermaid.ink/img/pako:eNp1Uk1vwjAM_StRTmGC_QAOkxCcJm2HlaFJaw8hcWlGmlT52ECI_z6nKR-qRA-26_ec2C8-UWEl0Dmttf0TDXeBrFelIfjxrmPsHeJX4Hqj_GRyTZPZ7IUId_SIsGX2pAguihAdjHiet50GVvSOLLbeui4oa0a0X6tjC2zTO1JoJcCNKK2VUfNUy96uIVlbqweij9ud411z64V8QK0MtGBCpghVfz89o63uhxiwYZIEKOORh7bKGEZ5mgb0gbEiuYsiYOTo_qV1DkTf3mdQWgUFPlPiFo9teajutBmKs0Ap3ew1stBWjy5YqbqOHnUSPARwyuwywRxS22irO1EfHTKICJIsDNdHr4YeO-D7dEzvq5H0mXL77-FgO4UFyQ0tpxChGfnhhjP2ivZOLjqlLbiWK4mbd0rpkoYGn6mkcwwl1DzqUNLSnJHKY7DF0Qg6x2eFKY2dxK5XiuMY7SUJUgXr3vIy9zt9_geLEe9f?type=png)](https://mermaid.live/edit#pako:eNp1Uk1vwjAM_StRTmGC_QAOkxCcJm2HlaFJaw8hcWlGmlT52ECI_z6nKR-qRA-26_ec2C8-UWEl0Dmttf0TDXeBrFelIfjxrmPsHeJX4Hqj_GRyTZPZ7IUId_SIsGX2pAguihAdjHiet50GVvSOLLbeui4oa0a0X6tjC2zTO1JoJcCNKK2VUfNUy96uIVlbqweij9ud411z64V8QK0MtGBCpghVfz89o63uhxiwYZIEKOORh7bKGEZ5mgb0gbEiuYsiYOTo_qV1DkTf3mdQWgUFPlPiFo9teajutBmKs0Ap3ew1stBWjy5YqbqOHnUSPARwyuwywRxS22irO1EfHTKICJIsDNdHr4YeO-D7dEzvq5H0mXL77-FgO4UFyQ0tpxChGfnhhjP2ivZOLjqlLbiWK4mbd0rpkoYGn6mkcwwl1DzqUNLSnJHKY7DF0Qg6x2eFKY2dxK5XiuMY7SUJUgXr3vIy9zt9_geLEe9f)

### Screenshots

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/1e44c60b-a9bf-45f2-902a-96bfc98355c5)

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/757ea243-56fd-4c2d-928c-c02099a3e824)

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/fd86ad60-797a-4574-9b7b-75bbddf8bd35)

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/7be57ab6-e928-4ef3-b7bf-e9ae7fb4311c)


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
