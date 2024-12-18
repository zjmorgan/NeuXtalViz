# ![](https://github.com/zjmorgan/NeuXtalViz/blob/main/src/icons/neuxtalviz_logo.svg)

Neutron Scattering Single Crystal Visualization

*new-crystal-vis*

```git pull```

NeuXtalViz is an application for three-dimensional (3d) visualization of neutron scattering data.
It brings together two main libraries; PyVista is the main tool for displaying the 3d data while Mantid serves as the main library for working with reduced single crystal neutron diffraction data.

- [PyVista](https://pyvista.org/)
- [Mantid](https://github.com/mantidproject/mantid/)

The application also relises on several other packages.

- [Matplotlib](https://matplotlib.org/)
- [NumPy](https://numpy.org/)
- [SciPy](https://scipy.org/)
- [scikit-learn](https://scikit-learn.org/stable/)

- [PyVista](https://pyvista.org/)
  
## Description

The current features of NeuXtalViz including the following tools

- Crystal Structure: load a CIF, display and modify crystal information.
- Sample Tools: load a UB, define a sample shape, and calculate absorption information.
- Modulation Information: load an indexed peaks workspace and use DBSCAN to determine modulation offsets.
- Volume Slicer: load a NXS histogram file and display slices from volume data.

[![](https://mermaid.ink/img/pako:eNp1Uk1vwjAM_StRTmGC_QAOkxCcJm2HlaFJaw8hcWlGmlT52ECI_z6nKR-qRA-26_ec2C8-UWEl0Dmttf0TDXeBrFelIfjxrmPsHeJX4Hqj_GRyTZPZ7IUId_SIsGX2pAguihAdjHiet50GVvSOLLbeui4oa0a0X6tjC2zTO1JoJcCNKK2VUfNUy96uIVlbqweij9ud411z64V8QK0MtGBCpghVfz89o63uhxiwYZIEKOORh7bKGEZ5mgb0gbEiuYsiYOTo_qV1DkTf3mdQWgUFPlPiFo9teajutBmKs0Ap3ew1stBWjy5YqbqOHnUSPARwyuwywRxS22irO1EfHTKICJIsDNdHr4YeO-D7dEzvq5H0mXL77-FgO4UFyQ0tpxChGfnhhjP2ivZOLjqlLbiWK4mbd0rpkoYGn6mkcwwl1DzqUNLSnJHKY7DF0Qg6x2eFKY2dxK5XiuMY7SUJUgXr3vIy9zt9_geLEe9f?type=png)](https://mermaid.live/edit#pako:eNp1Uk1vwjAM_StRTmGC_QAOkxCcJm2HlaFJaw8hcWlGmlT52ECI_z6nKR-qRA-26_ec2C8-UWEl0Dmttf0TDXeBrFelIfjxrmPsHeJX4Hqj_GRyTZPZ7IUId_SIsGX2pAguihAdjHiet50GVvSOLLbeui4oa0a0X6tjC2zTO1JoJcCNKK2VUfNUy96uIVlbqweij9ud411z64V8QK0MtGBCpghVfz89o63uhxiwYZIEKOORh7bKGEZ5mgb0gbEiuYsiYOTo_qV1DkTf3mdQWgUFPlPiFo9teajutBmKs0Ap3ew1stBWjy5YqbqOHnUSPARwyuwywRxS22irO1EfHTKICJIsDNdHr4YeO-D7dEzvq5H0mXL77-FgO4UFyQ0tpxChGfnhhjP2ivZOLjqlLbiWK4mbd0rpkoYGn6mkcwwl1DzqUNLSnJHKY7DF0Qg6x2eFKY2dxK5XiuMY7SUJUgXr3vIy9zt9_geLEe9f)

### Screenshots

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/6ba6c8a7-0c6c-425c-9e00-5a34e69973b9)

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/3bef441d-db49-446f-b9f8-914e8eae6d18)

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/8db3f589-9127-4c2c-a6c5-a99660365258)

![image](https://github.com/zjmorgan/NeuXtalViz/assets/13754794/0a1a96fa-5d07-4494-8120-8d1bd80d9266)


### Getting started

Create conda environment
`conda env create -f environment.yml`

Activate garnet environment
`conda activate nxv`

Install in editable mode for developlment
`python -m pip install -e .`

Run the GUI
`python src/NeuXtalViz.py`

### Additional information

The application is designed with a model-view-presenter (MVP) pattern.
This makes it possible to replace the view with a different one for other applications platforms. 
This could include [trame](https://kitware.github.io/trame/) for web-based applications.

[![](https://mermaid.ink/img/pako:eNpdUMEKgzAM_ZWSk4L-gIydvAqOjZ16CTbOgq1S04mI_75W3Q7L4fHIey8kWaEZFEEBbT_MTYeOxaOUVoQyaFkrIfL8Kkww9UkVMT3UernxLr01zckzwCmwQ0M_5RwVg-ISm6OjiSyTS5L6S9MzGgN_LmkhA0POoFZhyzUaJXBHhiQUgSpq0fcsQdotWNHzcF9sAwU7Txn4USFTqfEV1oKixX4KXVKaB1cdl-8P2D5mNFW4?type=png)](https://mermaid.live/edit#pako:eNpdUMEKgzAM_ZWSk4L-gIydvAqOjZ16CTbOgq1S04mI_75W3Q7L4fHIey8kWaEZFEEBbT_MTYeOxaOUVoQyaFkrIfL8Kkww9UkVMT3UernxLr01zckzwCmwQ0M_5RwVg-ISm6OjiSyTS5L6S9MzGgN_LmkhA0POoFZhyzUaJXBHhiQUgSpq0fcsQdotWNHzcF9sAwU7Txn4USFTqfEV1oKixX4KXVKaB1cdl-8P2D5mNFW4)
