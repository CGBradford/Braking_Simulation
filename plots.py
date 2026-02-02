"""
Visualization Utilities for Simulation Data

This module provides helper functions for visualizing numerical
simulation output using Matplotlib. It includes tools for rendering
both 3D scatter plots and 2D line plots from regularly sampled data.

The 3D plotting utility interprets a 2D array of scalar values as
samples on a uniform X–Y grid and renders them as a 3D scatter plot,
with point colors determined by a nonlinear sigmoid mapping of the
Z-values. Perceptually uniform sequential colormaps are supported
to improve visual interpretability.

The 2D plotting utility generates a standard X–Y line plot from a
one-dimensional dataset using a specified starting value and step
size for the independent variable.

Functions
---------
graph3d(data, xmin, xcount, ymin, ycount, color='plasma',
        Xtitle=None, Ytitle=None, Ztitle=None)
    Create a 3D scatter plot from a 2D data array sampled on a regular
    grid.

graph2d(data, x_min, x_count, xaxisname=None, yaxisname=None)
    Create a 2D line plot from a sequence of data values.

Notes
-----
- All plots are displayed immediately using Matplotlib and are not
  returned as objects.
- Color scaling in 3D plots is nonlinear due to the sigmoid mapping
  applied to Z-values.
- Only Matplotlib perceptually uniform sequential colormaps are
  supported for 3D plotting.
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def graph3d(data, xmin, xcount, ymin, ycount, color = 'plasma', Xtitle = None, Ytitle = None, Ztitle = None):
    """
    Create a 3D scatter plot from a 2D data array.

    This function interprets a 2D array of scalar values as Z-data
    sampled on a regular X–Y grid. X and Y coordinates are generated
    using the provided minimum values and step sizes. Each data point
    is plotted as a 3D scatter point.

    Color values are assigned using a modified sigmoid function applied
    to the Z-values, and mapped through a Matplotlib perceptually uniform
    sequential colormap.

    Parameters
    ----------
    data : list of lists or 2D array-like
        2D array of scalar values to be plotted as Z coordinates.
        Each element data[i][j] corresponds to a point at
        (xmin + i * xcount, ymin + j * ycount).
    xmin : float
        Minimum X-axis value.
    xcount : float
        Step size between successive X values.
    ymin : float
        Minimum Y-axis value.
    ycount : float
        Step size between successive Y values.
    color : str, optional
        Name of a Matplotlib perceptually uniform sequential colormap.
        Accepted values are 'viridis', 'plasma', 'inferno', 'magma',
        and 'cividis'. Defaults to 'plasma'.
    Xtitle : str, optional
        Label for the X-axis.
    Ytitle : str, optional
        Label for the Y-axis.
    Ztitle : str, optional
        Label for the Z-axis.

    Returns
    -------
    None
        This function displays a 3D plot using Matplotlib and does not
        return any values.

    Notes
    -----
    - The color scaling is nonlinear due to the sigmoid mapping.
    - If an unsupported colormap is provided, the function defaults
      to 'plasma' and prints a warning message.
    """
    e = 2.718  # Defining e = Euler's Number
    
    # Setting empty lists to hold X, Y, and Z data (rearranging data)
    X = []
    Y = []
    Z = []

    # Simple for loop for collecting X, Y, and Z data
    for i in range(len(data)):         # Numbering each row i
        for j in range(len(data[i])):  # Numbering each element (in row i) j
            X += [xcount * i + xmin]   # Append row number to X
            Y += [ycount * j + ymin]   # Append column number to Y
            Z += [data[i][j]]          # Append ij-entry to Z
    
    M = max(Z)   # Setting variable to hold max Z value
    m = min(Z)   # Setting variable to hold min Z value

    b = (2 * e) / (M - (sum(Z) / len(Z)))   # Setting modifier for the Sigmoid squishification function
    c = (-1) * b * (sum(Z) / len(Z))        # Setting modifier for the Sigmoid squishification function

    # Using modified sigmoid squishification function to apply a color value to each Z element
    colors = np.array([(e**(b * i + c))/(e**(b * i + c) + 1) for i in Z])  # Setting as a numpy array

    X = np.array(X)  # Setting X as a numpy array
    Y = np.array(Y)  # Setting Y as a numpy array
    Z = np.array(Z)  # Setting Z as a numpy array
    

    fig = plt.figure()  # Defining the matplotlib figure for 3D graphing
    ax = fig.add_subplot(111, projection='3d')  # Creating subplot within matplotlib

    # Create colormap validation check with an autocorrect to 'plasma' if user-defined
    # colormap is not recognized
    if not color in ['viridis', 'plasma', 'inferno', 'magma', 'cividis']:  # If color is not in the predetermined list
        color = 'plasma'                                                   # Redefine as 'plasma'
        print(f"Only Matplotlib Perpetual Uniform Sequential colormaps are accepted.")  # Print error string


    ax.scatter(X, Y, Z, c = colors, cmap = 'plasma')  # Define a scatterplot in ax using X,Y,Z data and colormap

    # Defining axis labels dependent on optional parameters
    if not Xtitle is None:
        ax.set_xlabel(str(Xtitle))
    if not Ytitle is None:
        ax.set_ylabel(str(Ytitle))
    if not Ztitle is None:
        ax.set_zlabel(str(Ztitle))

    plt.show()  # Utilize Matplotlib graph printing software to print 3D graph

def graph2d(data, x_min, x_count, xaxisname = None, yaxisname = None):
    """
    Create a 2D line plot from a sequence of data.

    This function generates X-axis values starting at a specified
    minimum and incrementing by a fixed step size. The provided data
    values are plotted against these X coordinates using Matplotlib.

    Parameters
    ----------
    data : list or 1D array-like
        Sequence of Y-values to be plotted.
    x_min : float
        Minimum value of the X-axis.
    x_count : float
        Step size between successive X values.
    xaxisname : str, optional
        Label for the X-axis.
    yaxisname : str, optional
        Label for the Y-axis.

    Returns
    -------
    None
        This function displays a 2D plot using Matplotlib and does not
        return any values.
    """
    X = [x_min + x_count * i for i in range(len(data))]  # Setting X data (for accurate numberline on x-axis) as a list
    Y = data                                             # Setting Y data as a list

    fig = plt.plot(np.array(X), np.array(Y))  # Defining a plot of numpy arrays for X and Y using Matplotlib

    # Defining axis labels dependent on optional parameters
    if not xaxisname is None:
        plt.xlabel(xaxisname)
    if not yaxisname is None:
        plt.ylabel(yaxisname)

    plt.show()  # Utilize Matplotlib graph printing software to print 2D graph
