"""
Braking Simulation Analysis and Visualization Utilities

This module provides helper functions for running parameter-sweep
braking simulations and visualizing the resulting stopping-distance
data. It acts as a bridge between the low-level physics simulation
(`simulate_stopping`) and higher-level plotting routines (`plots`).

Core functionality includes:
- Basic statistical utilities (mean, standard deviation)
- Matrix manipulation helpers (transpose)
- Automated 3D parameter sweeps over exactly two independent variables
- 3D visualization of stopping-distance surfaces
- Reduction of 3D simulation data into 2D plots via aggregation
  (minimum, mean, or standard deviation)

The primary workflow is:
1. Specify exactly two simulation inputs as `range` objects.
2. Run `get_3d_data()` to compute stopping distances over the parameter grid.
3. Visualize results in 3D using a perceptually uniform colormap.
4. Optionally reduce the 3D surface to a 2D curve using
   `translate_3d_2d()` for further analysis.

Notes
-----
- Exactly two simulation inputs must be provided as `range` objects.
- All other inputs are treated as constants during the sweep.
- Data spacing is assumed to be uniform along both independent axes.
- Plotting is handled immediately via Matplotlib; figures are not returned.
- Statistical functions return 0 for insufficient data as a deliberate
  design choice, differing from NumPy conventions.
"""
import simulate_stopping as ss
import plots
import math

def clear_terminal():
    """
    Clear the terminal screen using an ANSI escape sequence.

    This function prints the ANSI escape code '\\033[2J', which clears
    the terminal display on terminals that support ANSI control codes.
    Cursor position is not reset.

    Notes
    -----
    - This may not work in all environments (e.g., some IDE consoles).
    - No value is returned.
    """
    print('\033[2J')  # Print ANSI Escape code to clear terminal

def mean(data):
    """
    Compute the arithmetic mean of a list of numeric values.

    Parameters
    ----------
    data : list of float or int
        Input data values.

    Returns
    -------
    float
        Arithmetic mean of the input values. Returns 0 if the input
        list is empty.

    Notes
    -----
    - This function does not perform type checking.
    - Returning 0 for empty input is a design choice and differs from
      NumPy's behavior (which returns NaN).
    """
    if len(data) == 0:  # Making system robust to 0-length list inputs
        return 0
    return sum(data) / len(data)  # Sum formula

def transpose(data):
    """
    Compute the transpose of a 2D data matrix.

    Parameters
    ----------
    data : list of lists
        A 2D list representing a matrix. Rows are allowed to have
        different lengths; missing values are filled with zeros.

    Returns
    -------
    list of lists
        Transposed matrix where rows and columns are swapped.

    Notes
    -----
    - If rows have unequal lengths, missing values are padded with zeros.
    - This function does not modify the input data.
    """
    nrow = len(data)  # Define the number of rows in the data-matrix
    ncol = max([len(row) for row in data])  # Define the number of columns in the data matrix

    output_matrix = [[0 for j in range(ncol)] for i in range(nrow)]  # Setting an empty list to hold the output

    for i in range(len(data)):   # Numbering each row i
        for j in range(len(data[i])):  # Numbering each datapoint (in row i) j
            output_matrix[j][i] = data[i][j]  # Setting the transpose of each ij-entry (placed at ji)
    
    return output_matrix

def standard_deviation(data):
    """
    Compute the sample standard deviation of a dataset.

    This function calculates the Bessel-corrected (sample) standard
    deviation using the formula:

        sqrt( sum((x - mean)^2) / (n - 1) )

    Parameters
    ----------
    data : list of float or int
        Input data values.

    Returns
    -------
    float
        Sample standard deviation of the data. Returns 0 if the
        dataset contains fewer than two elements.

    Notes
    -----
    - Uses n - 1 in the denominator (sample standard deviation).
    - Returning 0 for n < 2 is a design choice and differs from NumPy,
      which returns NaN.
    """
    n = len(data)   # Set n to be the number of datapoints
    if n < 2:       # If there is insufficient data for a SD, return SD = 0
        return 0
    mean_x = mean(data)   # Deterine the mean of the data
    SD = math.sqrt((sum([(x - mean_x)**2 for x in data])) / (n - 1))  # Calculate standard deviation

    return SD

def get_3d_data(F_percent, Total_Torque, Wheelbase, COG, Wheel_diameter, Initial_Velocity, COF, DCOF, Mass, step, updates = False, return_data = False, colormap = None, switch_order = False):
    """
    Generate and plot 3D stopping-distance data from a braking simulation.

    This function performs a parameter sweep over exactly two independent
    input variables (specified as Python `range` objects) while holding
    all other inputs constant. For each pair of independent-variable values,
    a braking simulation is executed and the resulting stopping distance
    is recorded. The collected data are plotted as a 3D scatter plot.

    Parameters
    ----------
    F_percent : float or range
        Percentage of total braking torque applied to the front axle.
        May be specified as a constant or a range for parameter sweeping.
    Total_Torque : float or range
        Total applied stopping torque (Nm).
    Wheelbase : float or range
        Vehicle wheelbase length (m).
    COG : list or tuple or range
        Center-of-gravity location, formatted as [x, y] in meters,
        or a range if used as an independent variable.
    Wheel_diameter : float or range
        Wheel diameter (m).
    Initial_Velocity : float or range
        Initial vehicle velocity (m/s).
    COF : float or range
        Static coefficient of friction.
    DCOF : float or range
        Dynamic coefficient of friction.
    Mass : float or range
        Vehicle mass (kg).
    step : float or range
        Simulation time step (s).

    updates : bool, optional
        If True, prints progress updates to the terminal during simulation.
        Default is False.
    return_data : bool, optional
        If True, returns the computed stopping-distance data as a 2D list.
        Default is False.
    colormap : str, optional
        Matplotlib perceptually uniform sequential colormap name.
        Accepted values are 'viridis', 'plasma', 'inferno', 'magma',
        and 'cividis'. Defaults to 'plasma' if invalid or None.
    switch_order : bool, optional
        If True, swaps the order of the two independent variables on the
        X and Y axes. Default is False.

    Returns
    -------
    datas : list of lists, optional
        2D list of stopping distances (meters), where rows correspond
        to the first independent variable and columns correspond to
        the second. Only returned if `return_data=True`.

    Raises
    ------
    Exception
        If the number of independent variables (range inputs) is not
        exactly two.

    Notes
    -----
    - Exactly two inputs must be specified as `range` objects.
    - All other inputs are treated as constants during the sweep.
    - This function calls `ss.simulate_stop()` to compute stopping distance.
    - The resulting data are visualized using `plots.graph3d()`.
    """
    #######################################################################
    ### DETERMINING THE X AND Y VARIABLE INPUTS FOR DATA-PLOTTING SETUP ###
    #######################################################################
    # Defining input names to automate the axis labels
    input_names = ['% Front Stop-Torque',
                   'Total Stopping Torque (Nm)',
                   'Wheelbase Length (m)',
                   'COG ([m,m])',
                   'Wheel Diameter (m)',
                   'Initial Velocity (m/s)',
                   'Static Coefficient of Friction',
                   'Dynamic Coefficient of Friction',
                   'Mass (kg)',
                   'Step (s)']
    
    # Setting a list of all the simulation inputs
    input_list = [F_percent, Total_Torque, Wheelbase, 
                  COG, Wheel_diameter, Initial_Velocity, 
                  COF, DCOF, Mass, step]

    positions = []  # Setting an empty list to record the list positions of the independent variables
    ranges = []     # Setting an empty list to record the ranges of the independent variables

    for n, i in enumerate(input_list):  # For each input and its position value;
        if type(i) == type(range(0)):       # If it is a range input;
            ranges += [i]                       # Append its range to 'ranges'
            positions += [n]                    # Append its list position to 'positions'
    
    # Raising an error if more or less than 2 variable inputs were given
    if len(positions) != 2:
        raise Exception("get_3d_data only accepts 2 independent variables.")
    
    # Simple boolean to allow a second range input to be x-axis
    if switch_order:
        positions = [positions[1], positions[0]]  # Reverse positions order
        ranges = [ranges[1], ranges[0]]           # Reverse ranges order
    
    x_axis_title = input_names[positions[0]]  # Set x-axis title to be the X-position value of the axis title list
    y_axis_title = input_names[positions[1]]  # Set y-axis title to be the Y-position value of the axis title list
    
    ####################################
    ### COLLECTING DATA FOR PLOTTING ###
    ####################################
    M = max(ranges[0]) / (ranges[0][1] - ranges[0][0])  # Set a maximum value for the x-axis (for simulation update)
    count = 0  # Setting an iteration count (for simulation update)
    
    datas = []  # Setting an empty list to collect stopping distance data rows
    for iv1 in ranges[0]:   # For each x-value (from the x-range)

        if updates:         # If the user wants printed simulation updates
            clear_terminal()  # Clear the terminal
            print(f'\nSimulation ~{round(100 * count / M, 2)}% complete.\n')  # Print a simulation update

        data = []  # Set an empty list to collect this row's stopping distances
        for iv2 in ranges[1]:   # For each y-value (from the y-range)

            inputs = []        # Create an empty list to collect simulation inputs
            p1_used = False    # Create a simple switch to determine whether x-value is used yet
            for i in range(10):  # 10 = number of simulation inputs

                if i in positions:   # If this simulation input is one of the independent variables
                    if p1_used:      # If the x-value is already collected
                        inputs += [iv2]   # Append the y-value to the inputs list
                    else:            # Otherwise
                        inputs += [iv1]   # Append the x-value to the inputs list
                        p1_used = True    # Switch the x-value boolean switch to True

                else:                # If this simulation input is one of the constant inputs
                    inputs += [input_list[i]]    # Append the constant value to the inputs list

            # Find a stop distance for this particular set of inputs
            this_data = ss.simulate_stop(inputs[0], inputs[1], inputs[2],
                                         inputs[3], inputs[4], inputs[5],
                                         inputs[6], inputs[7], inputs[8],
                                         inputs[9])
            
            data += [this_data]  # Append this stop distance to this row's data list
        

        datas += [data]  # Append this row of data to the data matrix
        count += 1       # Increment the count by 1
    
    ###################################################
    ### DETERMINING ADDITIONAL PLOTTING INFORMATION ###
    ###################################################

    minx = min(ranges[0])   # Find minimum x-value
    countx = ranges[0][1] - ranges[0][0]  # Determine the count-rate of the X-variable
    miny = min(ranges[1])   # Find minimum y-value
    county = ranges[1][1] - ranges[1][0]  # Determine the count-rate of the Y-variable


    if colormap in ['viridis', 'plasma', 'inferno', 'magma', 'cividis']:  # Checking that colormap input is allowed
        cmap = colormap   # If it is, set cmap plotting input to be colormap input
    else:
        cmap = 'plasma'   # Otherwse, set the colormap to 'plasma'

    ############################################
    ### PRINTING 3D GRAPH AND RETURNING DATA ###
    ############################################
    plots.graph3d(datas, minx, countx,   # Plot data with plots.py
                  miny, county, cmap,
                  x_axis_title, y_axis_title,
                  'Stopping Distance (m)')

    if return_data:    # If user wants data output, return the collected data
        return datas
    else:
        data = None    # Ensure garbage collection from Python (helps with running many simulations)
    

def translate_3d_2d(get_3d_data_output, x_min, x_count, xaxisname, yaxisname, stay_axis = 'x', translation = 'mean'):
    """
    Reduce 3D simulation data to a 2D plot by collapsing one axis.

    This function converts a 2D matrix of simulation output values
    (representing a 3D surface) into a 1D dataset by applying a
    reduction operation (min, mean, or standard deviation) along
    one axis. The resulting data is plotted as a 2D graph.

    Parameters
    ----------
    get_3d_data_output : list of lists of float
        2D matrix of simulation output values (Z data).
    x_min : float
        Minimum value of the independent axis.
    x_count : float
        Step size between consecutive values on the independent axis.
    xaxisname : str
        Label for the x-axis.
    yaxisname : str
        Label for the y-axis.
    stay_axis : {'x', 'y'}, optional
        Axis to preserve as the independent variable.
        - 'x' reduces across rows
        - 'y' reduces across columns
    translation : {'min', 'mean', 'sd'}, optional
        Reduction operation applied to collapse the data.

    Returns
    -------
    None
        This function produces a plot but does not return data.

    Notes
    -----
    - If an invalid translation is provided, the mean is used by default.
    - Standard deviation uses a sample (n - 1) denominator.
    - Data spacing is assumed to be uniform.
    """
    if translation in ['min', 'mean', 'sd']:  # If the translation type is defined
        # Set a variable function 'func' to be the associated translation function
        if translation == 'min':
            func = min
        elif translation == 'mean':
            func = mean
        else:
            func = standard_deviation
    # Otherwise set translation type to be mean, and print an error code
    else:
        func = mean
        print(f"Translation type {translation} is not defined.")

    if stay_axis in ['x', 'y']:  # If stay_axis input is an allowed input
        if stay_axis == 'x':     # If the axis to stay is the x-axis
            datas = []               # Set an empty list to collect 2D data
            for row in get_3d_data_output:  # For each row in the matrix output
                datas += [func(row)]            # Append the function of the row to the datas list
        else:                    # If the axis to stay is the y-axis
            # Transpose the data (switch x and y)
            new_data = transpose(get_3d_data_output)
            datas = []               # Set an empty list to collect 2D data
            for row in new_data:     # For each row in the matrix output
                datas += [func(row)]      # Append the function of the row to the datas list
        
        # Plot the new 2D data in the X-Y plane
        plots.graph2d(datas, x_min, x_count, xaxisname, yaxisname)
    
    # Otherwise raise an error for improper stay-axis input
    else:
        raise Exception("translate_3d_2d requires [stay_axis] paramater be \'x\' or \'y\'.")



