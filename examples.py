"""
Interactive Braking Simulation Examples

This script provides a simple command-line interface for running
predefined braking simulation scenarios using the `get_data` module.
Each example performs a 3D parameter sweep over two independent
variables, visualizes the resulting stopping-distance surface, and
then reduces the 3D data into a 2D plot using a specified aggregation
method.

The user selects from several predefined examples that explore
different parameter combinations, such as:
- Front braking force distribution vs. total stopping torque
- Wheelbase length vs. stopping torque
- Effects of aggregation choice (mean vs. minimum stopping distance)

Workflow
--------
1. User selects an example number (1–4) from the terminal.
2. A 3D braking simulation is executed via `get_3d_data()`.
3. The resulting surface is plotted using a perceptually uniform colormap.
4. The 3D data are collapsed into a 2D curve using `translate_3d_2d()`.
5. The process repeats until the user exits.

Notes
-----
- Exactly two simulation inputs must be specified as `range` objects
  within each example definition.
- All other parameters are treated as constants during the sweep.
- This script is intended for exploratory analysis and visualization,
  not as a reusable library module.
- Figures are generated immediately and are not saved automatically.
"""
import get_data as GD

# Shows front force % on the x-axis and stopping torque on the y-axis, then keeps the front force percent on the x axis
# and uses mean value to collapse the y-axis
Example1 = [
    range(1,100,1),         # Front force percent
    range(50, 1000, 50),    # Total applied torque
    2,                      # Wheelbase
    [0.4, 0.3],             # Center of mass
    0.2,                    # Wheel diameter
    40,                     # Initial velocity
    0.8,                    # Static coefficient of friction
    0.2,                    # Dynamic coefficient of friction
    100,                    # Mass
    0.02,                   # Time step (dt)
    True,                   # Print simulation updates
    True,                   # Return simulation data
    'plasma',               # Colormap
    False,                  # Switch X and Y
    1,                      # Minimum X-value
    1,                      # X count
    'Front Force %',        # X-axis label
    'Braking Distance (m)', # Y-axis label
    'x',                    # Stay variable
    'mean'                  # 2D collapse type
]

# Shows front force % on the x-axis and stopping torque on the y-axis, then keeps the front force percent on the x axis
# and uses min value to collapse the y-axis
Example2 = [
    range(1,100,1),         # Front force percent
    range(50, 1000, 50),    # Total applied torque
    2,                      # Wheelbase
    [0.4, 0.3],             # Center of mass
    0.2,                    # Wheel diameter
    40,                     # Initial velocity
    0.8,                    # Static coefficient of friction
    0.2,                    # Dynamic coefficient of friction
    100,                    # Mass
    0.02,                   # Time step (dt)
    True,                   # Print simulation updates
    True,                   # Return simulation data
    'plasma',               # Colormap
    False,                  # Switch X and Y
    1,                      # Minimum X-value
    1,                      # X count
    'Front Force %',        # X-axis label
    'Braking Distance (m)', # Y-axis label
    'x',                    # Stay variable
    'min'                   # 2D collapse type
]

# Shows front force % on the x-axis and stopping torque on the y-axis, then keeps the front force percent on the x axis
# and uses mean value to collapse the y-axis
Example3 = [
    range(1,95,1),          # Front force percent
    range(20, 500, 10),     # Total applied torque
    2,                      # Wheelbase
    [0.4, 0.3],             # Center of mass
    0.2,                    # Wheel diameter
    40,                     # Initial velocity
    0.8,                    # Static coefficient of friction
    0.2,                    # Dynamic coefficient of friction
    100,                    # Mass
    0.01,                   # Time step (dt)
    True,                   # Print simulation updates
    True,                   # Return simulation data
    'plasma',               # Colormap
    False,                  # Switch X and Y
    1,                      # Minimum X-value
    1,                      # X count
    'Front Force %',        # X-axis label
    'Braking Distance (m)', # Y-axis label
    'x',                    # Stay variable
    'mean'                  # 2D collapse type
]

# Shows total stopping torque on the y-axis and wheelbase on the x-axis, then keeps the wheelbase on the x axis
# and uses mean value to collapse the y-axis
Example4 = [
    70,                     # Front force percent
    range(20, 1000, 20),    # Total applied torque
    range(1, 10, 1),        # Wheelbase
    [0.4, 0.3],             # Center of mass
    0.2,                    # Wheel diameter
    40,                     # Initial velocity
    0.8,                    # Static coefficient of friction
    0.2,                    # Dynamic coefficient of friction
    100,                    # Mass
    0.01,                   # Time step (dt)
    True,                   # Print simulation updates
    True,                   # Return simulation data
    'plasma',               # Colormap
    True,                   # Switch X and Y
    1,                      # Minimum X-value
    1,                      # X count
    'Wheelbase (m)',        # X-axis label
    'Braking Distance (m)', # Y-axis label
    'x',                    # Stay variable
    'mean'                  # 2D collapse type
]

while True:
    example_choice_input = input('Enter \"q\" to quit.\nChoose an example number 1-4 (inclusive)\n-> ')
    if example_choice_input in [str(i) for i in range(1, 5)]:
        example_choice = int(example_choice_input)
    elif example_choice_input == 'q':
        break
    if example_choice == 1:
        ed = Example1
    elif example_choice == 2:
        ed = Example2
    elif example_choice == 3:
        ed = Example3
    elif example_choice == 4:
        ed = Example4
    else: continue

    three = GD.get_3d_data(ed[0], ed[1], ed[2], ed[3], ed[4], ed[5], ed[6], ed[7], ed[8], ed[9], ed[10], ed[11], ed[12], ed[13])
    GD.translate_3d_2d(three, ed[14], ed[15], ed[16], ed[17], ed[18], ed[19])
