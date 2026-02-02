"""
Vehicle Braking Simulation Utilities

This module provides core physics and math utilities for simulating
vehicle braking behavior and stopping distance. It includes functions
for vector cross products and a forward-time numerical simulation of
vehicle deceleration under braking forces.

The braking simulation models:
- Distribution of braking torque between front and rear wheels
- Static and dynamic friction limits
- Longitudinal deceleration due to braking forces
- Torque-induced load transfer between axles
- Forward-flip detection when rear normal force becomes negative

All calculations are performed using basic Newtonian mechanics and
simple numerical integration with a fixed time step.

Functions
---------
cross_prod(A, B)
    Compute the 3D cross product of two vectors (with zero-padding for
    lower-dimensional inputs).

simulate_stop(F_percent, Total_Torque, Wheelbase, COG, Wheel_diameter,
              Initial_Velocity, COF, DCOF, Mass, step)
    Simulate vehicle braking and return the stopping distance.

Notes
-----
- This module does not depend on external physics engines or NumPy.
- Numerical accuracy depends on the chosen integration time step.
- The braking model assumes planar motion and rigid-body behavior.
- A sentinel stopping distance is returned if a forward flip is detected.
"""
import math

def cross_prod(A,B):
    """
    Compute the 3D cross product of two vectors.

    This function computes A × B for vectors of dimension 1, 2, or 3.
    Vectors with fewer than 3 components are zero-padded to 3D.
    Inputs with dimension greater than 3 are not allowed.

    Parameters
    ----------
    A : list or sequence of float
        First input vector (length ≤ 3).
    B : list or sequence of float
        Second input vector (length ≤ 3).

    Returns
    -------
    list of float
        A 3-element list representing the cross product A × B.

    Raises
    ------
    Exception
        If either input vector has dimension greater than 3.
    """
    if max(len(A), len(B)) > 3:
        raise Exception("cross_prod is only defined for up to and including 3-dimensional inputs")

    a = [i for i in A]    # Define a variable a = A to avoid function side-effects
    b = [i for i in B]    # Define a variable b = B to avoid function side-effects

    while len(a) < 3: a += [0]   # Append [0] to a until a is a 3-dimensional vector
    
    while len(b) < 3: b += [0]   # Append [0] to b until b is a 3-dimensional vector
    
    ab0 = a[1]*b[2] - a[2]*b[1]  # Define the x dimension of A cross B
    ab1 = a[2]*b[0] - a[0]*b[2]  # Define the y dimension of A cross B
    ab2 = a[0]*b[1] - a[1]*b[0]  # Define the z dimension of A cross B

    return [ab0, ab1, ab2]


def simulate_stop(F_percent, Total_Torque, Wheelbase, COG, Wheel_diameter, Initial_Velocity, COF, DCOF, Mass, step):
    """
    Simulate the stopping distance of a vehicle under braking.

    This function performs a forward-time numerical simulation of
    a vehicle decelerating due to braking forces applied at the
    front and rear wheels. Braking forces are limited by static
    and dynamic friction and produce a torque about the vehicle's
    center of gravity, which is used to update normal forces.

    The simulation proceeds until the vehicle velocity reaches zero
    or the rear normal force becomes negative (indicating a forward flip).

    Parameters
    ----------
    F_percent : float
        Percentage of total braking torque applied to the front wheels.
    Total_Torque : float
        Total braking torque applied to the vehicle (Nm).
    Wheelbase : float
        Distance between front and rear wheels (m).
    COG : list or tuple of float
        Center of gravity position [x, y] relative to the rear axle (m).
    Wheel_diameter : float
        Diameter of the wheels (m).
    Initial_Velocity : float
        Initial vehicle velocity (m/s).
    COF : float
        Static coefficient of friction.
    DCOF : float
        Dynamic coefficient of friction.
    Mass : float
        Vehicle mass (kg).
    step : float
        Time step used for numerical integration (s).

    Returns
    -------
    float
        Total stopping distance (m). If rear normal force becomes
        negative, a large sentinel value is returned.
    """

    # Use the F_percent parameter and the Total_Torque paramater to calculate front and rear stopping torques
    Front_Wheel_stop_torque = (F_percent / 100) * Total_Torque
    Rear_Wheel_stop_torque = ((100 - F_percent) / 100) * Total_Torque

    G = 9.81                  # Constant acceleration due to gravity (m/s^2)
    FN = G * Mass             # Defining the initial normal force of the entire car (N)

    rFx = Wheelbase - COG[0]  # Defining the position vector from the COG to the front wheel (x dimension)(m)
    rFy = 0 - COG[1]          # Defining the position vector from the COG to the front wheel (y dimension)(m)
    rF = [rFx, rFy]           # Defining the position vector from the COG to the front wheel (m)

    rRx = COG[0]              # Defining the position vector from the COG to the rear wheel (x dimension)(m)
    rRy = 0 - COG[1]          # Defining the position vector from the COG to the rear wheel (y dimension)(m)
    rR = [rRx, rRy]           # Defining the position vector from the COG to the rear wheel (m)

    FNF0 = FN * (rFx/(rFx + rRx))  # Defining the front wheel initial normal force modulus value (N)
    FNR0 = FN * (rRx/(rRx + rFx))  # Defining the rear wheel initial normal force modulus value (N)

    FNF = FNF0    # Setting front normal force modulus value (variable) to be the initial value (N)
    FNR = FNR0    # Setting rear normal force modulus value (variable) to be the initial value (N)

    FFapp = Rear_Wheel_stop_torque / (Wheel_diameter  / 2) # Calculating braking force applied by the front wheel (N)
    FRapp = Front_Wheel_stop_torque / (Wheel_diameter / 2) # Calculating braking force applied by the rear wheel (N)

    velo = Initial_Velocity  # Setting the velo (variable) to initially be the initial velocity value (m/s)
    FT = 0                   # Setting the moment on the vehicle to be 0 initially (Nm)

    time = 0       # Starting a clock for passing time as the simulation proceeds (s)
    count = 0      # Setting a count for the number of iterations the simualtion has taken
    distance = 0   # Starting a count for the distance travelled while stopping (m)

    while velo > 0:   # While the kart is still in motion

        count += 1    # Iterate the count

        # Calculate the braking force of each of the front and rear wheel
            # The braking force is taken to be the minimum of the applied braking force, and
            # the static coefficient of friction multiplied by the modulus of the normal force
            # for each wheel.
        FBF = [(-1) * min(COF * FNF, FFapp), 0]   # [N, N]
        FBR = [(-1) * min(COF * FNR, FRapp), 0]   # [N, N]

        # For each wheel, if the static coefficient of friction was taken to be the minimum in
        # the previous block, then we now apply the dynamic coefficient of friction to calculate the
        # braking force.
        if COF * FNF < FFapp:
            FBF = [(-1) * DCOF * FNF, 0]  # [N, N]
        if COF * FNR < FRapp:
            FBR = [(-1) * DCOF * FNR, 0]  # [N, N]

        FB = FBF[0] + FBR[0]  # Total braking force is the sum of the front and rear braking forces (N)


        FTF = cross_prod(rF, FBF)     # The moment applied to the COM of the vehicle by the front braking force (Nm)
        FTR = cross_prod(rR, FBR)     # The moment applied to the COM of the vehicle by the rear braking force (Nm)


        ft = [FTF[i] + FTR[i] for i in range(min(len(FTF), len(FTR)))]  # Intermedite vector sum of front and rear moments (Nm)
        FT = abs(ft[2])   # Modulus of the moment on the COM of the vehicle (Nm)


        A = FB/Mass    # Calculating the x-dimension of acceleration on the vehicle (m/s^2)
        velo += step * A   # Taking an approximate integral with dt = step for the velocity value (m/s)
        distance += step * velo  # Taking an approximate integral with dt = step for the distance value (m)

        FNF = FNF0 + FT / math.sqrt(rFx**2 + rFy**2)  # Recalculating the front normal force with the torque on the vehicle included (N)
        FNR = FNR0 - FT / math.sqrt(rFx**2 + rFy**2)  # Recalculating the rear normal force with the torque on the vehicle included (N)
        time += step            # Incrememnting the time clock

        # Adding a flag value and loop break in the event that the normal force
        # on the rear wheel is less than 0 (car is flipped forward)
        if FNR < 0:
            distance = 10000
            break

    return distance
