# Copyright 2022 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# See https://floris.readthedocs.io for documentation

import matplotlib.pyplot as plt
import numpy as np

from floris.tools import FlorisInterface
from floris.tools.optimization.yaw_optimization.yaw_optimizer_sr import YawOptimizationSR

import floris.tools.visualization as wakeviz
"""
This example demonstrates how to perform a yaw optimization for multiple wind directions
and 1 wind speed.

First, we initialize our Floris Interface, and then generate a 3 turbine wind farm.
Next, we create the yaw optimization object `yaw_opt` and perform the optimization using the
SerialRefine method. Finally, we plot the results.
"""
# Load the default example floris object
yaml_file = "inputs/cc.yaml" #"inputs/gch.yaml" #
# fi = FlorisInterface(yaml_file)  # GCH model matched to the default "legacy_gauss" of V2
fi = FlorisInterface("inputs/cc.yaml") # New CumulativeCurl model

# Reinitialize as a 2-turbine farm with range of WDs and 1 WS
D = 126.0  # Rotor diameter for the NREL 5 MW
ws = [8] # Wind speed
wd = [90] # Wind direction

nwty = 1 # y dir
mD = 5 # Distance turbines as multiplier rotor diameter

optPperWT = np.array(range(9))*1.1
for i in range(9):
    # Create the list using list comprehension
    nwt = i +2
    temp_x = [i * D * mD for i in range(nwt)]  # example [0.0, 5 * D, 10 * D, 15 * D, 20 * D]
    temp_y = [i * D * mD for i in range(nwty)]
    X, Y = np.meshgrid(temp_x, temp_y)
    layout_x = X.flatten()
    layout_y = Y.flatten()
    # layout_y = [0.0 for _ in layout_x]

    fi.reinitialize(
        layout_x=layout_x,
        layout_y=layout_y,
        wind_directions=wd,
        wind_speeds=ws,
    )
    # print(fi.floris.farm.rotor_diameters)

    # Initialize optimizer object and run optimization using the Serial-Refine method
    yaw_opt = YawOptimizationSR(fi)  # , exploit_layout_symmetry=False)
    df_opt = yaw_opt.optimize()

    # print("Optimization results:")
    # print(df_opt)

    # Split out the turbine results
    for t in range(nwt):
        df_opt['t%d' % t] = df_opt.yaw_angles_opt.apply(lambda x: x[t])

    farm_power_baseline = df_opt.farm_power_baseline[0]
    farm_power_baseline_str = "{:.2f}".format(farm_power_baseline / 10 ** 6)
    farm_power_opt = df_opt.farm_power_opt[0]
    farm_power_opt_str = "{:.2f}".format(farm_power_opt / 10 ** 6)
    power_increase_percent = (farm_power_opt - farm_power_baseline) / farm_power_baseline * 100
    power_increase_percent_str = "{:.2f}".format(power_increase_percent)
    print(power_increase_percent_str)
    optPperWT[i] = power_increase_percent

print(optPperWT)