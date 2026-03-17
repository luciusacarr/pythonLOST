import pickle
import sys
import types
import numpy as np
from scipy.spatial.transform import Rotation

# 1. MOCK THE CUSTOM CLASS
# This prevents a ModuleNotFoundError if you don't have the 'kw_ebs_star_tracking' package installed.
module_name = 'kw_ebs_star_tracking.core.camera'
kw_mock = types.ModuleType(module_name)
sys.modules['kw_ebs_star_tracking'] = types.ModuleType('kw_ebs_star_tracking')
sys.modules['kw_ebs_star_tracking.core'] = types.ModuleType('kw_ebs_star_tracking.core')
sys.modules[module_name] = kw_mock

# Create a dummy class to catch the data
class Camera:
    pass

kw_mock.Camera = Camera

# 2. LOAD THE FILE
file_path = 'ast_camera.camera_ast' # Make sure this matches your file name
with open(file_path, 'rb') as f:
    camera_data = pickle.load(f)

# 3. EXTRACT TIMES AND ROTATIONS
times = camera_data.times
rotations = camera_data.rotations 

# 4. CALCULATE RA, DEC, AND ROLL
# Assumption: The camera's boresight (forward direction) is the +Z axis [0, 0, 1].
# If your specific camera software uses +X as forward, change this to [1, 0, 0].
boresight = np.array([0, 0, 1])

# Apply all the rotations to the boresight vector at once
pointing_vectors = rotations.apply(boresight)

# Extract the X, Y, and Z components of the new pointing vectors
x = pointing_vectors[:, 0]
y = pointing_vectors[:, 1]
z = pointing_vectors[:, 2]

# Calculate Right Ascension (RA) and Declination (Dec)
ra = np.degrees(np.arctan2(y, x))
ra = np.mod(ra, 360) # Normalize RA to be between 0 and 360 degrees
dec = np.degrees(np.arcsin(z))

# Calculate Roll by converting the 3D rotation into Euler angles.
# 'zyx' extracts the rotations around the Z, Y, and X axes. 
# The rotation around the boresight (Z) is your roll.
euler_angles = rotations.as_euler('zyx', degrees=True)
roll = euler_angles[:, 0] 

# 5. PRINT THE RESULTS
print(f"{'Time':<15} | {'RA (deg)':<12} | {'Dec (deg)':<12} | {'Roll (deg)':<12}")
print("-" * 58)

# Print the first 20 keyframes as an example
for i in range(min(20, len(times))):
    print(f"{times[i]:<15.4f} | {ra[i]:<12.4f} | {dec[i]:<12.4f} | {roll[i]:<12.4f}")