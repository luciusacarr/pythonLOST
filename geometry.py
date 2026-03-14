# geometry.py
import math
from dataclasses import dataclass

# --- Basic Math & Vectors ---
@dataclass
class Vec2:
    x: float
    y: float
    
    def __add__(self, other):
        return Vec2(self.x + other.x, self.y + other.y)
        
    def __sub__(self, other):
        return Vec2(self.x - other.x, self.y - other.y)
        
    def __mul__(self, scalar: float):
        return Vec2(self.x * scalar, self.y * scalar)
        
    def __rmul__(self, scalar: float):
        return self.__mul__(scalar)
        
    def __truediv__(self, scalar: float):
        return Vec2(self.x / scalar, self.y / scalar)

@dataclass
class Vec3:
    x: float
    y: float
    z: float

@dataclass
class Quaternion:
    w: float
    x: float
    y: float
    z: float


    
    @classmethod
    def from_euler(cls, ra: float, dec: float, roll: float):
        cy = math.cos(ra * 0.5)
        sy = math.sin(ra * 0.5)
        cp = math.cos(dec * 0.5)
        sp = math.sin(dec * 0.5)
        cr = math.cos(roll * 0.5)
        sr = math.sin(roll * 0.5)

        w = cr * cp * cy + sr * sp * sy
        x = sr * cp * cy - cr * sp * sy
        y = cr * sp * cy + sr * cp * sy
        z = cr * cp * sy - sr * sp * cy

        return cls(w, x, y, z)
    
    @classmethod
    def from_axis_angle(cls, axis: tuple[float, float, float], angle: float):
        half_angle = angle * 0.5
        s = math.sin(half_angle)
        return cls(
            math.cos(half_angle),
            axis[0] * s,
            axis[1] * s,
            axis[2] * s
        )


    def conjugate(self):
        return Quaternion(self.w, -self.x, -self.y, -self.z)

    def __mul__(self, other):
        return Quaternion(
            self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
            self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w
        )

    def rotate(self, v: Vec3) -> Vec3:
        tx = 2 * (self.y * v.z - self.z * v.y)
        ty = 2 * (self.z * v.x - self.x * v.z)
        tz = 2 * (self.x * v.y - self.y * v.x)
        
        cx = self.y * tz - self.z * ty
        cy = self.z * tx - self.x * tz
        cz = self.x * ty - self.y * tx
        
        return Vec3(
            v.x + self.w * tx + cx,
            v.y + self.w * ty + cy,
            v.z + self.w * tz + cz
        )
        
    def angle(self) -> float:
            w_clamped = max(-1.0, min(1.0, self.w))
            return 2.0 * math.acos(w_clamped)
    

    def is_unit(self, tolerance : float) -> bool:
        return abs(self.x*self.x+self.y*self.y+self.z*self.z+self.w*self.w - 1) < tolerance
    
def spherical_to_quaternion(ra: float, dec: float, roll: float):
    if not (0.0 <= roll <= 2*math.pi) or not (0.0 <= ra <= 2*math.pi) or not (-math.pi <= dec <= math.pi):
        print("Please give RA, DEC, Roll in their principle form.")


    a = Quaternion.from_axis_angle((0.0, 0.0, 1.0), ra)
    b = Quaternion.from_axis_angle((0.0, 1.0, 0.0), -dec)
    c = Quaternion.from_axis_angle((1.0, 0.0, 0.0), -roll)
    
    result = (a * b * c).conjugate()

    assert result.is_unit(0.00001)
    return result

class Attitude:
    def __init__(self, quaternion: Quaternion):
        self.quaternion = quaternion
        
    def get_quaternion(self) -> Quaternion:
        return self.quaternion
        
    def rotate(self, spatial_vector: Vec3) -> Vec3:
        return self.quaternion.rotate(spatial_vector)

class Camera:
    def __init__(self, res_x: int, res_y: int, fov_degrees: float = 20.0):
        self.res_x = res_x
        self.res_y = res_y
        
        self.cx = res_x / 2.0
        self.cy = res_y / 2.0
        
        fov_radians = math.radians(fov_degrees)
        self.focal_length = self.cx / math.tan(fov_radians / 2.0)
        
    def x_resolution(self) -> int:
        return self.res_x
        
    def y_resolution(self) -> int:
        return self.res_y
        
    def spatial_to_camera(self, spatial: Vec3) -> Vec2:
        x_pixel = (-spatial.y / spatial.x) * self.focal_length + self.cx
        y_pixel = (-spatial.z / spatial.x) * self.focal_length + self.cy
        
        return Vec2(x_pixel, y_pixel)
        
    def in_sensor(self, coords: Vec2) -> bool:
        return 0 <= coords.x < self.res_x and 0 <= coords.y < self.res_y
        


class Star:
    def __init__(self, x: float, y: float, radius_x: float, radius_y: float, magnitude: float):
        self.position = Vec2(x, y)
        self.radius_x = radius_x
        self.radius_y = radius_y
        self.magnitude = magnitude # technically not a real magnitude since x_i > x_j if x_i is brighter than x_j. but ive not implemented this because i think its only used in centroiding?

@dataclass
class GeneratedStar:
    star: Star
    peak_brightness_per_time: float
    delta: Vec2

@dataclass
class StarIdentifier:
    input_index: int
    catalog_index: int

def mag_to_brightness(magnitude: float) -> float:
    return math.pow(10, -magnitude/250.0)

def centroid_imaging_probability(magnitude: float, cutoff_mag: float) -> float:
    brightness = mag_to_brightness(magnitude)
    cutoff_brightness = mag_to_brightness(cutoff_mag)
    stddev = cutoff_brightness/5.0

    return 1 - (0.5 * (1 + math.erf((cutoff_brightness-brightness)/(stddev*math.sqrt(2.0)))))


def static_pixel_brightness(pixel: Vec2, g_star: GeneratedStar, exposure_time: float, spread_std_dev: float) -> float:
    dx = g_star.star.position.x - pixel.x
    dy = g_star.star.position.y - pixel.y

    return g_star.peak_brightness_per_time * exposure_time * math.exp(-math.fma(dx,dx, dy*dy) / (2 * spread_std_dev**2))

def motion_blurred_pixel_brightness(sample_pos: Vec2, g_star: GeneratedStar, time_t: float, spread_std_dev: float) -> float:
    # i dont have motion blur added yet. placeholder.
    return static_pixel_brightness(sample_pos, g_star, time_t, spread_std_dev)