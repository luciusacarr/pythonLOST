# pipeline.py
import math
import numpy as np
import copy


from catalog import CatalogStar 

from geometry import (
     Vec2, Star, GeneratedStar, StarIdentifier, 
     mag_to_brightness, centroid_imaging_probability, 
     motion_blurred_pixel_brightness, static_pixel_brightness
)

class StarfieldSimulator:
    def __init__(self, catalog: list[CatalogStar], attitude, camera, rng: np.random.Generator,
                 centroids_only: bool, zero_mag_total_photons: float, 
                 star_spread_std_dev: float, saturation_photons: float, 
                 dark_current: float, read_noise_std_dev: float, 
                 motion_blur_direction, exposure_time: float, 
                 readout_time: float, shot_noise: bool, 
                 oversampling: int, num_false_stars: int, 
                 false_star_min_magnitude: int, false_star_max_magnitude: int, 
                 cutoff_mag: int, perturbation_stddev: float):
        
        assert false_star_max_magnitude <= false_star_min_magnitude
        assert perturbation_stddev >= 0.0

        self.camera = camera
        self.attitude = attitude
        self.catalog = catalog
        
        self.width = camera.x_resolution()
        self.height = camera.y_resolution()
        
        assert oversampling >= 1
        oversampling_per_axis = math.ceil(math.sqrt(oversampling))
        if oversampling_per_axis * oversampling_per_axis != oversampling:
            print(f"WARNING: oversampling was not a perfect square. Rounding up to {oversampling_per_axis**2}.")
            
        assert exposure_time > 0
        
        motion_blur_enabled = abs(motion_blur_direction.get_quaternion().angle()) > 0.001
        motion_blur_direction_q = motion_blur_direction.get_quaternion()
        current_attitude = attitude.get_quaternion()
        future_attitude = motion_blur_direction_q * current_attitude 

        
        generated_stars = []
        self.expected_stars = []
        self.expected_star_ids = []
        self.input_stars = []
        self.input_star_ids = []
        self.image = None


        zero_mag_peak_photon_density = zero_mag_total_photons / (2 * math.pi * star_spread_std_dev**2)

        catalog_with_false = copy.copy(catalog) 

        for _ in range(num_false_stars):
            ra = rng.uniform(0.0, 1.0) * 2 * math.pi
            de = math.asin(rng.uniform(0.0, 1.0) * 2 - 1)
            magnitude = rng.integers(false_star_max_magnitude, false_star_min_magnitude + 1)
            catalog_with_false.append(CatalogStar(ra=ra, dec=de, star_id=-1, flag="", magnitude=magnitude))

        for i, catalog_star in enumerate(catalog_with_false):
            is_true_star = i < len(catalog)

            rotated = attitude.rotate(catalog_star.spatial)
            if rotated.x <= 0:
                continue

            cam_coords = camera.spatial_to_camera(rotated)

            if camera.in_sensor(cam_coords):
                future_spatial = future_attitude.rotate(catalog_star.spatial)
                delta = camera.spatial_to_camera(future_spatial) - cam_coords
                if not motion_blur_enabled:
                    delta = Vec2(0, 0)

                peak_brightness_per_time = zero_mag_peak_photon_density * mag_to_brightness(catalog_star.magnitude)
                interesting_threshold = 0.05 
                
                val = -math.log(interesting_threshold / peak_brightness_per_time / exposure_time) * 2 * math.pi * star_spread_std_dev**2
                radius = math.ceil(math.sqrt(max(0, val))) 
                
                
                star = Star(cam_coords.x, cam_coords.y, radius, radius, -catalog_star.magnitude)
                generated_stars.append(GeneratedStar(star, peak_brightness_per_time, delta))

                self.expected_stars.append(star)
                if is_true_star:
                    self.expected_star_ids.append(StarIdentifier(len(self.expected_stars) - 1, i))

                input_star = copy.deepcopy(star)
                if perturbation_stddev > 0.0:
                    px = rng.normal(0.0, perturbation_stddev)
                    py = rng.normal(0.0, perturbation_stddev)
                    input_star.position.x += max(min(px, 2 * perturbation_stddev), -2 * perturbation_stddev)
                    input_star.position.y += max(min(py, 2 * perturbation_stddev), -2 * perturbation_stddev)

                if camera.in_sensor(input_star.position):
                    if (cutoff_mag >= 10000 or not is_true_star or 
                        rng.random() < centroid_imaging_probability(catalog_star.magnitude, cutoff_mag)):
                        self.input_stars.append(input_star)
                        if is_true_star:
                            self.input_star_ids.append(StarIdentifier(len(self.input_stars) - 1, i))

        if centroids_only:
            return

        photons_buffer = np.zeros(self.width * self.height, dtype=np.float64)
        oversampling_brightness_factor = oversampling_per_axis**2

        for g_star in generated_stars:
            star = g_star.star
            
            earliest_pos = star.position - g_star.delta * (exposure_time / 2.0 + readout_time / 2.0)
            latest_pos = star.position + g_star.delta * (exposure_time / 2.0 + readout_time / 2.0)
            
            x_min = max(0, int(min(earliest_pos.x - star.radius_x, latest_pos.x - star.radius_x)))
            x_max = min(self.width - 1, int(max(earliest_pos.x + star.radius_x, latest_pos.x + star.radius_x)))
            y_min = max(0, int(min(earliest_pos.y - star.radius_y, latest_pos.y - star.radius_y)))
            y_max = min(self.height - 1, int(max(earliest_pos.y + star.radius_y, latest_pos.y + star.radius_y)))

            for x_pixel in range(x_min, x_max + 1):
                for y_pixel in range(y_min, y_max + 1):
                    readout_offset = readout_time * (y_pixel - self.height / 2.0) / self.height
                    t_start = -exposure_time / 2.0 + readout_offset
                    t_end = exposure_time / 2.0 + readout_offset

                    for x_sample in range(oversampling_per_axis):
                        for y_sample in range(oversampling_per_axis):
                            x = x_pixel + (x_sample + 0.5) / oversampling_per_axis
                            y = y_pixel + (y_sample + 0.5) / oversampling_per_axis

                            sample_pos = Vec2(x, y)
                            if motion_blur_enabled:
                                cur_photons = (
                                    motion_blurred_pixel_brightness(sample_pos, g_star, t_end, star_spread_std_dev) - 
                                    motion_blurred_pixel_brightness(sample_pos, g_star, t_start, star_spread_std_dev)
                                ) / oversampling_brightness_factor
                            else:
                                cur_photons = static_pixel_brightness(sample_pos, g_star, exposure_time, star_spread_std_dev) / oversampling_brightness_factor

                            assert 0.0 <= cur_photons
                            photons_buffer[x_pixel + y_pixel * self.width] += cur_photons

        k_max_brightness = 255

        base_brightness = np.full(self.width * self.height, dark_current, dtype=np.float64)
        read_noise = rng.normal(0.0, read_noise_std_dev, size=self.width * self.height)

        if shot_noise:
            if np.max(photons_buffer) > (2**31 - 1): 
                raise ValueError("One of the pixels had too many photons. Generated image would not be physically accurate.")
            quantized_photons = rng.poisson(photons_buffer)
        else:
            quantized_photons = np.round(photons_buffer)

        cur_brightness = base_brightness + read_noise + (quantized_photons / saturation_photons)

        clamped_brightness = np.clip(cur_brightness, 0.0, 1.0)
        self.image_data = np.floor(clamped_brightness * k_max_brightness).astype(np.uint8)
        
        self.image = self.image_data.reshape((self.height, self.width))