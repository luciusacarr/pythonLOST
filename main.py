# main.py
import argparse
import time
import numpy as np
from PIL import Image
import math

from catalog import load_catalog
from pipeline import StarfieldSimulator
from geometry import Camera, Attitude, Quaternion, Vec3, spherical_to_quaternion

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    parser = argparse.ArgumentParser(description="Starfield Image Generator Pipeline")

    parser.add_argument("--output-name", type=str, default="starfield.png")

    parser.add_argument("--generate-x-resolution", type=int, default=1024)
    parser.add_argument("--generate-y-resolution", type=int, default=1024)
    
    parser.add_argument("--generate-centroids-only", type=str2bool, default=False)
    
    parser.add_argument("--generate-zero-mag-photons", type=float, default=20000)
    parser.add_argument("--generate-saturation-photons", type=float, default=150.0)
    parser.add_argument("--generate-spread-stddev", type=float, default=1.0)
    parser.add_argument("--generate-shot-noise", type=str2bool, default=True)
    parser.add_argument("--generate-dark-current", type=float, default=.1)
    parser.add_argument("--generate-read-noise-stddev", type=float, default=0.05)
    
    parser.add_argument("--generate-ra", type=float, default=88.0)
    parser.add_argument("--generate-de", type=float, default=7.0)
    parser.add_argument("--generate-roll", type=float, default=0.0)
    parser.add_argument("--generate-random-attitudes", type=str2bool, default=False)

    parser.add_argument("--generate-blur-ra", type=float, default=0.0)
    parser.add_argument("--generate-blur-de", type=float, default=0.0)
    parser.add_argument("--generate-blur-roll", type=float, default=0.0)

    parser.add_argument("--generate-exposure", type=float, default=0.2)
    parser.add_argument("--generate-readout-time", type=float, default=0.0)
    parser.add_argument("--generate-oversampling", type=int, default=4)

    parser.add_argument("--generate-false-stars", type=int, default=0)
    parser.add_argument("--generate-false-min-mag", type=float, default=8.0)
    parser.add_argument("--generate-false-max-mag", type=float, default=1.0)
    parser.add_argument("--generate-perturb-centroids", type=float, default=0.0)
    parser.add_argument("--generate-cutoff-mag", type=float, default=6.0)

    parser.add_argument("--generate-seed", type=int, default=394859)
    parser.add_argument("--generate-time-based-seed", type=str2bool, default=False)

    parser.add_argument("--fov", type=float, default=20.0, help="Camera field of view in degrees")

    args = parser.parse_args()

    print("Loading star catalog...")
    catalog = load_catalog("stars.tsv") 
    
    seed = int(time.time()) if args.generate_time_based_seed else args.generate_seed
    rng = np.random.default_rng(seed=seed) 

    camera = Camera(res_x=args.generate_x_resolution, res_y=args.generate_y_resolution, fov_degrees=args.fov)
    
    ra_rad = math.radians(args.generate_ra)
    dec_rad = math.radians(args.generate_de)
    roll_rad = math.radians(args.generate_roll)

    print(ra_rad)

    target_quat = spherical_to_quaternion(ra_rad, dec_rad, roll_rad)
    attitude = Attitude(quaternion=target_quat)

    print(f"\n--- DEBUG INFO ---")
    print(f"Input Angles: RA={args.generate_ra}°, Dec={args.generate_de}°, Roll={args.generate_roll}°")
    print(f"Resulting Quat: w={target_quat.w:.4f}, x={target_quat.x:.4f}, y={target_quat.y:.4f}, z={target_quat.z:.4f}")

    identity_quat = Quaternion(w=1.0, x=0.0, y=0.0, z=0.0)
    motion_blur = Attitude(quaternion=identity_quat)

    print("Initializing Starfield Simulator...")
    pipeline = StarfieldSimulator(
        catalog=catalog,
        attitude=attitude,
        camera=camera,
        rng=rng,
        centroids_only=args.generate_centroids_only,
        zero_mag_total_photons=args.generate_zero_mag_photons,
        star_spread_std_dev=args.generate_spread_stddev,
        saturation_photons=args.generate_saturation_photons,
        dark_current=args.generate_dark_current,
        read_noise_std_dev=args.generate_read_noise_stddev,
        motion_blur_direction=motion_blur,
        exposure_time=args.generate_exposure,
        readout_time=args.generate_readout_time,
        shot_noise=args.generate_shot_noise,
        oversampling=args.generate_oversampling,
        num_false_stars=args.generate_false_stars,
        false_star_min_magnitude=args.generate_false_min_mag,
        false_star_max_magnitude=args.generate_false_max_mag,
        cutoff_mag=args.generate_cutoff_mag,
        perturbation_stddev=args.generate_perturb_centroids
    )

    if not args.generate_centroids_only:
        print("Simulation complete! Generating image...")
        final_image = Image.fromarray(pipeline.image)
        output_filename = (args.output_name if args.output_name[-4::1] == ".png" else args.output_name + ".png")
        final_image.save(output_filename)

        # post effects can go here. i personally think pre-effects are easier to implement, for any future reference.
        
        print(f"Image saved successfully to {output_filename}!")
    else:
        print(f"Simulation complete! Generated {len(pipeline.expected_stars)} centroids.")

if __name__ == "__main__":
    main()