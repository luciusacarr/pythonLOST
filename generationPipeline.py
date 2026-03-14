# direct implementation of the c++ generation pipeline in python.
# not necessarily as fast, but keeps things separate.
# obviously if we wish to utilize real images this serves no purpose,
# but real images are not necessarily the goal.
"""GeneratedPipelineInput::GeneratedPipelineInput(const Catalog &catalog,
                                               Attitude attitude,
                                               Camera camera,
                                               std::default_random_engine *rng,

                                               bool centroidsOnly,
                                               decimal zeroMagTotalPhotons,
                                               decimal starSpreadStdDev,
                                               decimal saturationPhotons,
                                               decimal darkCurrent,
                                               decimal readNoiseStdDev,
                                               Attitude motionBlurDirection, // applied on top of the attitude
                                               decimal exposureTime,
                                               decimal readoutTime, // zero for no rolling shutter
                                               bool shotNoise,
                                               int oversampling,
                                               int numFalseStars,
                                               int falseStarMinMagnitude,
                                               int falseStarMaxMagnitude,
                                               int cutoffMag,
                                               decimal perturbationStddev)
    : camera(camera), attitude(attitude), catalog(catalog) {

    assert(falseStarMaxMagnitude <= falseStarMinMagnitude);
    assert(perturbationStddev >= DECIMAL(0.0));



    image.width = camera.XResolution();
    image.height = camera.YResolution();
    // number of true photons each pixel receives.

    assert(oversampling >= 1);
    int oversamplingPerAxis = DECIMAL_CEIL(DECIMAL_SQRT(oversampling));
    if (oversamplingPerAxis*oversamplingPerAxis != oversampling) {
        std::cerr << "WARNING: oversampling was not a perfect square. Rounding up to "
                  << oversamplingPerAxis*oversamplingPerAxis << "." << std::endl;
    }
    assert(exposureTime > 0);
    bool motionBlurEnabled = abs(motionBlurDirection.GetQuaternion().Angle()) > DECIMAL(0.001);
    Quaternion motionBlurDirectionQ = motionBlurDirection.GetQuaternion();
    // attitude at the middle of exposure time
    Quaternion currentAttitude = attitude.GetQuaternion();
    // attitude 1 time unit after middle of exposure
    Quaternion futureAttitude = motionBlurDirectionQ*currentAttitude;
    std::vector<GeneratedStar> generatedStars;

    // a star with 1 photon has peak density 1/(2pi sigma^2), because 2d gaussian formula. Then just
    // multiply up proportionally!
    decimal zeroMagPeakPhotonDensity = zeroMagTotalPhotons / (2*DECIMAL_M_PI * starSpreadStdDev*starSpreadStdDev);

    // TODO: Is it 100% correct to just copy the standard deviation in both dimensions?
    std::normal_distribution<decimal> perturbation1DDistribution(DECIMAL(0.0), perturbationStddev);

    Catalog catalogWithFalse = catalog;

    std::uniform_real_distribution<decimal> uniformDistribution(DECIMAL(0.0), DECIMAL(1.0));
    std::uniform_int_distribution<int> magnitudeDistribution(falseStarMaxMagnitude, falseStarMinMagnitude);
    for (int i = 0; i < numFalseStars; i++) {
        decimal ra = uniformDistribution(*rng) * 2*DECIMAL_M_PI;
        // to be uniform around sphere. Borel-Kolmogorov paradox is calling
        decimal de = DECIMAL_ASIN(uniformDistribution(*rng)*2 - 1);
        decimal magnitude = magnitudeDistribution(*rng);

        catalogWithFalse.push_back(CatalogStar(ra, de, magnitude, -1));
    }

    for (int i = 0; i < (int)catalogWithFalse.size(); i++) {
        bool isTrueStar = i < (int)catalog.size();

        const CatalogStar &catalogStar = catalogWithFalse[i];
        Vec3 rotated = attitude.Rotate(catalogWithFalse[i].spatial);
        if (rotated.x <= 0) {
            continue;
        }
        Vec2 camCoords = camera.SpatialToCamera(rotated);

        if (camera.InSensor(camCoords)) {
            Vec3 futureSpatial = futureAttitude.Rotate(catalogWithFalse[i].spatial);
            Vec2 delta = camera.SpatialToCamera(futureSpatial) - camCoords;
            if (!motionBlurEnabled) {
                delta = {0, 0}; // avoid decimaling point funny business
            }
            // radiant intensity, in photons per time unit per pixel, at the center of the star.
            decimal peakBrightnessPerTime = zeroMagPeakPhotonDensity * MagToBrightness(catalogStar.magnitude);
            decimal interestingThreshold = DECIMAL(0.05); // we don't need to check pixels that are expected to
                                               // receive this many photons or fewer.
            // inverse of the function defining the Gaussian distribution: Find out how far from the
            // mean we'll have to go until the number of photons is less than interestingThreshold
            decimal radius = DECIMAL_CEIL(DECIMAL_SQRT(-DECIMAL_LOG(interestingThreshold/peakBrightnessPerTime/exposureTime)*2*DECIMAL_M_PI*starSpreadStdDev*starSpreadStdDev));
            Star star = Star(camCoords.x, camCoords.y,
                             radius, radius,
                             // important to invert magnitude here, so that centroid magnitude becomes larger for brighter stars.
                             // It's possible to make it so that the magnitude is always positive too, but allowing weirder magnitudes helps keep star-id algos honest about their assumptions on magnitude.
                             // we don't use its magnitude anywhere else in generation; peakBrightness was already calculated.
                             -catalogStar.magnitude);
            generatedStars.push_back(GeneratedStar(star, peakBrightnessPerTime, delta));

            // Now add the star to the input and expected lists.
            // We do actually want to add false stars as well, because:
            // a) A centroider isn't any worse because it picks up a false star that looks exactly like a normal star, so why should we exclude them from compare-centroids?
            // b) We want to feed false stars into star-ids to evaluate their false-star resilience without running centroiding.

            // Add all stars to expected, cause that's how we roll
            expectedStars.push_back(star);
            // and provide expected identifications for all of them
            if (isTrueStar) {
                expectedStarIds.push_back(StarIdentifier(expectedStars.size()-1, i));
            }

            // for input, though, add perturbation and stuff.
            Star inputStar = star;
            if (perturbationStddev > DECIMAL(0.0)) {
                // clamp to within 2 standard deviations for some reason:
                inputStar.position.x += std::max(std::min(perturbation1DDistribution(*rng), 2*perturbationStddev), -2*perturbationStddev);
                inputStar.position.y += std::max(std::min(perturbation1DDistribution(*rng), 2*perturbationStddev), -2*perturbationStddev);
            }
            // If it got perturbed outside of the sensor, don't add it.
            if (camera.InSensor(inputStar.position)
                // and also don't add it if it's too dim.
                && (cutoffMag >= 10000 // but always add the star if the cutoff is very high
                    || !isTrueStar // and always add the false stars
                    || std::bernoulli_distribution(CentroidImagingProbability(catalogStar.magnitude, cutoffMag))(*rng))) {
                inputStars.push_back(inputStar);
                if (isTrueStar) {
                    inputStarIds.push_back(StarIdentifier(inputStars.size()-1, i));
                }
            }
        }
    }

    if (centroidsOnly) {
        // we're outta here
        return;
    }

    std::vector<decimal> photonsBuffer(image.width*image.height, 0);

    for (const GeneratedStar &star : generatedStars) {
        // delta will be exactly (0,0) when motion blur disabled
        Vec2 earliestPosition = star.position - star.delta*(exposureTime/DECIMAL(2.0) + readoutTime/DECIMAL(2.0));
        Vec2 latestPosition = star.position + star.delta*(exposureTime/DECIMAL(2.0) + readoutTime/DECIMAL(2.0));
        int xMin = std::max(0, (int)std::min(earliestPosition.x - star.radiusX, latestPosition.x - star.radiusX));
        int xMax = std::min(image.width-1, (int)std::max(earliestPosition.x + star.radiusX, latestPosition.x + star.radiusX));
        int yMin = std::max(0, (int)std::min(earliestPosition.y - star.radiusX, latestPosition.y - star.radiusX));
        int yMax = std::min(image.height-1, (int)std::max(earliestPosition.y + star.radiusX, latestPosition.y + star.radiusX));

        // peak brightness is measured in photons per time unit per pixel, so if oversampling, we
        // need to convert units to photons per time unit per sample
        decimal oversamplingBrightnessFactor = oversamplingPerAxis*oversamplingPerAxis;

        // the star.x and star.y refer to the pixel whose top left corner the star should appear at
        // (and fractional amounts are relative to the corner). When we color a pixel, we ideally
        // would integrate the intensity of the star over that pixel, but we can make do by sampling
        // the intensity of the star at the /center/ of the pixel, ie, star.x+.5 and star.y+.5
        for (int xPixel = xMin; xPixel <= xMax; xPixel++) {
            for (int yPixel = yMin; yPixel <= yMax; yPixel++) {
                // offset of beginning & end of readout compared to beginning & end of readout for
                // center row
                decimal readoutOffset = readoutTime * (yPixel - image.height/DECIMAL(2.0)) / image.height;
                decimal tStart = -exposureTime/DECIMAL(2.0) + readoutOffset;
                decimal tEnd = exposureTime/DECIMAL(2.0) + readoutOffset;

                // loop through all samples in the current pixel
                for (int xSample = 0; xSample < oversamplingPerAxis; xSample++) {
                    for (int ySample = 0; ySample < oversamplingPerAxis; ySample++) {
                        decimal x = xPixel + (xSample+DECIMAL(0.5))/oversamplingPerAxis;
                        decimal y = yPixel + (ySample+DECIMAL(0.5))/oversamplingPerAxis;

                        decimal curPhotons;
                        if (motionBlurEnabled) {
                            curPhotons =
                                (MotionBlurredPixelBrightness({x, y}, star, tEnd, starSpreadStdDev)
                                 - MotionBlurredPixelBrightness({x, y}, star, tStart, starSpreadStdDev))
                                / oversamplingBrightnessFactor;
                        } else {
                            curPhotons = StaticPixelBrightness({x, y}, star, exposureTime, starSpreadStdDev)
                                / oversamplingBrightnessFactor;
                        }

                        assert(DECIMAL(0.0) <= curPhotons);

                        photonsBuffer[xPixel + yPixel*image.width] += curPhotons;
                    }
                }
            }
        }
    }

    std::normal_distribution<decimal> readNoiseDist(DECIMAL(0.0), readNoiseStdDev);

    // convert from photon counts to observed pixel brightnesses, applying noise and such.
    imageData = std::vector<unsigned char>(image.width*image.height);
    image.image = imageData.data();
    for (int i = 0; i < image.width * image.height; i++) {
        decimal curBrightness = 0;

        // dark current (Constant)
        curBrightness += darkCurrent;

        // read noise (Gaussian)
        curBrightness += readNoiseDist(*rng);

        // shot noise (Poisson), and quantize
        long quantizedPhotons;
        if (shotNoise) {
            // with GNU libstdc++, it keeps sampling from the distribution until it's within the min-max
            // range. This is problematic if the mean is far above the max long value, because then it
            // might have to sample many many times (and furthermore, the results won't be useful
            // anyway)
            decimal photons = photonsBuffer[i];
            if (photons > DECIMAL(LONG_MAX) - DECIMAL(3.0) * DECIMAL_SQRT(LONG_MAX)) {
                std::cout << "ERROR: One of the pixels had too many photons. Generated image would not be physically accurate, exiting." << std::endl;
                exit(1);
            }
            std::poisson_distribution<long> shotNoiseDist(photonsBuffer[i]);
            quantizedPhotons = shotNoiseDist(*rng);
        } else {
            quantizedPhotons = round(photonsBuffer[i]);
        }
        curBrightness += quantizedPhotons / saturationPhotons;

        // std::clamp not introduced until C++17, so we avoid it.
        decimal clampedBrightness = std::max(std::min(curBrightness, DECIMAL(1.0)), DECIMAL(0.0));
        imageData[i] = floor(clampedBrightness * kMaxBrightness); // TODO: off-by-one, 256?
    }
}"""


import argparse
from PIL import Image

parser = argparse.ArgumentParser(description='Generate an image based on width and height.')
parser.add_argument('--width', type=int, default=400, help='Width of the image')
parser.add_argument('--height', type=int, default=300, help='Height of the image')
parser.add_argument('--fov', type=int, default=300, help='Field of view')
parser.add_argument('--ra', type=int, default=0, help='Right ascension')
parser.add_argument('--dec', type=int, default=0, help='Declination')
parser.add_argument('--roll', type=int, default=0, help='Roll angle')
parser.add_argument('--exposure', type=int, default=1, help='Exposure time')
parser.add_argument('--output', type=str, default='generated_image.png', help='Output filename')


class GenerationPipeline:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.image = Image.new('RGB', (width, height), color='white')

    def generate(self):

        for x in range(self.width):
            for y in range(self.height):
                r = int((x / self.width) * 255)
                g = int((y / self.height) * 255)
                b = 128  
                self.image.putpixel((x, y), (r, g, b))

    def save(self, filename):
        self.image.save(filename)

if __name__ == "__main__":
    # user can prompt tags like python3 generationPipeline.py --width 400 --height 300 --output generated_image.png
    args = parser.parse_args()

    pipeline = GenerationPipeline(args.width, args.height)
    pipeline.generate()
    pipeline.save(args.output)
