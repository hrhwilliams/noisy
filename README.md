# noisy
This project contains implementations of gradient noise algorithms such as Perlin noise and Simplex noise.
It also contains a simple PNG encoder to write a heightmap as the output of the algorithm.

Only Simplex noise is implemented right now.

Running ./noise will prompt for: width, height, scale, octaves, persistence, and lacunarity. It will generate
a PNG file with the seed of the function as its filename. If a number is passed to the program, e.g. './noise 5',
it will use that number as the seed.

Scale around 80, 3 octaves, persistence of 0.333 and lacunarity of 2.5 generates an okay height map

Terms -

Width: width of the image in pixels

Height: height of the image in pixels

Scale: Scale of the noise function
Octaves: Number of times to reapply the noise function to itself.

Persistence: Factor to scale how much successive calls contribute to the overall function. Values less than
             one will have each octave contribute less and less, and values more than one will cause each
             octave to contribute more than the previous one.
             
Lacunarity: Factor to scale the frequency of each successive call by. Higher values mean each octave higher
            will have a smaller scale, and lower values vice-versa.

### Examples
Image width: 500
Image height: 500
Scale: 75
Octaves: 2
Persistence: 0.5
Lacunarity: 2

![](examples/noise.png)

(wow!)
