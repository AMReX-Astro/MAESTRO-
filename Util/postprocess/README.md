# Rotating fully convective star diagnostics

Process a fully convective star problem to produce a set of diagnostic measures
including:

- **Convection speed**: The root-mean-squared radial velocity

<img src="https://render.githubusercontent.com/render/math?math=\large v_c(r) = \sqrt{\frac{1}{4\pi}\int_{0}^{2\pi}d\phi \int_0^{\pi}d\theta \sin\theta v_r^2}">

- **Ratio of mean rotation rates to surface rotation**: 

<img src="https://render.githubusercontent.com/render/math?math=\large \text{ratio}(r) = \dfrac{\int_V dV \Omega(r,\theta,\phi)}{\int_V dV \Omega(\text{surface},\theta,\phi)}">

Note that the denominator is currently set at the surface, but it does not strictly need to be at the surface; it just needs to be at a fixed radius.

- **Meridional circulation**: There are two components, radial and polar circulation

<img src="https://render.githubusercontent.com/render/math?math=\large u_r(r,\theta,t) = \frac{1}{2\pi t} \int_0^t dt \int_0^{2\pi}d\phi v_r(r,\theta,\phi,t)">, 

<img src="https://render.githubusercontent.com/render/math?math=\large u_{\theta}(r,\theta,t) = \frac{1}{2\pi t} \int_0^t dt \int_0^{2\pi}d\phi v_{\theta}(r,\theta,\phi,t)">

where <img src="https://render.githubusercontent.com/render/math?math=v_r"> and <img src="https://render.githubusercontent.com/render/math?math=v_{\theta}"> are the radial and polar velocities.
The instantaneous circulation components are computed by neglecting the time-averaging terms.

- **Baroclinity**: Averaged over the azimuthal angle

<img src="https://render.githubusercontent.com/render/math?math=\large \xi(r,\theta) = \frac{1}{2\pi} \int_0^{2\pi}d\phi\:\: \dfrac{\hat{\phi}\cdot (\nabla\ln P \times \nabla s)}{|\nabla \ln P|\: |\nabla s|}">.


- **Latitudinal shear**: Expressed as the projection of the rotation profile onto the (*l=2*,*m=0*) spherical harmonic, averaged over the surface. This can be computed in Cartesian coordinates

<img src="https://render.githubusercontent.com/render/math?math=\large \Omega_2(r) = \int_V dV \Omega(x,y,z)\: Y_{2,0}(x,y,z)\: K(r%27-r,\delta r)">

using a normalized Gaussian kernel 

<img src="https://render.githubusercontent.com/render/math?math=K(r%27-r,\delta r) = \dfrac{1}{\sqrt{\pi\delta r}} \exp[-(r%27-r)^2/\delta r \quad,\quad r%27=\sqrt{x^2%2By^2%2Bz^2}">

with kernel width

<img src="https://render.githubusercontent.com/render/math?math=\delta r = f(r) \Delta x^2 \quad,\quad f(r) = \min\{r/L,0.5\}">

where *L* is the half length of the problem domain.
Here,

<img src="https://render.githubusercontent.com/render/math?math=Y_{2,0}(x,y,z) = \dfrac{1}{4}\sqrt{\dfrac{5}{\pi}} \:\dfrac{2z^2-x^2-y^2}{r%27^2}">

- **Brünt-Väisälä frequency**:

<img src="https://render.githubusercontent.com/render/math?math=\large N^2 = -\dfrac{\gamma-1}{\gamma}\:\mathbf{g}\cdot\nabla s">

where <img src="https://render.githubusercontent.com/render/math?math=\mathbf{g}"> is the acceleration due to gravity and *s* is a dimensionless entropy.


## Building & running

The n-dimensional diagnostic can be built by executing `make DIM=n`. This will
produce the executable `radial_nd.exe`. Note that `DIM=3` is default setting and
the only one that is currently supported. To run, the executable must be provided
with the name of the plotfile to be analyzed:
```
./radial_3d.exe infile=plotfile_name
```

Additional arguments are as follows:

- **model file diagnostics**: extra argument `modelfile=modelfile_name` needs to be provided
- **exact solution test problem**: if `plotfile_name` is not provided, the code will write and evaluate a test problem with an exact solution.
- **time-averaged values**: two extra arguments need to be provided
  1. `dt=timesteps_between_plotfiles`, which is an integer value
  2. `nfiles=num_plotfiles` indicates the number of *additional* plotfiles needed, excluding the `infile` plotfile.
