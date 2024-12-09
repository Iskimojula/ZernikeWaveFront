import aberrationRendering
from aberrationRendering import Zernike
import numpy
import pylab


# Wavefront representation
array_size = 512
n_radial_orders = 5

z = Zernike(print_indices = True)

r = z.r
theta = z.theta
indices = z.index_mapping
zernike_matrix = z.zernike_matrix

'''
We can represent any wavefront error by computing the matrix multiplication of the Zernike matrix with the Zernike
 coefficient vector.
'''
zernike_coefficients = numpy.zeros(zernike_matrix.shape[0])
zernike_coefficients[4] = 0.2e-6
#zernike_coefficients[8] = 0.3e-6
wavefront_error = aberrationRendering.wavefront(zernike_coefficients,zernike_matrix)

pylab.figure(1)
pylab.imshow(wavefront_error, interpolation='nearest', cmap='gray')

'''
we specify the pupil amplitude function
'''
amplitude_function = aberrationRendering.amplitudeFunctionCircle(array_size)

pylab.figure(2)
pylab.imshow(amplitude_function, interpolation='nearest', cmap='gray')


'''
 check that there are a sufficient number of pixels across the wavefront
'''
oversampling = 4
wavefront_gradient_check,psf_fov_check = \
aberrationRendering.checkSampling(wavefront_error,amplitude_function,array_size,oversampling)


#PSF generation
'''
The PSF is the Fourier transform of the complex pupil function
'''
wavelength = 550e-9
pupil_radius = 2.5e-3
psf = aberrationRendering.PSF(wavefront_error, amplitude_function, wavelength=wavelength, oversampling= oversampling)

pylab.figure(3)
pylab.imshow(psf, interpolation='nearest', cmap='gray')

psf_diffraction_limited = aberrationRendering.PSF(zernike_matrix[0]*0.0, amplitude_function, wavelength=wavelength, oversampling= oversampling)

#Convolution of the PSF with an input intensity pattern
pixel_scale = aberrationRendering.psfPixelScale(wavelength=wavelength,oversampling=oversampling,pupil_radius=pupil_radius)
pixel_scale_arcminutes = pixel_scale * 180 * 60/ numpy.pi

field_of_view_arcminutes = 120.0
n_pixels = aberrationRendering.numberPixels(field_of_view_arcminutes,
wavelength=wavelength,oversampling = oversampling,pupil_radius=pupil_radius)

square = numpy.ones((n_pixels,n_pixels),numpy.uint8) * 255
square[int(n_pixels/4):int(n_pixels*3/4),int(n_pixels/4):int(n_pixels*3/4)] = 0


pylab.figure(5)
pylab.imshow(square, interpolation='nearest', cmap='gray')


square_full = numpy.ones((psf.shape),numpy.uint8)*255
square_full[:square.shape[0],:square.shape[1]] = square

output_intensity_pattern = aberrationRendering.convolveImage(square_full,psf)

pylab.figure(6)
pylab.imshow(output_intensity_pattern, interpolation='nearest', cmap='gray')
pylab.show()