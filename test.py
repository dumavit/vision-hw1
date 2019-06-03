from uwimg import *
im = load_image("data/dogsmall.jpg")
a = nn_resize(im, im.w*4, im.h*4)
save_image(a, "out/dog4x-nn")


a = bilinear_resize(im, im.w*4, im.h*4)
save_image(a, "out/dog4x-bl")


im = load_image("data/dog.jpg")
a = nn_resize(im, im.w//7, im.h//7)
save_image(a, "out/dog7th-bl")


im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
save_image(blur, "out/dog-box7")


im = load_image("data/dog.jpg")
f = make_box_filter(7)
blur = convolve_image(im, f, 1)
thumb = nn_resize(blur, blur.w//7, blur.h//7)
save_image(thumb, "out/dogthumb")


im = load_image("data/melisa.png")
f = make_highpass_filter()
highpass = convolve_image(im, f, 1)
clamp_image(highpass)
save_image(highpass, "out/melisa-highpass")


im = load_image("data/melisa.png")
f = make_sharpen_filter()
sharpen = convolve_image(im, f, 1)
clamp_image(sharpen)
save_image(sharpen, "out/melisa-sharpen")


im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
blur = convolve_image(im, f, 1)
save_image(blur, "out/dog-gauss")


im = load_image("data/dog.jpg")
f = make_gaussian_filter(2)
lfreq = convolve_image(im, f, 1)
hfreq = im - lfreq
reconstruct = lfreq + hfreq
save_image(lfreq, "out/low-frequency")
save_image(hfreq, "out/high-frequency")
save_image(reconstruct, "out/reconstruct")

im = load_image("data/dog.jpg")
res = sobel_image(im)
mag = res[0]
feature_normalize(mag)
save_image(mag, "out/magnitude")

im = load_image("data/dog.jpg")
coloured_sobel = colorize_sobel(im)
f = make_gaussian_filter(1)
blured = convolve_image(coloured_sobel, f, 1)
save_image(blured, "out/coloured_sobel")

