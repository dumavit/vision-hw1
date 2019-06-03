#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float sum = 0;

    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
            sum += get_pixel(im, i, j, 0);

    // Use only 2D filters
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float pixel = get_pixel(im, i, j, 0);
            set_pixel(im, i, j, 0, pixel / sum);
        }
}

image make_box_filter(int w)
{
    // TODO
    image im = make_image(w, w, 1);
    for (int i = 0; i < w; i++)
        for (int j = 0; j < w; j++)
            set_pixel(im, i, j, 0, 1.0 / (w * w));

    return im;
}

float get_conv(image filter, int filter_channel, image im, int x, int y, int im_channel)
{
    float result = 0;
    float half_filter_size = filter.w / 2.0;

    for (int i = 0; i < filter.w; i++)
        for (int j = 0; j < filter.h; j++)
        {
            float filter_value = get_pixel(filter, i, j, filter_channel);
            float image_color = get_pixel(im, x - half_filter_size / 2 + i, y - half_filter_size + j, im_channel);
            result += filter_value * image_color;
        }

    return result;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    assert(filter.c == im.c || filter.c == 1);
    if (preserve == 1)
        assert(filter.c == 1);

    int output_channels = preserve == 1 ? im.c : 1;

    image out = make_image(im.w, im.h, output_channels);

    float conv_result;

    // Iterate through image
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            if (preserve == 1)
            {
                // Iterate through channels
                for (int c = 0; c < im.c; c++)
                {
                    conv_result = get_conv(filter, 0, im, i, j, c);
                    set_pixel(out, i, j, c, conv_result);
                }
            }
            else
            {
                conv_result = 0;
                for (int c = 0; c < im.c; c++)
                    conv_result += get_conv(filter, c, im, i, j, c);

                set_pixel(out, i, j, 0, conv_result);
            }
        }

    return out;
}

image make_highpass_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);

    filter.data[0] = 0;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 4;
    filter.data[5] = -1;
    filter.data[6] = 0;
    filter.data[7] = -1;
    filter.data[8] = 0;

    return filter;
}

image make_sharpen_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);

    filter.data[0] = 0;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 5;
    filter.data[5] = -1;
    filter.data[6] = 0;
    filter.data[7] = -1;
    filter.data[8] = 0;

    return filter;
}

image make_emboss_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);

    filter.data[0] = -2;
    filter.data[1] = -1;
    filter.data[2] = 0;
    filter.data[3] = -1;
    filter.data[4] = 1;
    filter.data[5] = 1;
    filter.data[6] = 0;
    filter.data[7] = 1;
    filter.data[8] = 2;

    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Highpass filter shouldn't preserve channels

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Have to clamp images after convolutions

float gauss(float sigma, int x, int y)
{
    float PI = 3.14159265359;
    float sigma_squared = powf(sigma, 2);
    float exp_argument = -(x * x + y * y) / (2 * sigma_squared);

    return expf(exp_argument) / (2 * PI * sigma_squared);
}

image make_gaussian_filter(float sigma)
{
    // TODO
    image filter = make_image(3, 3, 1);

    filter.data[0] = gauss(sigma, -1, -1);
    filter.data[1] = gauss(sigma, -1, 0);
    filter.data[2] = gauss(sigma, -1, 1);

    filter.data[3] = gauss(sigma, 0, -1);
    filter.data[4] = gauss(sigma, 0, 0);
    filter.data[5] = gauss(sigma, 0, 1);

    filter.data[6] = gauss(sigma, 1, -1);
    filter.data[7] = gauss(sigma, 1, 0);
    filter.data[8] = gauss(sigma, 1, 1);

    l1_normalize(filter);
    return filter;
}

image add_image(image a, image b)
{
    assert(a.c == b.c && a.h == b.h && a.w == b.w);
    image out = make_image(a.w, a.h, a.c);

    for (int i = 0; i < a.w; i++)
        for (int j = 0; j < a.h; j++)
            for (int k = 0; k < a.c; k++)
            {
                float pixel_a = get_pixel(a, i, j, k);
                float pixel_b = get_pixel(b, i, j, k);
                set_pixel(out, i, j, k, pixel_a + pixel_b);
            }

    return out;
}

image sub_image(image a, image b)
{
    assert(a.c == b.c && a.h == b.h && a.w == b.w);
    image out = make_image(a.w, a.h, a.c);

    for (int i = 0; i < a.w; i++)
        for (int j = 0; j < a.h; j++)
            for (int k = 0; k < a.c; k++)
            {
                float pixel_a = get_pixel(a, i, j, k);
                float pixel_b = get_pixel(b, i, j, k);
                set_pixel(out, i, j, k, pixel_a - pixel_b);
            }

    return out;
}

image make_gx_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);

    filter.data[0] = -1;
    filter.data[1] = 0;
    filter.data[2] = 1;
    filter.data[3] = -2;
    filter.data[4] = 0;
    filter.data[5] = 2;
    filter.data[6] = -1;
    filter.data[7] = 0;
    filter.data[8] = 1;

    return filter;
}

image make_gy_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);

    filter.data[0] = -1;
    filter.data[1] = -2;
    filter.data[2] = -1;
    filter.data[3] = 0;
    filter.data[4] = 0;
    filter.data[5] = 0;
    filter.data[6] = 1;
    filter.data[7] = 2;
    filter.data[8] = 1;

    return filter;
}

void feature_normalize(image im)
{
    for (int channel = 0; channel < im.c; channel++)
    {
        float min_pixel = INFINITY;
        float max_pixel = -INFINITY;

        // Find minimum and maximum pixel value in this channel
        for (int i = 0; i < im.w; i++)
            for (int j = 0; j < im.h; j++)
            {
                float pixel = get_pixel(im, i, j, channel);
                if (pixel > max_pixel)
                    max_pixel = pixel;
                if (pixel < min_pixel)
                    min_pixel = pixel;
            }

        // Set normalized value
        for (int i = 0; i < im.w; i++)
            for (int j = 0; j < im.h; j++)
            {
                float pixel = get_pixel(im, i, j, channel);
                float normalized_pixel = max_pixel - min_pixel > 0
                                             ? (pixel - min_pixel) / (max_pixel - min_pixel)
                                             : 0;

                set_pixel(im, i, j, channel, normalized_pixel);
            }
    }
}

image *sobel_image(image im)
{
    // TODO
    image *out = calloc(2, sizeof(image));

    image magnitude = make_image(im.w, im.h, 1);
    image direction = make_image(im.w, im.h, 1);

    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    image gx_image = convolve_image(im, gx_filter, 0);
    image gy_image = convolve_image(im, gy_filter, 0);

    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float gx_pixel = get_pixel(gx_image, i, j, 0);
            float gy_pixel = get_pixel(gy_image, i, j, 0);

            float m = sqrtf(gy_pixel * gy_pixel + gx_pixel * gx_pixel);
            set_pixel(magnitude, i, j, 0, m);

            float d = atanf(gy_pixel / gx_pixel);
            // Check nan value

            if (isnan(d))
                set_pixel(direction, i, j, 0, 0);
            else
                set_pixel(direction, i, j, 0, d);
        }

    out[0] = magnitude;
    out[1] = direction;
    return out;
}

image colorize_sobel(image im)
{
    image *sobel = sobel_image(im);
    image magnitude = sobel[0];
    image direction = sobel[1];
    image out = make_image(im.w, im.h, im.c);

    feature_normalize(magnitude);
    feature_normalize(direction);

    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float m = get_pixel(magnitude, i, j, 0);
            float d = get_pixel(direction, i, j, 4);

            set_pixel(out, i, j, 0, d);
            set_pixel(out, i, j, 1, m);
            set_pixel(out, i, j, 2, m);
        }

    hsv_to_rgb(out);
    return out;
}
