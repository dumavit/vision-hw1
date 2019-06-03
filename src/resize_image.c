#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    x = round(x);
    y = round(y);

    return get_pixel(im, x, y, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image out = make_image(w, h, im.c);
    float x_mult = (float)im.w / w;
    float y_mult = (float)im.h / h;

    for (int i = 0; i < w; i++)
        for (int j = 0; j < h; j++)
            for (int k = 0; k < im.c; k++)
            {
                float interpolated = nn_interpolate(im, i * x_mult, j * y_mult, k);
                set_pixel(out, i, j, k, interpolated);
            }
    return out;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    float x_floor = floorf(x);
    float y_floor = floorf(y);

    float a1 = (x_floor+1-x) * (y_floor+1-y);
    float a2 = (x - x_floor) * (y_floor+1-y);
    float a3 = (x_floor+1-x) * (y - y_floor);
    float a4 = (x - x_floor) * (y - y_floor);

    float v1 = get_pixel(im, x_floor, y_floor, c);
    float v2 = get_pixel(im, x_floor+1, y_floor, c);
    float v3 = get_pixel(im, x_floor, y_floor+1, c);
    float v4 = get_pixel(im, x_floor+1, y_floor+1, c);

    return a1 * v1 + a2 * v2 + a3 * v3 + a4 * v4;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image out = make_image(w, h, im.c);
    float x_mult = (float)im.w / w;
    float y_mult = (float)im.h / h;

    for (int i = 0; i < w; i++)
        for (int j = 0; j < h; j++)
            for (int k = 0; k < im.c; k++)
            {
                float interpolated = bilinear_interpolate(im, i * x_mult, j * y_mult, k);
                set_pixel(out, i, j, k, interpolated);
            }

    return out;
}
