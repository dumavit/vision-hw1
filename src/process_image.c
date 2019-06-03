#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    // Clamp padding strategy
    x = x < 0 ? 0 : (x < im.w ? x : im.w - 1);
    y = y < 0 ? 0 : (y < im.h ? y : im.h - 1);
    c = c < 0 ? 0 : (c < im.c ? c : im.c - 1);

    int index = x + y * im.w + c * im.w * im.h;
    if (index < 0 || index >= im.w * im.h * im.c)
    {
        return 0;
    }
    return im.data[index];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    int index = x + y * im.w + c * im.w * im.h;

    if (index >= 0 && index < im.w * im.h * im.c)
    {
        im.data[index] = v;
    }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
            for (int k = 0; k < im.c; k++)
                set_pixel(im, i, j, k, get_pixel(im, i, j, k));

    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in

    for (int i = 0; i < im.w; i++)
    {
        for (int j = 0; j < im.h; j++)
        {
            float r = get_pixel(im, i, j, 0);
            float g = get_pixel(im, i, j, 1);
            float b = get_pixel(im, i, j, 2);

            float y = 0.299 * r + 0.587 * g + 0.114 * b;
            set_pixel(gray, i, j, 0, y);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float color = get_pixel(im, i, j, c);
            set_pixel(im, i, j, c, color + v);
        }
}

void clamp_image(image im)
{
    // TODO Fill this in
    // TODO Fill this in
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
            for (int k = 0; k < im.c; k++)
            {
                float color = get_pixel(im, i, j, k);
                color = color < 0 ? 0 : (color > 1 ? 1 : color);
                set_pixel(im, i, j, k, color);
            }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float r = get_pixel(im, i, j, 0);
            float g = get_pixel(im, i, j, 1);
            float b = get_pixel(im, i, j, 2);

            float V = three_way_max(r, g, b);
            float m = three_way_min(r, g, b);
            float C = V - m;
            float S = V == 0 ? 0 : C / V;

            float H;
            if (C == 0)
            {
                H = 0;
            }
            else
            {
                float H_prime;
                if (V == r)
                {
                    H_prime = (g - b) / C;
                }
                else if (V == g)
                {
                    H_prime = (b - r) / C + 2;
                }
                else
                {
                    H_prime = (r - g) / C + 4;
                }
                H = H_prime < 0 ? (H_prime / 6 + 1) : H_prime / 6;
            }

            set_pixel(im, i, j, 0, H);
            set_pixel(im, i, j, 1, S);
            set_pixel(im, i, j, 2, V);
        }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float H = get_pixel(im, i, j, 0);
            float S = get_pixel(im, i, j, 1);
            float V = get_pixel(im, i, j, 2);

            float H_prime = 6 * H;

            float C = V * S;
            float X = C * (1 - fabs(fmod(H_prime, 2) - 1));
            float m = V - C;
            int H_i = (int)6 * H;

            float r, g, b;
            if (H_i == 0)
            {
                r = C, g = X, b = 0;
            }
            else if (H_i == 1)
            {
                r = X, g = C, b = 0;
            }
            else if (H_i == 2)
            {
                r = 0, g = C, b = X;
            }
            else if (H_i == 3)
            {
                r = 0, g = X, b = C;
            }
            else if (H_i == 4)
            {
                r = X, g = 0, b = C;
            }
            else
            {
                r = C, g = 0, b = X;
            }

            set_pixel(im, i, j, 0, r + m);
            set_pixel(im, i, j, 1, g + m);
            set_pixel(im, i, j, 2, b + m);
        }
}

void scale_image(image im, int c, float v)
{
    for (int i = 0; i < im.w; i++)
        for (int j = 0; j < im.h; j++)
        {
            float color = get_pixel(im, i, j, c);
            set_pixel(im, i, j, c, color * v);
        }
}