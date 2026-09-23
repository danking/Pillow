from __future__ import annotations

import math
import random

import pytest

from PIL import Image, ImageFilter

from .helper import assert_image_equal, hopper

MODES = (
    "L",
    "LA",
    "La",
    "I",
    "I;16",
    "I;16B",
    "I;16L",
    "I;16N",
    "RGB",
    "RGBA",
    "RGBa",
    "RGBX",
    "CMYK",
)


def medianfilter_reference(im: Image.Image, filter_size: int) -> Image.Image:
    size = im.size
    margin = filter_size // 2
    rank = filter_size * filter_size // 2
    expected_bands = []
    for band in im.split():
        expected = []
        for y in range(size[1]):
            for x in range(size[0]):
                window = [
                    band.getpixel(
                        (
                            min(max(x + dx, 0), size[0] - 1),
                            min(max(y + dy, 0), size[1] - 1),
                        )
                    )
                    for dy in range(-margin, margin + 1)
                    for dx in range(-margin, margin + 1)
                ]
                expected.append(sorted(window)[rank])
        expected_bands.append(Image.new("L", size))
        expected_bands[-1].putdata(expected)

    return Image.merge(im.mode, expected_bands)


def medianfilter_3x3_reference(im: Image.Image) -> Image.Image:
    return medianfilter_reference(im, 3)


def medianfilter_5x5_reference(im: Image.Image) -> Image.Image:
    return medianfilter_reference(im, 5)


@pytest.mark.parametrize(
    "filter_to_apply",
    (
        ImageFilter.BLUR,
        ImageFilter.CONTOUR,
        ImageFilter.DETAIL,
        ImageFilter.EDGE_ENHANCE,
        ImageFilter.EDGE_ENHANCE_MORE,
        ImageFilter.EMBOSS,
        ImageFilter.FIND_EDGES,
        ImageFilter.SMOOTH,
        ImageFilter.SMOOTH_MORE,
        ImageFilter.SHARPEN,
        ImageFilter.MaxFilter,
        ImageFilter.MedianFilter,
        ImageFilter.MinFilter,
        ImageFilter.ModeFilter,
        ImageFilter.GaussianBlur,
        ImageFilter.GaussianBlur(0),
        ImageFilter.GaussianBlur(5),
        ImageFilter.GaussianBlur((2, 5)),
        ImageFilter.BoxBlur(0),
        ImageFilter.BoxBlur(5),
        ImageFilter.BoxBlur((2, 5)),
        ImageFilter.UnsharpMask,
        ImageFilter.UnsharpMask(10),
    ),
)
@pytest.mark.parametrize("mode", MODES)
def test_sanity(
    filter_to_apply: ImageFilter.Filter | type[ImageFilter.Filter], mode: str
) -> None:
    im = hopper(mode)
    if mode[0] != "I" or (
        callable(filter_to_apply)
        and issubclass(filter_to_apply, ImageFilter.BuiltinFilter)
    ):
        out = im.filter(filter_to_apply)
        assert out.mode == im.mode
        assert out.size == im.size


@pytest.mark.parametrize("mode", MODES)
def test_sanity_error(mode: str) -> None:
    im = hopper(mode)
    with pytest.raises(TypeError):
        im.filter("hello")  # type: ignore[arg-type]


def test_noop_on_small_images() -> None:
    # If image is smaller than the kernel size, return it as-is.
    kernel_size: tuple[int, int] = ImageFilter.SMOOTH_MORE.filterargs[0]
    kernel_w, kernel_h = kernel_size
    for w in range(1, kernel_w):
        for h in range(1, kernel_h):
            im = hopper("RGB").resize((w, h))
            # Precondition for the below equality test:
            # filter is larger or equal to image.
            assert im.size < kernel_size
            assert_image_equal(im.filter(ImageFilter.SMOOTH_MORE), im)


@pytest.mark.parametrize(
    "mode, expected",
    (
        ("1", (4, 0)),
        ("L", (4, 0)),
        ("P", (4, 0)),
        ("RGB", ((4, 0, 0), (0, 0, 0))),
    ),
)
def test_modefilter(
    mode: str,
    expected: tuple[int, int] | tuple[tuple[int, int, int], tuple[int, int, int]],
) -> None:
    im = Image.new(mode, (3, 3), None)
    im.putdata(list(range(9)))
    # image is:
    #   0 1 2
    #   3 4 5
    #   6 7 8
    mod = im.filter(ImageFilter.ModeFilter).getpixel((1, 1))
    im.putdata([0, 0, 1, 2, 5, 1, 5, 2, 0])  # mode=0
    mod2 = im.filter(ImageFilter.ModeFilter).getpixel((1, 1))
    assert (mod, mod2) == expected


@pytest.mark.parametrize(
    "mode, expected",
    (
        ("1", (0, 4, 8)),
        ("L", (0, 4, 8)),
        ("RGB", ((0, 0, 0), (4, 0, 0), (8, 0, 0))),
        ("I", (0, 4, 8)),
        ("F", (0.0, 4.0, 8.0)),
    ),
)
def test_rankfilter(
    mode: str,
    expected: (
        tuple[float, float, float]
        | tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
    ),
) -> None:
    im = Image.new(mode, (3, 3), None)
    im.putdata(list(range(9)))
    # image is:
    #   0 1 2
    #   3 4 5
    #   6 7 8
    minimum = im.filter(ImageFilter.MinFilter).getpixel((1, 1))
    med = im.filter(ImageFilter.MedianFilter).getpixel((1, 1))
    maximum = im.filter(ImageFilter.MaxFilter).getpixel((1, 1))
    assert (minimum, med, maximum) == expected


@pytest.mark.parametrize(
    "filter", (ImageFilter.MinFilter, ImageFilter.MedianFilter, ImageFilter.MaxFilter)
)
def test_rankfilter_error(filter: ImageFilter.RankFilter) -> None:
    with pytest.raises(ValueError):
        im = Image.new("P", (3, 3), None)
        im.putdata(list(range(9)))
        # image is:
        #   0 1 2
        #   3 4 5
        #   6 7 8
        im.filter(filter).getpixel((1, 1))


def test_rankfilter_properties() -> None:
    rankfilter = ImageFilter.RankFilter(3, 2)

    assert rankfilter.size == 3
    assert rankfilter.rank == 2

    with pytest.raises(ValueError, match="bad filter size"):
        ImageFilter.RankFilter(2, 1)
    with pytest.raises(ValueError, match="bad filter size"):
        ImageFilter.MaxFilter(2)
    with pytest.raises(ValueError, match="bad filter size"):
        ImageFilter.MedianFilter(2)
    with pytest.raises(ValueError, match="bad filter size"):
        ImageFilter.MinFilter(2)

    with pytest.raises(ValueError, match="filter size too large"):
        ImageFilter.RankFilter(23171, 1)
    im = Image.new("1", (1, 1))
    with pytest.raises(ValueError, match="filter size too large"):
        im.im.expand(23171)

    with pytest.raises(ValueError, match="bad rank value"):
        ImageFilter.RankFilter(1, 1)


def test_rankfilter_overflow() -> None:
    # Large margins used to overflow the ImagingExpand overflow guard itself (SIGFPE),
    # by mutating RankFilter.size after construction, bypassing __init__'s validation.
    im = Image.new("L", (16, 16))
    rankfilter = ImageFilter.RankFilter(3, 0)

    for size in (2**31, 2**32 - 1):  # margins of 2**30 and INT_MAX
        rankfilter.size = size
        with pytest.raises(ValueError, match="filter size too large"):
            im.filter(rankfilter)


EXPAND_MODES = (
    "1",
    "L",
    "P",
    "LA",
    "La",
    "PA",
    "I",
    "F",
    "I;16",
    "I;16B",
    "I;16L",
    "I;16N",
    "RGB",
    "RGBA",
    "RGBa",
    "RGBX",
    "CMYK",
    "YCbCr",
    "LAB",
    "HSV",
)


def expand_test_image(mode: str, size: tuple[int, int]) -> Image.Image:
    im = Image.new(mode, size)
    bands = Image.getmodebands(mode)
    values: list[float | int | tuple[int, ...]] = []
    for i in range(size[0] * size[1]):
        if mode == "F":
            values.append((i - 7) / 3)
        elif mode == "I":
            values.append(i * 100003 - 200000)
        elif mode == "1":
            values.append(255 if i % 3 else 0)
        elif bands == 1:
            values.append((i * 977 + 31) % 65536)
        else:
            values.append(
                tuple((i * 73 + band * 41 + 17) % 256 for band in range(bands))
            )
    im.putdata(values)
    if mode == "P":
        im.putpalette([value % 256 for value in range(768)])
    return im


@pytest.mark.parametrize("mode", EXPAND_MODES)
def test_expand_matches_clamped_border_oracle(mode: str) -> None:
    rng = random.Random(f"expand-{mode}")
    dimensions = [(1, 1), (2, 3), (15, 2), (16, 5), (17, 3), (31, 4)]
    dimensions.extend((rng.randrange(1, 34), rng.randrange(1, 9)) for _ in range(10))

    for width, height in dimensions:
        im = expand_test_image(mode, (width, height))
        for margin in (0, 1, 2, 3):
            expanded = Image.Image()._new(im.im.expand(margin))
            assert expanded.mode == mode
            assert expanded.size == (width + 2 * margin, height + 2 * margin)
            for y in range(expanded.height):
                source_y = min(max(y - margin, 0), height - 1)
                for x in range(expanded.width):
                    source_x = min(max(x - margin, 0), width - 1)
                    assert expanded.getpixel((x, y)) == im.getpixel(
                        (source_x, source_y)
                    )
            if mode == "P":
                assert expanded.getpalette() == im.getpalette()


@pytest.mark.parametrize("mode", ("L", "RGB", "I", "F"))
@pytest.mark.parametrize("size", ((0, 0), (0, 3), (3, 0)))
def test_expand_zero_margin_empty_image(mode: str, size: tuple[int, int]) -> None:
    im = Image.new(mode, size)
    expanded = Image.Image()._new(im.im.expand(0))

    assert expanded.mode == im.mode
    assert expanded.size == im.size
    assert expanded.im is not im.im


def test_expand_returns_independent_copy() -> None:
    im = expand_test_image("L", (3, 2))
    original = im.tobytes()

    expanded = Image.Image()._new(im.im.expand(2))
    expanded.putpixel((2, 2), 255)

    assert im.tobytes() == original


@pytest.mark.parametrize(
    "mode", ("L", "LA", "La", "RGB", "RGBA", "RGBa", "RGBX", "CMYK")
)
@pytest.mark.parametrize("pattern", ("random", "constant", "low_cardinality", "duplicates"))
def test_medianfilter_3x3_uint8_patterns(mode: str, pattern: str) -> None:
    rng = random.Random(8675309)
    bands = Image.getmodebands(mode)
    size = (17, 13)

    if pattern == "random":
        values = [rng.randrange(256) for _ in range(size[0] * size[1] * bands)]
    elif pattern == "constant":
        values = [73] * (size[0] * size[1] * bands)
    elif pattern == "low_cardinality":
        palette = (0, 17, 17, 128, 240, 255)
        values = [
            palette[rng.randrange(len(palette))]
            for _ in range(size[0] * size[1] * bands)
        ]
    else:
        values = [(i // bands) % 3 * 85 for i in range(size[0] * size[1] * bands)]

    im = Image.frombytes(mode, size, bytes(values))
    result = im.filter(ImageFilter.MedianFilter(3))

    assert_image_equal(result, medianfilter_3x3_reference(im))


@pytest.mark.parametrize(
    "mode", ("L", "LA", "La", "RGB", "RGBA", "RGBa", "RGBX", "CMYK")
)
@pytest.mark.parametrize(
    "pattern", ("random", "constant", "low_cardinality", "duplicates", "extremes")
)
@pytest.mark.parametrize("size", ((17, 13), (21, 9), (18, 11)))
def test_medianfilter_5x5_uint8_patterns(
    mode: str, pattern: str, size: tuple[int, int]
) -> None:
    rng = random.Random(8675309)
    bands = Image.getmodebands(mode)

    if pattern == "random":
        values = [rng.randrange(256) for _ in range(size[0] * size[1] * bands)]
    elif pattern == "constant":
        values = [73] * (size[0] * size[1] * bands)
    elif pattern == "low_cardinality":
        palette = (0, 17, 17, 128, 240, 255)
        values = [
            palette[rng.randrange(len(palette))]
            for _ in range(size[0] * size[1] * bands)
        ]
    elif pattern == "duplicates":
        values = [(i // bands) % 5 * 51 for i in range(size[0] * size[1] * bands)]
    else:
        values = [
            (0, 255, 1, 254, 128)[(i // bands) % 5]
            for i in range(size[0] * size[1] * bands)
        ]

    im = Image.frombytes(mode, size, bytes(values))
    result = im.filter(ImageFilter.MedianFilter(5))

    assert_image_equal(result, medianfilter_5x5_reference(im))


@pytest.mark.parametrize("mode", ("L", "RGB", "RGBA"))
@pytest.mark.parametrize("pattern", ("random", "constant", "low_cardinality"))
@pytest.mark.parametrize("size", ((64, 7), (65, 9), (79, 6)))
def test_medianfilter_5x5_uint8_wide_patterns(
    mode: str, pattern: str, size: tuple[int, int]
) -> None:
    rng = random.Random(8675309)
    count = size[0] * size[1] * Image.getmodebands(mode)
    values = [rng.randrange(256) for _ in range(count)]
    if pattern == "constant":
        values = [73] * count
    elif pattern == "low_cardinality":
        palette = (0, 17, 128, 240, 255)
        values = [palette[value % len(palette)] for value in values]

    im = Image.frombytes(mode, size, bytes(values))

    assert_image_equal(
        im.filter(ImageFilter.MedianFilter(5)), medianfilter_5x5_reference(im)
    )


@pytest.mark.parametrize("mode", ("L", "RGB"))
@pytest.mark.parametrize("size", ((1, 1), (1, 5), (5, 1), (2, 2)))
def test_medianfilter_3x3_tiny_uint8_images(mode: str, size: tuple[int, int]) -> None:
    im = hopper(mode).resize(size)
    assert_image_equal(
        im.filter(ImageFilter.MedianFilter(3)), medianfilter_3x3_reference(im)
    )


@pytest.mark.parametrize("mode", ("L", "RGB"))
@pytest.mark.parametrize("size", ((1, 1), (1, 5), (5, 1), (2, 2), (4, 7), (7, 4)))
def test_medianfilter_5x5_tiny_uint8_images(mode: str, size: tuple[int, int]) -> None:
    im = hopper(mode).resize(size)
    assert_image_equal(
        im.filter(ImageFilter.MedianFilter(5)), medianfilter_5x5_reference(im)
    )


@pytest.mark.parametrize("mode", ("I", "F"))
def test_medianfilter_3x3_numeric_modes(mode: str) -> None:
    im = Image.new(mode, (3, 3))
    im.putdata([5, 1, 5, 7, 3, 3, 9, 1, 5])
    assert im.filter(ImageFilter.MedianFilter(3)).getpixel((1, 1)) == 5


@pytest.mark.parametrize("mode", ("I", "F"))
def test_medianfilter_5x5_numeric_modes(mode: str) -> None:
    values = [
        17,
        5,
        5,
        2,
        9,
        1,
        21,
        21,
        3,
        3,
        8,
        13,
        13,
        13,
        34,
        55,
        1,
        2,
        89,
        144,
        0,
        0,
        233,
        377,
        610,
    ]
    im = Image.new(mode, (5, 5))
    im.putdata(values)
    assert im.filter(ImageFilter.MedianFilter(5)).getpixel((2, 2)) == sorted(values)[12]


def test_builtinfilter_p() -> None:
    builtin_filter = ImageFilter.BuiltinFilter()

    with pytest.raises(ValueError):
        builtin_filter.filter(hopper("P").im)


def test_kernel_not_enough_coefficients() -> None:
    with pytest.raises(ValueError):
        ImageFilter.Kernel((3, 3), (0, 0))


EMBOSS_MATRIX = {
    3: (
        -1, -1,  0,
        -1,  0,  1,
         0,  1,  1,
    ),
    5: (
        -1, -1, -1, -1,  0,
        -1, -1, -1,  0,  1,
        -1, -1,  0,  1,  1,
        -1,  0,  1,  1,  1,
         0,  1,  1,  1,  1,
    ),
}  # fmt: skip


@pytest.mark.parametrize("size", (3, 5))
def test_consistency(size: int) -> None:
    kernel = ImageFilter.Kernel((size, size), EMBOSS_MATRIX[size], 0.3)
    with Image.open("Tests/images/hopper.bmp") as source:
        with Image.open(f"Tests/images/hopper_emboss_{size}x{size}.bmp") as reference:
            assert_image_equal(source.filter(kernel), reference)


@pytest.mark.parametrize("size", (3, 5))
@pytest.mark.parametrize("mode", ("I;16", "I;16L", "I;16B", "I;16N"))
def test_consistency_i16(size: int, mode: str) -> None:
    kernel = ImageFilter.Kernel((size, size), EMBOSS_MATRIX[size], 0.3)
    reference = hopper("I").filter(kernel)
    result = hopper(mode).filter(kernel)
    assert_image_equal(result.convert("I"), reference)


@pytest.mark.parametrize("mode", ("I;16", "I;16L", "I;16B", "I;16N"))
def test_consistency_i16_high_byte(mode: str) -> None:
    # Exercise filters with a 16bpc image that has content in the high byte, too.
    im = Image.new(mode, (8, 8), 1000)
    # Ensure repeated smoothing retains the exact same color value,
    # rather than drifting due to rounding errors.
    for _ in range(5):
        im = im.filter(ImageFilter.SMOOTH)
        assert im.getpixel((4, 4)) == 1000


@pytest.mark.parametrize("size", (-1, math.nan, math.inf, 2**31))
def test_invalid_box_blur_filter(size: int) -> None:
    for radius in (size, (size, size), (size, 1), (1, size)):
        with pytest.raises(ValueError, match="radius"):
            ImageFilter.BoxBlur(radius)

        im = hopper()
        box_blur_filter = ImageFilter.BoxBlur(2)
        box_blur_filter.radius = radius
        with pytest.raises(ValueError, match="radius"):
            im.filter(box_blur_filter)


@pytest.mark.parametrize("radius", (math.nan, math.inf, 2**31))
def test_invalid_gaussian_blur_filter(radius: int) -> None:
    im = hopper()
    with pytest.raises(ValueError, match="radius"):
        im.filter(ImageFilter.GaussianBlur(radius))


def test_rankfilter_size_1() -> None:
    im = Image.new("L", (3, 3), 128)

    # Size 1 should not crash (margin is 0)
    assert im.filter(ImageFilter.MinFilter(1)).getpixel((1, 1)) == 128
    assert im.filter(ImageFilter.MaxFilter(1)).getpixel((1, 1)) == 128
    assert im.filter(ImageFilter.MedianFilter(1)).getpixel((1, 1)) == 128
    assert im.filter(ImageFilter.RankFilter(1, 0)).getpixel((1, 1)) == 128
