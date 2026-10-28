"""Generate the project banner used at the top of the README.

Run with:  python docs/make_banner.py

Produces docs/banner.png (1280 x 400). The artwork is generated
programmatically - it is a designed banner, not a screenshot.
"""

import math
import os

from PIL import Image, ImageDraw, ImageFont

W, H = 1280, 400
FONT_DIR = "/usr/share/fonts/truetype/dejavu"
OUT = os.path.join(os.path.dirname(__file__), "banner.png")


def lerp(a, b, t):
    return tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


def main():
    img = Image.new("RGB", (W, H), "#0F2027")
    draw = ImageDraw.Draw(img)

    # Vertical gradient background.
    top = (15, 32, 39)
    bottom = (22, 61, 72)
    for y in range(H):
        draw.line([(0, y), (W, y)], fill=lerp(top, bottom, y / H))

    # DNA double helix motif on the right.
    amp, period = 74, 145
    strand_a, strand_b = [], []
    for x in range(905, W):
        t = (x - 905) / 1.0
        y1 = H / 2 + amp * math.sin(t / period * 2 * math.pi)
        y2 = H / 2 + amp * math.sin(t / period * 2 * math.pi + math.pi)
        strand_a.append((x, y1))
        strand_b.append((x, y2))

    # Rungs between the strands.
    for i in range(0, len(strand_a), 26):
        draw.line([strand_a[i], strand_b[i]], fill=(46, 125, 91), width=2)

    draw.line(strand_a, fill=(126, 217, 173), width=5, joint="curve")
    draw.line(strand_b, fill=(94, 160, 200), width=5, joint="curve")

    # Text block.
    title_font = ImageFont.truetype(os.path.join(FONT_DIR, "DejaVuSans-Bold.ttf"), 58)
    sub_font = ImageFont.truetype(os.path.join(FONT_DIR, "DejaVuSans.ttf"), 26)
    tag_font = ImageFont.truetype(os.path.join(FONT_DIR, "DejaVuSans.ttf"), 20)

    draw.text((70, 120), "DNA Sequence", font=title_font, fill=(255, 255, 255))
    draw.text((70, 185), "Analysis Toolkit", font=title_font, fill=(126, 217, 173))
    draw.text(
        (72, 268),
        "Composition · Reverse complement · Translation · ORFs",
        font=sub_font,
        fill=(206, 224, 232),
    )
    draw.text(
        (72, 310),
        "Pure Python · command line and Streamlit interface · MIT licensed",
        font=tag_font,
        fill=(150, 180, 195),
    )

    img.save(OUT, "PNG")
    print("wrote", OUT, img.size)


if __name__ == "__main__":
    main()
