#!/opt/epd-7.2-2-rh5-x86_64/bin/python2.7

from PIL import Image
import sys

if not len(sys.argv) > 3:
    raise SystemExit("Usage: %s src1 [src2] .. dest" % sys.argv[0])

images = map(Image.open, sys.argv[1:-1])
w = sum(i.size[0] for i in images)
mh = max(i.size[1] for i in images)

result = Image.new("RGBA", (w, mh))

x = 0
for i in images:
    result.paste(i, (x, 0))
    x += i.size[0]

result.save(sys.argv[-1])

