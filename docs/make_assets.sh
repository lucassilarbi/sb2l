#!/usr/bin/env bash
#
# Regenerate every image of the README from the editor itself.
#
#   ./docs/make_assets.sh [path/to/sb2l_gui]
#
# Needs ImageMagick and gifsicle; a display is required, so the script runs
# the editor under xvfb-run when there is none. Output goes to docs/img/.
set -euo pipefail

GUI=${1:-build/gui/sb2l_gui}
OUT=$(cd "$(dirname "$0")" && pwd)/img
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"

if [ ! -x "$GUI" ]; then
    echo "editor not found: $GUI (build it with -DSB2L_BUILD_GUI=ON)" >&2
    exit 1
fi

# The editor needs a GL context; use the display if there is one, else xvfb.
if [ -n "${DISPLAY:-}" ]; then RUN=(); else RUN=(xvfb-run -a -s "-screen 0 1280x800x24"); fi
shot() { "${RUN[@]}" "$GUI" "$@" >/dev/null; }

# The controls panel occupies the left ~312 px; `canvas` keeps the drawing only.
canvas() { convert "$1" -crop 968x800+312+0 +repage -resize "${3:-820}x" "$2"; }

# Every single-configuration figure uses the affine parameter and the Taylor
# form, the pair that gives the narrowest enclosure. The 3x3 matrix below is
# the exception: the parameter set is what it varies.
echo "== feature captures =="
shot --shot "$TMP/2d_boxes.ppm"      --cs IR --ps Z
shot --shot "$TMP/2d_zono.ppm"       --cs Z  --ps Z  --gen 3
# Rational weights are shown on real control points: the quotient basis widens
# the enclosure sharply when it is combined with control boxes (see the note
# in the README), and the weight itself is what this figure is about.
shot --shot "$TMP/rational.ppm"      --cs R  --ps R  --ct CR
shot --shot "$TMP/3d_boxes.ppm"      --cs IR --ps Z  --dim 3
# Four generators and enlarged control sets: with the sets at their editing
# size the result elements are legible volumes rather than thin lenses, and
# the facet mesh of a control zonotope is actually visible.
shot --shot "$TMP/3d_zono.ppm"       --cs Z  --ps Z  --dim 3 --gen 4 --n 6 --d 10 --size 1.4

# Full window (panel included) for the two overview figures, canvas only for the rest.
convert "$TMP/2d_boxes.ppm" -resize 1000x "$OUT/gui_2d.png"
convert "$TMP/3d_boxes.ppm" -resize 1000x "$OUT/gui_3d.png"
canvas "$TMP/2d_zono.ppm"  "$OUT/zonotopes_2d.png"
canvas "$TMP/rational.ppm" "$OUT/rational.png"
canvas "$TMP/3d_zono.ppm"  "$OUT/zonotopes_3d.png"

echo "== 3x3 matrix: curve parameter x control points =="
for ps in R IR Z; do
    for cs in R IR Z; do
        shot --shot "$TMP/m_${ps}_${cs}.ppm" --ps "$ps" --cs "$cs" --gen 3
        canvas "$TMP/m_${ps}_${cs}.ppm" "$TMP/m_${ps}_${cs}.png" 420
    done
done
montage -label 'u:%[label]' \
    -background '#1a1a1c' -fill '#d0d0d0' -pointsize 20 -geometry +4+4 -tile 3x3 \
    -label 'u in R,  P in R'   "$TMP/m_R_R.png"   -label 'u in R,  P in IR'  "$TMP/m_R_IR.png" \
    -label 'u in R,  P in Z'   "$TMP/m_R_Z.png"   -label 'u in IR, P in R'   "$TMP/m_IR_R.png" \
    -label 'u in IR, P in IR'  "$TMP/m_IR_IR.png" -label 'u in IR, P in Z'   "$TMP/m_IR_Z.png" \
    -label 'u in Z,  P in R'   "$TMP/m_Z_R.png"   -label 'u in Z,  P in IR'  "$TMP/m_Z_IR.png" \
    -label 'u in Z,  P in Z'   "$TMP/m_Z_Z.png"   "$OUT/matrix.png"

echo "== animated tour =="
shot --tour "$TMP"
for f in "$TMP"/f0*.ppm; do convert "$f" -resize 760x "${f%.ppm}.png"; done
convert -delay 6 -loop 0 "$TMP"/f0*.png "$TMP/raw.gif"
gifsicle --optimize=3 --colors 96 --lossy=40 "$TMP/raw.gif" > "$OUT/tour.gif"

echo "done -> $OUT"
ls -la "$OUT"
