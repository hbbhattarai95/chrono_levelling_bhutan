"""Map helpers: Bhutan National Grid (TM, 90E, k0 = 1, FE 250000) and boundary reading."""
import json, struct
import numpy as np


def tm_forward(lon, lat, lon0=90.0, k0=1.0, fe=250000.0, a=6378137.0, f=1 / 298.257222101):
    e2 = f * (2 - f); ep2 = e2 / (1 - e2)
    phi = np.radians(lat); dl = np.radians(np.asarray(lon) - lon0)
    N = a / np.sqrt(1 - e2 * np.sin(phi) ** 2); T = np.tan(phi) ** 2
    Cc = ep2 * np.cos(phi) ** 2; A = dl * np.cos(phi)
    M = a * ((1 - e2 / 4 - 3 * e2 ** 2 / 64 - 5 * e2 ** 3 / 256) * phi
             - (3 * e2 / 8 + 3 * e2 ** 2 / 32 + 45 * e2 ** 3 / 1024) * np.sin(2 * phi)
             + (15 * e2 ** 2 / 256 + 45 * e2 ** 3 / 1024) * np.sin(4 * phi)
             - (35 * e2 ** 3 / 3072) * np.sin(6 * phi))
    x = k0 * N * (A + (1 - T + Cc) * A ** 3 / 6 + (5 - 18 * T + T ** 2 + 72 * Cc - 58 * ep2) * A ** 5 / 120) + fe
    y = k0 * (M + N * np.tan(phi) * (A ** 2 / 2 + (5 - T + 9 * Cc + 4 * Cc ** 2) * A ** 4 / 24
                                    + (61 - 58 * T + T ** 2 + 600 * Cc - 330 * ep2) * A ** 6 / 720))
    return x, y


def read_boundary(path):
    """Read a polygon boundary (.shp or .geojson) and return rings in Bhutan National Grid (m)."""
    rings = []
    if path.lower().endswith(".shp"):
        d = open(path, "rb").read(); pos = 100
        while pos < len(d):
            _, clen = struct.unpack(">ii", d[pos:pos + 8]); rec = d[pos + 8:pos + 8 + clen * 2]; pos += 8 + clen * 2
            if struct.unpack("<i", rec[:4])[0] not in (3, 5, 13, 15, 23, 25):
                continue
            npart, npts = struct.unpack("<ii", rec[36:44])
            parts = list(struct.unpack(f"<{npart}i", rec[44:44 + 4 * npart])) + [npts]; o = 44 + 4 * npart
            pts = np.frombuffer(rec[o:o + 16 * npts], dtype="<f8").reshape(npts, 2)
            rings += [pts[parts[k]:parts[k + 1]] for k in range(npart)]
    else:
        g = json.load(open(path))
        def geom(Gm):
            if Gm["type"] == "Polygon": rings.extend(np.array(r)[:, :2] for r in Gm["coordinates"])
            elif Gm["type"] == "MultiPolygon":
                for P in Gm["coordinates"]: rings.extend(np.array(r)[:, :2] for r in P)
            elif Gm["type"] == "GeometryCollection":
                for x in Gm["geometries"]: geom(x)
        if g["type"] == "FeatureCollection":
            for ft in g["features"]: geom(ft["geometry"])
        elif g["type"] == "Feature": geom(g["geometry"])
        else: geom(g)
    out = []
    for r in rings:
        if np.abs(r[:, 0]).max() <= 360 and np.abs(r[:, 1]).max() <= 90:
            x, y = tm_forward(r[:, 0], r[:, 1]); out.append(np.c_[x, y])
        else:
            out.append(r)
    return out
