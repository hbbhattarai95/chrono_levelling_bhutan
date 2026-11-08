"""Map utilities for Bhutan National Grid (Transverse Mercator, central meridian
90°E, scale factor 1, and false easting 250,000 m) and boundary-file reading.
"""

import json, struct
import numpy as np


def tm_forward(lon, lat, lon0=90.0, k0=1.0, fe=250000.0, a=6378137.0, f=1 / 298.257222101):
    """Convert longitude and latitude to projected coordinates using the
    Transverse Mercator projection.

    The default parameters define the Bhutan National Grid. The calculation
    uses the ellipsoidal Transverse Mercator series expansion.

    Parameters
    ----------
    lon, lat : array-like or float
        Geographic longitude and latitude in decimal degrees.
    lon0 : float
        Longitude of the central meridian, in degrees.
    k0 : float
        Scale factor at the central meridian.
    fe : float
        False easting, in metres.
    a : float
        Semi-major axis of the reference ellipsoid, in metres.
    f : float
        Flattening of the reference ellipsoid.

    Returns
    -------
    x, y : array-like or float
        Projected easting and northing coordinates, in metres.
    """
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
    """Read polygon boundaries from a shapefile or GeoJSON file.

    Extract the polygon rings and return their coordinates in Bhutan
    National Grid (metres). Geographic coordinates are projected using
    the Transverse Mercator transformation defined in tm_forward().
    Coordinates that are already projected are retained unchanged.

    Parameters
    ----------
    path : str
        Path to the input .shp or .geojson boundary file.

    Returns
    -------
    out : list of numpy.ndarray
        A list containing the coordinate arrays for each polygon ring.
        Each array contains projected easting and northing coordinates
        in metres.
    """
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