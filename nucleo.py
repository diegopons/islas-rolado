# -*- coding: utf-8 -*-
"""Nucleo de calculo de la aplicacion "Islas de conservacion y plan de rolado" (sin interfaz).
- lectura del CHM (GeoTIFF) y del lote (KML / KMZ / GeoJSON), llevados a UTM a 1 m
- clasificacion por altura, islas de conservacion (metodo MBGI potreros 4 y 6), estadisticos y plan de rolado
"""
import io, json, zipfile, re
import xml.etree.ElementTree as ET
import math
import numpy as np
from scipy import ndimage as ndi
from shapely import wkb as _wkb
from shapely.geometry import Polygon
from shapely.ops import unary_union

S8 = np.ones((3, 3))
CLASES = {1: "Suelo desnudo (<0,10 m)", 2: "Pastos (0,10-0,50 m)", 3: "Arbustal", 4: "Arboles"}


def _partes(g, amin=0.0):
    if g is None or g.is_empty: return []
    return [q for q in (g.geoms if hasattr(g, "geoms") else [g]) if q.geom_type == "Polygon" and q.area > amin]


def _cierre(m, r):
    """Cierre morfologico con disco de radio r (celdas), por transformada de distancia."""
    if not m.any(): return m.copy()
    dil = ndi.distance_transform_edt(~m) <= r
    return ndi.distance_transform_edt(dil) > r


def _apertura(m, r):
    if not m.any(): return m.copy()
    ero = ndi.distance_transform_edt(m) > r
    return ndi.distance_transform_edt(~ero) <= r if ero.any() else np.zeros_like(m)


def disenar(h, en, lote_geom, poligonizar, rasterizar, P, base_geom=None, log=lambda s: None):
    """Nucleo del metodo (sin QGIS). h: CHM a 1 m (nan = sin dato); en: mascara del lote; lote_geom: poligono (shapely);
    poligonizar(mask)->geometria shapely; rasterizar(geom)->mask. Devuelve dict con islas finales, priorizadas y tablas."""
    T = P["T"]; hh = np.nan_to_num(h, nan=0.0); enf = en & np.isfinite(h); A_lote = lote_geom.area
    c = np.select([hh > T, hh >= 0.5, hh >= 0.1], [4, 3, 2], 1)
    # 1-2. parches de arbol e islas por proximidad
    tr = (hh > T) & en; lab, n = ndi.label(tr, structure=S8); ar = np.bincount(lab.ravel(), minlength=n + 1)
    tr = np.isin(lab, [i for i in range(1, n + 1) if ar[i] >= P["AMIN_ARBOL"]])
    if not tr.any(): raise ValueError("No hay arboles de mas de %.1f m dentro del lote." % T)
    isl = _cierre(tr, P["DIST"] / 2.0) & en
    # 3. relleno por isla (huecos + entrantes), sin unir islas separadas salvo que se toquen
    li, ni = ndi.label(isl, structure=S8); rell = np.zeros_like(isl); R = P["RELLENO"]; pad = int(math.ceil(R)) + 3
    for i, sl in enumerate(ndi.find_objects(li), 1):
        r0, r1 = max(sl[0].start - pad, 0), min(sl[0].stop + pad, isl.shape[0]); c0, c1 = max(sl[1].start - pad, 0), min(sl[1].stop + pad, isl.shape[1])
        g = li[r0:r1, c0:c1] == i
        if R > 0:
            cc = ndi.binary_fill_holes(_cierre(g, R)); a = _apertura(cc, R / 3.0); g2 = _cierre(a | g, 3.0)
        else: g2 = g
        rell[r0:r1, c0:c1] |= g2
    rell = ndi.binary_fill_holes(rell) & en
    # 4. priorizacion por radio equivalente
    lr, nr = ndi.label(rell, structure=S8); ar = np.bincount(lr.ravel(), minlength=nr + 1); amin = math.pi * P["RADIO"] ** 2
    prio = np.isin(lr, [i for i in range(1, nr + 1) if ar[i] >= amin])
    log("Parches de arbol: %d | islas: %d | priorizadas: %d (%.2f ha, %.1f %% del lote)" % (int(ndi.label(tr, structure=S8)[1]), nr, int((ar[1:] >= amin).sum()), prio.sum() / 1e4, 100.0 * prio.sum() / max(en.sum(), 1)))
    base = (rasterizar(base_geom) & en) if base_geom is not None else prio
    if not base.any(): raise ValueError("No quedaron islas priorizadas: baje el radio minimo o aumente la distancia de union.")
    # 5. campo de prioridad y umbral que da la meta
    s = P["SIG_ARB"]; dens = ndi.gaussian_filter(tr.astype(float), s) / np.maximum(ndi.gaussian_filter(en.astype(float), s), 1e-3)
    F = ndi.distance_transform_edt(base) - ndi.distance_transform_edt(~base) + P["PESO_ARB"] * dens
    AP = P["A_MIN_PIEZA"]

    def construir(t):
        m = (F >= t) & en
        m = (ndi.gaussian_filter(m.astype(float), P["SIG_SUAV"]) >= 0.5) & en
        l2, n2 = ndi.label(m, structure=S8); a2 = np.bincount(l2.ravel(), minlength=n2 + 1); m = np.isin(l2, [i for i in range(1, n2 + 1) if a2[i] >= AP])
        if not m.any(): return [], 0.0
        l3, n3 = ndi.label(~m & en); a3 = np.bincount(l3.ravel(), minlength=n3 + 1); m |= np.isin(l3, [i for i in range(1, n3 + 1) if a3[i] < AP]) & en
        pg = poligonizar(m)                                                    # 6. simplificacion
        pg = pg.buffer(8, 8).buffer(-4, 8).simplify(P["TOL"]).buffer(0).intersection(lote_geom)
        ps = [Polygon(q.exterior, [r for r in q.interiors if Polygon(r).area >= AP]) for q in _partes(pg, AP)]
        return ps, sum(q.area for q in ps) / A_lote

    lo, hi = float(np.nanmin(F[en])) - 1.0, float(np.nanmax(F[en])) + 1.0; mejor = ([], 0.0)
    for _ in range(40):
        t = (lo + hi) / 2.0; ps, fr = construir(t)
        if abs(fr - P["META"]) < abs(mejor[1] - P["META"]): mejor = (ps, fr)
        if abs(fr - P["META"]) < 0.002 or hi - lo < 0.01: break
        if fr > P["META"]: lo = t
        else: hi = t
    ps, fr = mejor
    ps.sort(key=lambda q: -q.area)
    log("Islas finales: %d | %.2f ha = %.1f %% del lote (meta %.0f %%)" % (len(ps), sum(q.area for q in ps) / 1e4, 100 * fr, 100 * P["META"]))
    # estadisticos
    islas = []; fin = np.zeros_like(en)
    for i, q in enumerate(ps, 1):
        ra = rasterizar(q) & enf; fin |= ra
        islas.append(dict(geom=q, n=i, superficie_ha=round(q.area / 1e4, 3), pct_del_lote=round(100 * q.area / A_lote, 1), vertices=len(q.exterior.coords) - 1,
                          arboles_ha=round(float((ra & (c == 4)).sum()) / 1e4, 3), arbustal_ha=round(float((ra & (c == 3)).sum()) / 1e4, 3),
                          pastos_ha=round(float((ra & (c == 2)).sum()) / 1e4, 3), suelo_ha=round(float((ra & (c == 1)).sum()) / 1e4, 3)))
    tabla = []; I = max(int(fin.sum()), 1); tot = max(int(enf.sum()), 1)
    for k, nm in CLASES.items():
        tt = int(((c == k) & enf).sum()); ii = int(((c == k) & fin).sum())
        tabla.append(dict(clase=nm, lote_ha=round(tt / 1e4, 2), lote_pct=round(100.0 * tt / tot, 1), en_islas_ha=round(ii / 1e4, 2), pct_de_las_islas=round(100.0 * ii / I, 1),
                          pct_de_la_clase_en_islas=round(100.0 * ii / max(tt, 1), 1), fuera_ha=round((tt - ii) / 1e4, 2)))
    tabla.append(dict(clase="TOTAL", lote_ha=round(tot / 1e4, 2), lote_pct=100.0, en_islas_ha=round(int(fin.sum()) / 1e4, 2), pct_de_las_islas=100.0,
                      pct_de_la_clase_en_islas=round(100.0 * int(fin.sum()) / tot, 1), fuera_ha=round((tot - int(fin.sum())) / 1e4, 2)))
    prio_geoms = []
    if prio.any():
        for q in sorted(_partes(poligonizar(prio).simplify(0.5).intersection(lote_geom), 1.0), key=lambda q: -q.area):
            prio_geoms.append(dict(geom=q, superficie_ha=round(q.area / 1e4, 3), radio_equivalente_m=round(math.sqrt(q.area / math.pi), 1)))
    return dict(islas=islas, tabla=tabla, priorizadas=prio_geoms, pct=100 * fr)


# =========================== entrada / salida ===========================
import rasterio
from rasterio.io import MemoryFile
from rasterio.warp import reproject, Resampling, transform_geom
from rasterio.features import shapes as _shapes, rasterize as _rasterize
from rasterio.transform import from_origin
from shapely.geometry import shape, mapping, MultiPolygon

PARAMS = dict(META=0.50, T=4.0, DIST=15.0, RADIO=15.0, RELLENO=20.0, PESO_ARB=60.0, TOL=6.0, A_MIN_PIEZA=5000.0, AMIN_ARBOL=4.0, SIG_ARB=15.0, SIG_SUAV=8.0)


def leer_poligonos(nombre, datos):
    """Devuelve [(nombre, poligono shapely en WGS84)] desde KML, KMZ o GeoJSON."""
    ext = nombre.lower().rsplit(".", 1)[-1]
    if ext == "kmz":
        z = zipfile.ZipFile(io.BytesIO(datos)); k = [n for n in z.namelist() if n.lower().endswith(".kml")]
        if not k: raise ValueError("El KMZ no contiene ningun KML.")
        datos = z.read(k[0]); ext = "kml"
    out = []
    if ext == "kml":
        txt = re.sub(rb'\sxmlns(:\w+)?="[^"]+"', b"", datos)                      # sin espacios de nombres
        txt = re.sub(rb"<(/?)\w+:", rb"<\1", txt)
        root = ET.fromstring(txt)
        for pm in root.iter("Placemark"):
            nm = (pm.findtext("name") or "").strip()
            for pol in pm.iter("Polygon"):
                c = pol.find("outerBoundaryIs/LinearRing/coordinates")
                if c is None or not (c.text or "").strip(): continue
                pts = [tuple(float(v) for v in t.split(",")[:2]) for t in c.text.split()]
                if len(pts) >= 4:
                    g = Polygon(pts).buffer(0)
                    if not g.is_empty: out.append((nm or "Poligono %d" % (len(out) + 1), g))
    elif ext in ("geojson", "json"):
        j = json.loads(datos.decode("utf-8-sig")); fs = j["features"] if j.get("type") == "FeatureCollection" else [j]
        for f in fs:
            g = shape(f.get("geometry", f)).buffer(0); pr = f.get("properties") or {}
            nm = str(pr.get("nombre") or pr.get("Name") or pr.get("name") or pr.get("id") or "Poligono %d" % (len(out) + 1))
            for q in (g.geoms if isinstance(g, MultiPolygon) else [g]):
                if q.geom_type == "Polygon" and not q.is_empty: out.append((nm, q))
    else:
        raise ValueError("Formato de lote no admitido (use KML, KMZ o GeoJSON en WGS84).")
    if not out: raise ValueError("No se encontraron poligonos en el archivo del lote.")
    return out


def _utm(lon, lat):
    return "EPSG:%d" % ((32700 if lat < 0 else 32600) + int((lon + 180) // 6) + 1)


def preparar(chm, lote_wgs=None):
    """chm: ruta o bytes de un GeoTIFF. Devuelve dict(h, en, transform, crs, lote(UTM), entero) en una grilla UTM de 1 m."""
    fuentes = list(chm) if isinstance(chm, (list, tuple)) else [chm]
    ds = MemoryFile(fuentes[0]).open() if isinstance(fuentes[0], (bytes, bytearray)) else rasterio.open(fuentes[0])
    with ds:
        if ds.crs is None: raise ValueError("El CHM no tiene sistema de coordenadas.")
        if lote_wgs is None:                                   # sin lote: se usa la extension con dato del CHM
            b = ds.bounds; g = transform_geom(ds.crs, "EPSG:4326", mapping(shape({"type": "Polygon", "coordinates": [[(b.left, b.bottom), (b.right, b.bottom), (b.right, b.top), (b.left, b.top), (b.left, b.bottom)]]})))
            lote_wgs = shape(g); sin_lote = True
        else: sin_lote = False
        c = lote_wgs.centroid; crs = _utm(c.x, c.y)
        lote = shape(transform_geom("EPSG:4326", crs, mapping(lote_wgs))).buffer(0)
        x0, y0, x1, y1 = lote.bounds; m = 30 if not sin_lote else 0
        x0, y0, x1, y1 = math.floor(x0) - m, math.floor(y0) - m, math.ceil(x1) + m, math.ceil(y1) + m
        nx, ny = int(x1 - x0), int(y1 - y0)
        if nx * ny > 60e6: raise ValueError("El area es demasiado grande para 1 m (%.0f ha). Recorte el CHM o use un lote mas chico." % (nx * ny / 1e4))
        tr = from_origin(x0, y1, 1.0, 1.0); h = np.full((ny, nx), np.nan, "float32")
        entero = np.issubdtype(np.dtype(ds.dtypes[0]), np.integer)
        fino = max(abs(ds.res[0]), abs(ds.res[1])) < (0.7 if not ds.crs.is_geographic else 0.7 / 111000.0)
        metodo = Resampling.average if (fino and not entero) else Resampling.nearest
        reproject(rasterio.band(ds, 1), h, src_transform=ds.transform, src_crs=ds.crs, src_nodata=ds.nodata, dst_transform=tr, dst_crs=crs, dst_nodata=np.nan, resampling=metodo)
        nodata0 = ds.nodata
    for otra in fuentes[1:]:                                   # mas de una fuente (lote que cae en varios mosaicos): se completan los huecos
        with rasterio.open(otra) as d2:
            h2 = np.full((ny, nx), np.nan, "float32")
            reproject(rasterio.band(d2, 1), h2, src_transform=d2.transform, src_crs=d2.crs, src_nodata=d2.nodata, dst_transform=tr, dst_crs=crs, dst_nodata=np.nan, resampling=metodo)
            falta = ~np.isfinite(h); h[falta] = h2[falta]
    h = h.astype("float64"); h[(h < -100) | (h > 150)] = np.nan
    if entero and nodata0 is None: h[h == 255] = np.nan
    en = _rasterize([(lote, 1)], out_shape=h.shape, transform=tr, fill=0, dtype="uint8").astype(bool)
    if sin_lote:
        en &= np.isfinite(h)
        if not en.any(): raise ValueError("El CHM no tiene datos validos.")
        lote = unary_union([shape(g) for g, v in _shapes(en.astype("uint8"), mask=en, transform=tr)]).buffer(0)
    if not en.any() or not np.isfinite(h[en]).any(): raise ValueError("El CHM no cubre el lote indicado.")
    return dict(h=h, en=en, transform=tr, crs=crs, lote=lote, entero=bool(entero), cobertura=float(np.isfinite(h[en]).mean()))


def procesar(prep, P, log=lambda s: None):
    """Corre el diseno y arma todos los productos del tablero."""
    h, en, tr, lote = prep["h"], prep["en"], prep["transform"], prep["lote"]
    pol = lambda m: unary_union([shape(g) for g, v in _shapes(m.astype("uint8"), mask=m, transform=tr)]) if m.any() else Polygon()
    ras = lambda g: _rasterize([(g, 1)], out_shape=h.shape, transform=tr, fill=0, dtype="uint8").astype(bool)
    r = disenar(h, en, lote, pol, ras, P, log=log)
    hh = np.nan_to_num(h, nan=0.0); enf = en & np.isfinite(h)
    c = np.select([hh > P["T"], hh >= 0.5, hh >= 0.1], [4, 3, 2], 1).astype("uint8"); c[~enf] = 0
    isl = np.zeros(h.shape, bool)
    for it in r["islas"]: isl |= ras(it["geom"])
    isl &= enf
    arb = c == 4
    plan = np.zeros(h.shape, "uint8")            # 1 conservar (isla) | 2 arbol fuera de isla (se esquiva) | 3 rolar
    plan[enf & isl] = 1; plan[enf & ~isl & arb] = 2; plan[enf & ~isl & ~arb] = 3
    ha = lambda m: float(m.sum()) / 1e4
    tot = ha(enf)
    res = dict(lote_ha=tot, islas_ha=ha(plan == 1), arboles_fuera_ha=ha(plan == 2), rolar_ha=ha(plan == 3),
               arboles_ha=ha(arb & enf), arboles_en_islas_ha=ha(arb & isl), n_islas=len(r["islas"]),
               rolar_por_clase={k: ha((plan == 3) & (c == k)) for k in (1, 2, 3)},
               altura_media=float(np.nanmean(h[enf])), altura_p95=float(np.nanpercentile(h[enf], 95)), altura_max=float(np.nanmax(h[enf])))
    res["pct_islas"] = 100 * res["islas_ha"] / tot; res["pct_rolar"] = 100 * res["rolar_ha"] / tot
    res["pct_arboles_protegidos"] = 100 * res["arboles_en_islas_ha"] / res["arboles_ha"] if res["arboles_ha"] > 0 else 0.0
    return dict(prep, c=c, isl=isl, plan=plan, islas=r["islas"], priorizadas=r["priorizadas"], tabla=r["tabla"], res=res)


def a_wgs(geom, crs):
    return shape(transform_geom(crs, "EPSG:4326", mapping(geom)))


def exportar_geojson(R, nombre_lote):
    fs = []
    for it in R["islas"]:
        pr = {k: v for k, v in it.items() if k != "geom"}; pr["id"] = "I%02d" % it["n"]; pr["lote"] = nombre_lote
        fs.append({"type": "Feature", "properties": pr, "geometry": mapping(a_wgs(it["geom"], R["crs"]))})
    return json.dumps({"type": "FeatureCollection", "name": "islas_conservacion", "features": fs})


def exportar_kml(R, nombre_lote):
    def anillo(r): return " ".join("%.7f,%.7f,0" % (x, y) for x, y in r.coords)
    pm = ""
    for it in R["islas"]:
        g = a_wgs(it["geom"], R["crs"])
        inn = "".join("<innerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></innerBoundaryIs>" % anillo(r) for r in g.interiors)
        pm += ("<Placemark><name>I%02d</name><description>%s | %.2f ha (%.1f %% del lote) | arboles %.2f ha</description><styleUrl>#isla</styleUrl>"
               "<Polygon><tessellate>1</tessellate><outerBoundaryIs><LinearRing><coordinates>%s</coordinates></LinearRing></outerBoundaryIs>%s</Polygon></Placemark>\n"
               % (it["n"], nombre_lote, it["superficie_ha"], it["pct_del_lote"], it["arboles_ha"], anillo(g.exterior), inn))
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2"><Document><name>Islas de conservacion - %s</name>'
            '<Style id="isla"><LineStyle><color>ff00c800</color><width>3</width></LineStyle><PolyStyle><color>7000dc00</color></PolyStyle></Style>\n%s</Document></kml>\n' % (nombre_lote, pm))


# =========================== CHM satelital automatico (Meta / WRI, version 2) ===========================
import os
META_URL = os.environ.get("ISLAS_META_URL", "https://dataforgood-fb-data.s3.amazonaws.com/forests/v2/global/dinov3_global_chm_v2_ml3/chm/{qk}.tif")
META_ZOOM = 10
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR"); os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif")
os.environ.setdefault("GDAL_HTTP_TIMEOUT", "40"); os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "2")


def _tile(lon, lat, z):
    n = 2 ** z; lat = max(min(lat, 85.0511), -85.0511)
    return int((lon + 180.0) / 360.0 * n), int((1.0 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2.0 * n)


def _quadkey(x, y, z):
    q = ""
    for i in range(z, 0, -1):
        d = 0; m = 1 << (i - 1)
        if x & m: d += 1
        if y & m: d += 2
        q += str(d)
    return q


def fuentes_meta(lote_wgs):
    """Direcciones de los mosaicos de Meta/WRI v2 (1 m, global) que cubren el lote. No descarga nada: la lectura posterior
    trae solo la ventana del lote (GeoTIFF optimizado para la nube)."""
    x0, y0, x1, y1 = lote_wgs.bounds
    if (x1 - x0) > 0.6 or (y1 - y0) > 0.6: raise ValueError("El poligono es demasiado grande para trabajar a 1 m.")
    tx0, ty0 = _tile(x0, y1, META_ZOOM); tx1, ty1 = _tile(x1, y0, META_ZOOM)
    return [META_URL.format(qk=_quadkey(x, y, META_ZOOM)) for x in range(tx0, tx1 + 1) for y in range(ty0, ty1 + 1)]


def preparar_meta(lote_wgs):
    fs = fuentes_meta(lote_wgs); ok = []
    for f in fs:
        try:
            with rasterio.open(f): ok.append(f)
        except Exception:
            pass
    if not ok: raise ValueError("No hay datos de altura satelital (Meta/WRI v2) para esta zona, o no hay conexion a internet.")
    p = preparar(ok, lote_wgs); p["fuente"] = "Meta/WRI CHM v2 (satelital, 1 m)"; p["entero"] = True
    return p
