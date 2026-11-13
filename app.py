# -*- coding: utf-8 -*-
"""Islas de conservacion y plan de rolado - tablero web (Streamlit).
Ejecutar:  streamlit run app.py     (o doble clic en iniciar_app.bat)"""
import os, io, hashlib
import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from PIL import Image
import nucleo as N

AQUI = os.path.dirname(os.path.abspath(__file__)); EJ = os.path.join(AQUI, "datos_ejemplo")
LOCAL = os.name == "nt" and not os.environ.get("ISLAS_NUBE")      # en un servidor (Linux) se ocultan las opciones que dependen del disco local
COL_CLASE = {1: "#C4A484", 2: "#E3D96A", 3: "#7CBE5A", 4: "#14642C"}
COL_PLAN = {1: "#2E9E4F", 2: "#14642C", 3: "#E39A3B"}
NOM_PLAN = {1: "Conservar (islas)", 2: "Árboles fuera de islas (se esquivan)", 3: "Rolar"}

st.set_page_config(page_title="Islas de conservación y rolado", page_icon="🌳", layout="wide")
st.markdown("<style>.block-container{padding-top:2.2rem} div[data-testid='stImage'] img{border-radius:10px} div[data-testid='stMetricValue']{font-size:1.9rem}</style>", unsafe_allow_html=True)


def fnum(x, d=2):
    return ("{:,.%df}" % d).format(x).replace(",", "X").replace(".", ",").replace("X", ".")


def hexrgb(h):
    return tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))


@st.cache_data(show_spinner=False)
def leer_lotes(nombre, datos):
    return N.leer_poligonos(nombre, datos)


@st.cache_data(show_spinner=False)
def calcular(clave_chm, _chm, lote_wkt, params):
    from shapely import wkt
    prep = N.preparar_meta(wkt.loads(lote_wkt)) if _chm == "META" else N.preparar(_chm, wkt.loads(lote_wkt) if lote_wkt else None)
    msgs = []
    R = N.procesar(prep, dict(params), log=msgs.append)
    R["mensajes"] = msgs
    return R


# ----------------------------------------------------------------- barra lateral
sb = st.sidebar
sb.title("🌳 Islas y rolado")
_hay_ej = os.path.isdir(EJ) and any(f.lower().endswith(".tif") for f in os.listdir(EJ))
modo = sb.radio("Lote", ["Dibujar en el mapa", "Subir archivo"] + (["Ejemplo incluido"] if _hay_ej else []))
chm = None; clave = None; lote = None; nombre_lote = "Lote"; usar_meta = False; dibujar = False
if modo == "Ejemplo incluido":
    ejemplos = sorted(f for f in os.listdir(EJ) if f.lower().endswith(".tif")) if os.path.isdir(EJ) else []
    if not ejemplos:
        sb.warning("No hay datos de ejemplo en la carpeta datos_ejemplo.")
    else:
        e = sb.selectbox("Modelo de altura (CHM)", ejemplos, format_func=lambda s: s[:-4].replace("_", " "))
        chm = os.path.join(EJ, e); clave = "ej:" + e + str(os.path.getmtime(chm))
        num = "".join(ch for ch in e.split("Potrero_")[-1] if ch.isdigit())
        k = os.path.join(EJ, "Potrero_%s.kml" % num)
        if os.path.exists(k):
            pol = leer_lotes(os.path.basename(k), open(k, "rb").read()); nombre_lote, lote = pol[0]
else:
    if modo == "Dibujar en el mapa":
        dibujar = True
        nombre_lote = sb.text_input("Nombre del lote", "Lote dibujado")
        if st.session_state.get("lote_dibujado"):
            from shapely import wkt as _wkt
            lote = _wkt.loads(st.session_state["lote_dibujado"])
            if sb.button("Borrar el lote dibujado"):
                st.session_state.pop("lote_dibujado", None); st.session_state["mapa_n"] = st.session_state.get("mapa_n", 0) + 1; st.rerun()
    else:
        upl = sb.file_uploader("Límite del lote (KML, KMZ o GeoJSON)", type=["kml", "kmz", "geojson", "json"])
        if upl is not None:
            try:
                pol = leer_lotes(upl.name, upl.getvalue())
                i_ = sb.selectbox("Polígono", range(len(pol)), format_func=lambda i: pol[i][0]) if len(pol) > 1 else 0
                nombre_lote, lote = pol[i_]
            except Exception as ex:
                sb.error(str(ex))
    fuente = sb.radio("Altura de la vegetación", ["Satélite (Meta/WRI v2, automático)", "CHM propio (dron)"])
    if fuente.startswith("Satélite"):
        usar_meta = True
        sb.caption("Se lee por internet solo la ventana del lote, a 1 m. Cobertura global; alturas en metros enteros.")
    else:
        up = sb.file_uploader("CHM (GeoTIFF, hasta 200 MB)", type=["tif", "tiff"])
        ruta = sb.text_input("…o ruta del CHM en esta computadora", placeholder=r"H:\carpeta\CHM.tif") if LOCAL else ""
        if up is not None:
            chm = up.getvalue(); clave = "up:" + hashlib.md5(chm).hexdigest()
        elif ruta.strip():
            ruta = ruta.strip().strip('"')
            if os.path.exists(ruta): chm = ruta; clave = "ruta:" + ruta + str(os.path.getmtime(ruta))
            else: sb.error("No existe el archivo indicado.")
        if modo == "Subir archivo" and lote is None: sb.caption("Sin límite: se usa toda el área con dato del CHM.")

sb.subheader("Costo")
moneda = sb.text_input("Moneda", "$")
costo_ha = sb.number_input("Costo del rolado por hectárea", min_value=0.0, value=0.0, step=1000.0, format="%.2f")

sb.subheader("Diseño")
meta = sb.slider("Superficie a conservar en islas (% del lote)", 10, 90, 50, 5)
T = sb.number_input("Altura mínima de árbol (m)", 1.0, 20.0, 4.0, 0.5)
with sb.expander("Parámetros avanzados"):
    dist = st.number_input("Distancia de unión entre árboles (m)", 2.0, 60.0, 15.0, 1.0)
    radio = st.number_input("Radio equivalente mínimo de isla (m)", 0.0, 60.0, 15.0, 1.0)
    rell = st.number_input("Radio de relleno de huecos y entrantes (m)", 0.0, 60.0, 20.0, 1.0)
    peso = st.number_input("Peso de la densidad de árboles", 0.0, 500.0, 60.0, 10.0)
    tol = st.number_input("Simplificación del contorno (m)", 0.0, 30.0, 6.0, 1.0)
    amin = st.number_input("Superficie mínima de isla y de hueco (m²)", 0.0, 50000.0, 5000.0, 500.0)

# ----------------------------------------------------------------- portada (fotos/portada*.jpg)
import base64
_F = os.path.join(AQUI, "fotos")
_port = sorted(os.path.join(_F, f) for f in os.listdir(_F) if f.lower().startswith("portada") and f.lower().endswith((".jpg", ".jpeg", ".png", ".webp"))) if os.path.isdir(_F) else []


@st.cache_data(show_spinner=False)
def portada_b64(ruta, mtime):
    im = Image.open(ruta).convert("RGB"); im.thumbnail((2000, 2000)); buf = io.BytesIO(); im.save(buf, "JPEG", quality=82)
    return base64.b64encode(buf.getvalue()).decode()


if _port:
    st.markdown(
        "<div style=\"border-radius:14px;overflow:hidden;margin-bottom:14px;background:linear-gradient(90deg,rgba(12,40,20,.88) 0%%,rgba(12,40,20,.55) 45%%,rgba(12,40,20,.05) 100%%),"
        "url('data:image/jpeg;base64,%s') center/cover;padding:46px 38px 40px 38px\">"
        "<div style='color:#fff;font-size:2.5rem;font-weight:700;line-height:1.15'>Islas de conservación<br>y plan de rolado</div>"
        "<div style='color:#E6F2E6;font-size:1.05rem;margin-top:10px'>Manejo de bosques con ganadería integrada · conservar el monte, rolar con criterio</div></div>"
        % portada_b64(_port[0], os.path.getmtime(_port[0])), unsafe_allow_html=True)
else:
    st.title("Islas de conservación y plan de rolado")

# ----------------------------------------------------------------- fotos (carpeta fotos/)
FOTOS = os.path.join(AQUI, "fotos")
CATEG = [("monte", "Monte"), ("vacas", "Vacas"), ("rolado", "Rolado"), ("pasturas", "Pasturas")]
EXT_FOTO = (".jpg", ".jpeg", ".png", ".webp")


@st.cache_data(show_spinner=False)
def miniatura(ruta, mtime, ancho=640, alto=400):
    im = Image.open(ruta).convert("RGB"); r = max(ancho / im.width, alto / im.height)
    im = im.resize((max(ancho, int(im.width * r + 0.5)), max(alto, int(im.height * r + 0.5))))
    x = (im.width - ancho) // 2; y = (im.height - alto) // 2
    return im.crop((x, y, x + ancho, y + alto))


def fotos_de(cat):
    if not os.path.isdir(FOTOS): return []
    return sorted(os.path.join(FOTOS, f) for f in os.listdir(FOTOS) if f.lower().startswith(cat) and f.lower().endswith(EXT_FOTO))


cols = st.columns(4)
for col, (cat, titulo) in zip(cols, CATEG):
    fs = fotos_de(cat)
    with col:
        if fs:
            i = st.session_state.get("foto_" + cat, 0) % len(fs)
            st.image(miniatura(fs[i], os.path.getmtime(fs[i])), caption=titulo, width="stretch")
            if len(fs) > 1 and st.button("Otra foto", key="b_" + cat): st.session_state["foto_" + cat] = i + 1; st.rerun()
        else:
            st.markdown("<div style='height:150px;border:1px dashed #BBB;border-radius:10px;background:#FAFAF7;display:flex;align-items:center;justify-content:center;"
                        "color:#888;font-size:0.85rem;text-align:center;padding:8px'>%s<br>sin foto todavía</div>" % titulo, unsafe_allow_html=True)
if LOCAL:
    with st.expander("Agregar fotos (monte, vacas, rolado, pasturas)"):
        st.caption("Las fotos se guardan en la carpeta «fotos» de la aplicación. También puede copiarlas a mano con nombres que empiecen por monte, vacas, rolado o pasturas.")
        ca, cb = st.columns([1, 3])
        cat_sel = ca.selectbox("Categoría", [c for c, _ in CATEG], format_func=lambda c: dict(CATEG)[c])
        nuevas = cb.file_uploader("Fotos (JPG o PNG)", type=["jpg", "jpeg", "png", "webp"], accept_multiple_files=True, key="sube_fotos")
        if nuevas and st.button("Guardar fotos"):
            os.makedirs(FOTOS, exist_ok=True)
            for f in nuevas:
                limpio = "".join(ch if (ch.isalnum() or ch in "._-") else "_" for ch in os.path.basename(f.name))
                with open(os.path.join(FOTOS, "%s_%s" % (cat_sel, limpio)), "wb") as fh: fh.write(f.getvalue())
            st.success("Se guardaron %d foto(s)." % len(nuevas)); st.rerun()

# ----------------------------------------------------------------- mapa para dibujar el lote
if dibujar:
    import folium
    from folium.plugins import Draw
    from streamlit_folium import st_folium
    from shapely.geometry import shape as _shape, mapping as _mapping
    with st.expander("Dibujar el lote en el mapa", expanded=lote is None):
        st.caption("Use la herramienta de polígono (o rectángulo) de la izquierda del mapa, marque los vértices y cierre la figura. Puede buscar el lugar con las coordenadas.")
        cc1, cc2, cc3 = st.columns([1, 1, 2])
        lat0 = cc1.number_input("Latitud", -90.0, 90.0, st.session_state.get("lat0", -31.3100), format="%.5f")
        lon0 = cc2.number_input("Longitud", -180.0, 180.0, st.session_state.get("lon0", -65.5150), format="%.5f")
        st.session_state["lat0"], st.session_state["lon0"] = lat0, lon0
        if lote is not None:
            cen = lote.centroid; centro = [cen.y, cen.x]
        else: centro = [lat0, lon0]
        m = folium.Map(location=centro, zoom_start=15, tiles=None, control_scale=True)
        folium.TileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", attr="Esri, Maxar, Earthstar Geographics", name="Satélite", max_zoom=19).add_to(m)
        folium.TileLayer("OpenStreetMap", name="Calles").add_to(m)
        if lote is not None:
            folium.GeoJson(_mapping(lote), style_function=lambda f: dict(color="#B00020", weight=3, fillOpacity=0.05)).add_to(m)
            b_ = lote.bounds; m.fit_bounds([[b_[1], b_[0]], [b_[3], b_[2]]])
        Draw(export=False, draw_options=dict(polyline=False, circle=False, marker=False, circlemarker=False, polygon=dict(allowIntersection=False, showArea=True), rectangle=True), edit_options=dict(edit=False, remove=True)).add_to(m)
        folium.LayerControl(collapsed=True).add_to(m)
        sal = st_folium(m, height=520, width=None, returned_objects=["last_active_drawing"], key="mapa_%d" % st.session_state.get("mapa_n", 0))
        d_ = (sal or {}).get("last_active_drawing")
        if d_ and d_.get("geometry", {}).get("type") == "Polygon":
            g_ = _shape(d_["geometry"]).buffer(0)
            if not g_.is_empty and g_.wkt != st.session_state.get("lote_dibujado"):
                st.session_state["lote_dibujado"] = g_.wkt; st.rerun()
    if lote is None:
        st.info("Dibuje el lote en el mapa para calcular las islas, la superficie a rolar y el costo."); st.stop()

if usar_meta and lote is not None:
    chm = "META"; clave = "meta:" + hashlib.md5(lote.wkt.encode()).hexdigest()
if chm is None:
    st.info("Falta el lote." if usar_meta else "Elija un ejemplo o cargue un modelo de altura de la vegetación (CHM) en la barra lateral.")
    st.stop()

params = dict(N.PARAMS, META=meta / 100.0, T=float(T), DIST=float(dist), RADIO=float(radio), RELLENO=float(rell), PESO_ARB=float(peso), TOL=float(tol), A_MIN_PIEZA=float(amin))
try:
    with st.spinner("Leyendo la altura satelital del lote y calculando islas…" if chm == "META" else "Calculando islas…"):
        R = calcular(clave, chm, lote.wkt if lote is not None else "", tuple(sorted(params.items())))
except Exception as ex:
    st.error("No se pudo calcular: %s" % ex); st.stop()

res = R["res"]; entero = R["entero"]
costo_total = res["rolar_ha"] * costo_ha
st.caption("%s · %s ha · resolución de cálculo 1 m · %s" % (nombre_lote, fnum(res["lote_ha"]), "alturas en metros enteros (producto satelital)" if entero else "alturas continuas (dron)"))
if entero:
    st.warning("Este modelo de altura tiene valores enteros (típico de productos satelitales como Meta/WRI). Sirve para ubicar árboles, "
               "pero **no distingue suelo, pasto y arbustal**: esas clases se muestran agrupadas como vegetación baja.")
if R["cobertura"] < 0.98:
    st.warning("El CHM cubre solo el %s %% del lote; los resultados corresponden a la parte con dato." % fnum(100 * R["cobertura"], 0))
if abs(res["pct_islas"] - meta) > 2:
    st.info("Las islas quedaron en %s %% (meta %d %%): se descartan islas y huecos menores a la superficie mínima." % (fnum(res["pct_islas"], 1), meta))

# ----------------------------------------------------------------- indicadores
c1, c2, c3, c4 = st.columns(4)
c1.metric("Superficie del lote", fnum(res["lote_ha"]) + " ha")
c2.metric("A conservar", fnum(res["islas_ha"]) + " ha", "%s %% del lote" % fnum(res["pct_islas"], 1), delta_color="off", delta_arrow="off")
c3.metric("A rolar", fnum(res["rolar_ha"]) + " ha", "%s %% del lote" % fnum(res["pct_rolar"], 1), delta_color="off", delta_arrow="off")
c4.metric("Costo total del rolado", "%s %s" % (moneda, fnum(costo_total, 0)) if costo_ha > 0 else "—", "%s %s por ha" % (moneda, fnum(costo_ha, 0)) if costo_ha > 0 else "cargue el costo por ha", delta_color="off", delta_arrow="off")


def reloj(valor, titulo, color, ref=None):
    g = dict(axis=dict(range=[0, 100], ticksuffix=" %"), bar=dict(color=color, thickness=0.75), bgcolor="#F1F1EC", borderwidth=0)
    if ref is not None: g["threshold"] = dict(line=dict(color="#333", width=3), thickness=0.9, value=ref)
    f = go.Figure(go.Indicator(mode="gauge+number", value=valor, number=dict(suffix=" %", valueformat=".1f"), title=dict(text=titulo, font=dict(size=15)), gauge=g))
    f.update_layout(height=230, margin=dict(l=25, r=25, t=55, b=10), separators=",."); return f


g1, g2, g3 = st.columns(3)
g1.plotly_chart(reloj(res["pct_islas"], "Lote en conservación", COL_PLAN[1], meta), width="stretch")
g2.plotly_chart(reloj(res["pct_arboles_protegidos"], "Árboles dentro de islas", COL_PLAN[2]), width="stretch")
g3.plotly_chart(reloj(res["pct_rolar"], "Lote a rolar", COL_PLAN[3]), width="stretch")

t_mapa, t_veg, t_islas, t_rec, t_desc = st.tabs(["Mapa", "Vegetación", "Islas", "Recomendación y costo", "Descargas"])

# ----------------------------------------------------------------- mapa
with t_mapa:
    prod = st.radio("Producto", ["Plan de rolado", "Clases de altura", "Altura de la vegetación"], horizontal=True)
    h, en, c, plan = R["h"], R["en"] & np.isfinite(R["h"]), R["c"], R["plan"]; ny, nx = h.shape
    rgba = np.zeros((ny, nx, 4), np.uint8)
    if prod == "Plan de rolado":
        for k, col in COL_PLAN.items(): rgba[plan == k] = hexrgb(col) + (255,)
        ley = [(NOM_PLAN[k], COL_PLAN[k]) for k in (1, 2, 3)]
    elif prod == "Clases de altura":
        if entero:
            rgba[en & (c < 4)] = hexrgb(COL_CLASE[3]) + (255,); rgba[c == 4] = hexrgb(COL_CLASE[4]) + (255,)
            ley = [("Vegetación baja (≤ %g m)" % T, COL_CLASE[3]), ("Árboles (> %g m)" % T, COL_CLASE[4])]
        else:
            for k, col in COL_CLASE.items(): rgba[c == k] = hexrgb(col) + (255,)
            ley = [("Suelo desnudo (< 0,10 m)", COL_CLASE[1]), ("Pastos (0,10–0,50 m)", COL_CLASE[2]), ("Arbustal (0,50–%g m)" % T, COL_CLASE[3]), ("Árboles (> %g m)" % T, COL_CLASE[4])]
    else:
        v = np.clip(np.nan_to_num(h) / 8.0, 0, 1); a = np.array([247, 244, 233]); b = np.array([0, 68, 27])
        rgb = (a[None, None, :] * (1 - v[..., None]) + b[None, None, :] * v[..., None]).astype(np.uint8)
        rgba[..., :3] = rgb; rgba[..., 3] = np.where(en, 255, 0)
        ley = [("0 m", "#F7F4E9"), ("4 m", "#7BAA7C"), ("8 m o más", "#00441B")]
    paso = max(1, int(np.ceil(max(nx, ny) / 1400)))
    img = Image.fromarray(rgba[::paso, ::paso])
    tr = R["transform"]; x0, y1 = tr.c, tr.f
    fig = go.Figure()
    fig.add_layout_image(dict(source=img, xref="x", yref="y", x=x0, y=y1, sizex=nx, sizey=ny, sizing="stretch", layer="below"))
    xs, ys = R["lote"].exterior.xy
    fig.add_trace(go.Scatter(x=list(xs), y=list(ys), mode="lines", line=dict(color="#B00020", width=2), name="Límite del lote", hoverinfo="skip"))
    for it in R["islas"]:
        g = it["geom"]; xs, ys = g.exterior.xy
        fig.add_trace(go.Scatter(x=list(xs), y=list(ys), mode="lines", line=dict(color="#111", width=2.5), showlegend=False,
                                 text="Isla I%02d<br>%s ha (%s %% del lote)<br>árboles %s ha" % (it["n"], fnum(it["superficie_ha"]), fnum(it["pct_del_lote"], 1), fnum(it["arboles_ha"])), hoverinfo="text"))
        for r_ in g.interiors:
            xs, ys = r_.xy; fig.add_trace(go.Scatter(x=list(xs), y=list(ys), mode="lines", line=dict(color="#111", width=1.5, dash="dot"), showlegend=False, hoverinfo="skip"))
        p = g.representative_point()
        fig.add_annotation(x=p.x, y=p.y, text="<b>I%02d</b><br>%s ha" % (it["n"], fnum(it["superficie_ha"])), showarrow=False, font=dict(size=12, color="#111"), bgcolor="rgba(255,255,255,0.75)")
    for nom, col in ley:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers", marker=dict(size=13, color=col, symbol="square", line=dict(color="#666", width=0.5)), name=nom))
    fig.update_xaxes(range=[x0, x0 + nx], showgrid=False, zeroline=False, title="Este UTM (m)", tickformat=".0f")
    fig.update_yaxes(range=[y1 - ny, y1], showgrid=False, zeroline=False, scaleanchor="x", scaleratio=1, title="Norte UTM (m)", tickformat=".0f")
    fig.update_layout(height=680, margin=dict(l=10, r=10, t=10, b=10), legend=dict(orientation="h", y=1.04, x=0), plot_bgcolor="white", dragmode="pan")
    st.plotly_chart(fig, width="stretch", config=dict(scrollZoom=True))
    st.caption("Contorno negro: islas de conservación. Pase el cursor sobre una isla para ver su superficie. Coordenadas %s." % R["crs"])

# ----------------------------------------------------------------- vegetacion
tab = pd.DataFrame(R["tabla"])
with t_veg:
    if entero:
        baja = tab[tab.clase.isin(["Suelo desnudo (<0,10 m)", "Pastos (0,10-0,50 m)", "Arbustal"])][["lote_ha", "en_islas_ha", "fuera_ha"]].sum()
        arb = tab[tab.clase == "Arboles"].iloc[0]
        filas = [("Vegetación baja (≤ %g m)" % T, baja.lote_ha, baja.en_islas_ha, baja.fuera_ha, COL_CLASE[3]), ("Árboles (> %g m)" % T, arb.lote_ha, arb.en_islas_ha, arb.fuera_ha, COL_CLASE[4])]
    else:
        nom = {"Suelo desnudo (<0,10 m)": ("Suelo desnudo (< 0,10 m)", 1), "Pastos (0,10-0,50 m)": ("Pastos (0,10–0,50 m)", 2), "Arbustal": ("Arbustal (0,50–%g m)" % T, 3), "Arboles": ("Árboles (> %g m)" % T, 4)}
        filas = [(nom[r.clase][0], r.lote_ha, r.en_islas_ha, r.fuera_ha, COL_CLASE[nom[r.clase][1]]) for r in tab.itertuples() if r.clase in nom]
    df = pd.DataFrame(filas, columns=["Clase", "Lote (ha)", "En islas (ha)", "Fuera de islas (ha)", "color"])
    df["% del lote"] = 100 * df["Lote (ha)"] / df["Lote (ha)"].sum()
    df["% de la clase en islas"] = np.where(df["Lote (ha)"] > 0, 100 * df["En islas (ha)"] / df["Lote (ha)"], 0)
    a, b = st.columns([1, 1.3])
    with a:
        st.subheader("Composición del lote")
        f = go.Figure(go.Pie(labels=df["Clase"], values=df["Lote (ha)"], hole=0.55, marker=dict(colors=df["color"], line=dict(color="white", width=2)), sort=False,
                             texttemplate="%{percent:.1%}", hovertemplate="%{label}<br>%{value:.2f} ha<extra></extra>"))
        f.update_layout(height=360, margin=dict(l=10, r=10, t=10, b=10), legend=dict(orientation="h", y=-0.08),
                        annotations=[dict(text="%s ha" % fnum(res["lote_ha"], 1), showarrow=False, font=dict(size=20))])
        st.plotly_chart(f, width="stretch")
    with b:
        st.subheader("Dentro y fuera de las islas")
        f = go.Figure()
        f.add_trace(go.Bar(y=df["Clase"], x=df["En islas (ha)"], orientation="h", name="En islas (se conserva)", marker_color=COL_PLAN[1], text=[fnum(v) for v in df["En islas (ha)"]], textposition="auto"))
        f.add_trace(go.Bar(y=df["Clase"], x=df["Fuera de islas (ha)"], orientation="h", name="Fuera de islas", marker_color=COL_PLAN[3], text=[fnum(v) for v in df["Fuera de islas (ha)"]], textposition="auto"))
        f.update_layout(barmode="stack", height=360, margin=dict(l=10, r=10, t=10, b=10), xaxis_title="hectáreas", legend=dict(orientation="h", y=1.12), yaxis=dict(autorange="reversed"))
        st.plotly_chart(f, width="stretch")
    st.dataframe(df.drop(columns="color").style.format({"Lote (ha)": fnum, "En islas (ha)": fnum, "Fuera de islas (ha)": fnum, "% del lote": lambda v: fnum(v, 1), "% de la clase en islas": lambda v: fnum(v, 1)}), width="stretch", hide_index=True)
    k1, k2, k3 = st.columns(3)
    k1.metric("Altura media", fnum(res["altura_media"]) + " m"); k2.metric("Altura percentil 95", fnum(res["altura_p95"]) + " m"); k3.metric("Altura máxima", fnum(res["altura_max"]) + " m")
    hv = R["h"][R["en"] & np.isfinite(R["h"])]
    hist, bordes = np.histogram(np.clip(hv, 0, 12), bins=24, range=(0, 12))
    f = go.Figure(go.Bar(x=(bordes[:-1] + bordes[1:]) / 2, y=hist / 1e4, marker_color=[COL_CLASE[4] if x > T else COL_CLASE[3] for x in bordes[:-1]], width=0.45, hovertemplate="%{x:.1f} m<br>%{y:.2f} ha<extra></extra>"))
    f.add_vline(x=T, line=dict(color="#333", dash="dash"), annotation_text="árbol > %g m" % T, annotation_position="top right")
    f.update_layout(height=280, margin=dict(l=10, r=10, t=30, b=10), xaxis_title="altura de la vegetación (m; 12 = 12 o más)", yaxis_title="hectáreas", title=dict(text="Distribución de alturas", font=dict(size=15)))
    st.plotly_chart(f, width="stretch")

# ----------------------------------------------------------------- islas
with t_islas:
    di = pd.DataFrame([{"Isla": "I%02d" % it["n"], "Superficie (ha)": it["superficie_ha"], "% del lote": it["pct_del_lote"], "Árboles (ha)": it["arboles_ha"],
                        "Arbustal (ha)": it["arbustal_ha"], "Pastos (ha)": it["pastos_ha"], "Suelo (ha)": it["suelo_ha"], "Vértices": it["vertices"]} for it in R["islas"]])
    if di.empty:
        st.warning("No se generaron islas con estos parámetros.")
    else:
        a, b = st.columns([1, 1.4])
        a.metric("Cantidad de islas", len(di)); a.metric("Isla más grande", fnum(di["Superficie (ha)"].max()) + " ha"); a.metric("Superficie media", fnum(di["Superficie (ha)"].mean()) + " ha")
        f = go.Figure(go.Bar(x=di["Isla"], y=di["Superficie (ha)"], marker_color=COL_PLAN[1], text=[fnum(v) for v in di["Superficie (ha)"]], textposition="outside"))
        f.update_layout(height=320, margin=dict(l=10, r=10, t=30, b=10), yaxis_title="hectáreas", title=dict(text="Superficie por isla", font=dict(size=15)))
        b.plotly_chart(f, width="stretch")
        if entero: di = di.drop(columns=["Pastos (ha)", "Suelo (ha)"]).rename(columns={"Arbustal (ha)": "Vegetación baja (ha)"})
        st.dataframe(di.style.format({c_: fnum for c_ in di.columns if "(ha)" in c_} | {"% del lote": lambda v: fnum(v, 1)}), width="stretch", hide_index=True)
    for m in R["mensajes"]: st.caption(m)

# ----------------------------------------------------------------- recomendacion
with t_rec:
    st.subheader("Recomendación")
    st.markdown(
        "- **Conservar %s ha (%s %% del lote)** en %d isla%s, que contienen el **%s %%** de los árboles de más de %g m.\n"
        "- **Rolar %s ha (%s %% del lote)**: todo lo que queda fuera de las islas, salvo las copas de árboles.\n"
        "- **Esquivar %s ha** de árboles de más de %g m que quedan fuera de las islas: no se rolan."
        % (fnum(res["islas_ha"]), fnum(res["pct_islas"], 1), res["n_islas"], "" if res["n_islas"] == 1 else "s", fnum(res["pct_arboles_protegidos"], 1), T,
           fnum(res["rolar_ha"]), fnum(res["pct_rolar"], 1), fnum(res["arboles_fuera_ha"]), T))
    a, b = st.columns([1.2, 1])
    with a:
        partes = [(NOM_PLAN[1], res["islas_ha"], COL_PLAN[1]), (NOM_PLAN[2], res["arboles_fuera_ha"], COL_PLAN[2]), (NOM_PLAN[3], res["rolar_ha"], COL_PLAN[3])]
        f = go.Figure()
        for nom_, v, col in partes:
            f.add_trace(go.Bar(y=["Lote"], x=[v], orientation="h", name=nom_, marker_color=col, text="%s ha" % fnum(v), textposition="inside", insidetextanchor="middle"))
        f.update_layout(barmode="stack", height=230, margin=dict(l=10, r=10, t=40, b=10), xaxis_title="hectáreas", legend=dict(orientation="h", y=-0.75, traceorder="normal"), title=dict(text="Destino de la superficie del lote", font=dict(size=15)))
        st.plotly_chart(f, width="stretch")
        if not entero:
            rc = res["rolar_por_clase"]
            f = go.Figure(go.Bar(x=["Suelo desnudo", "Pastos", "Arbustal"], y=[rc[1], rc[2], rc[3]], marker_color=[COL_CLASE[1], COL_CLASE[2], COL_CLASE[3]], text=[fnum(rc[k]) for k in (1, 2, 3)], textposition="outside"))
            f.update_layout(height=280, margin=dict(l=10, r=10, t=40, b=10), yaxis_title="hectáreas", title=dict(text="Qué se rola, por tipo de vegetación", font=dict(size=15)))
            st.plotly_chart(f, width="stretch")
    with b:
        st.subheader("Costo del rolado")
        if costo_ha > 0:
            st.metric("Monto total", "%s %s" % (moneda, fnum(costo_total, 0)))
            st.markdown("%s ha × %s %s por ha" % (fnum(res["rolar_ha"]), moneda, fnum(costo_ha, 2)))
            st.metric("Costo por hectárea de lote", "%s %s" % (moneda, fnum(costo_total / res["lote_ha"], 0)))
            esc = pd.DataFrame([{"Variación del costo por ha": ("%+d %%" % v) if v else "sin cambio", "Costo por ha": costo_ha * (1 + v / 100), "Monto total": costo_total * (1 + v / 100)} for v in (-20, -10, 0, 10, 20)])
            st.dataframe(esc.style.format({"Costo por ha": lambda v: moneda + " " + fnum(v, 0), "Monto total": lambda v: moneda + " " + fnum(v, 0)}), width="stretch", hide_index=True)
        else:
            st.info("Cargue el costo del rolado por hectárea en la barra lateral para ver el monto total.")
    st.caption("La recomendación surge de la meta de conservación elegida (%d %% del lote) y de la ubicación de los árboles; es un apoyo para la decisión técnica, no la reemplaza. "
               "El costo considera la superficie efectivamente rolada, sin traslados ni tiempos muertos." % meta)

# ----------------------------------------------------------------- descargas
with t_desc:
    base = nombre_lote.replace(" ", "_")
    resumen = pd.DataFrame([("Lote", nombre_lote), ("Superficie del lote (ha)", round(res["lote_ha"], 2)), ("Islas (cantidad)", res["n_islas"]), ("A conservar (ha)", round(res["islas_ha"], 2)), ("A conservar (% del lote)", round(res["pct_islas"], 1)),
                            ("Árboles fuera de islas (ha)", round(res["arboles_fuera_ha"], 2)), ("A rolar (ha)", round(res["rolar_ha"], 2)), ("A rolar (% del lote)", round(res["pct_rolar"], 1)),
                            ("Árboles dentro de islas (%)", round(res["pct_arboles_protegidos"], 1)), ("Costo por ha (%s)" % moneda, costo_ha), ("Costo total (%s)" % moneda, round(costo_total, 2)),
                            ("Altura mínima de árbol (m)", T), ("Meta de conservación (%)", meta)], columns=["Indicador", "Valor"])
    a, b, c_, d = st.columns(4)
    a.download_button("Islas (KML)", N.exportar_kml(R, nombre_lote), file_name="islas_%s.kml" % base, mime="application/vnd.google-earth.kml+xml")
    b.download_button("Islas (GeoJSON)", N.exportar_geojson(R, nombre_lote), file_name="islas_%s.geojson" % base, mime="application/geo+json")
    c_.download_button("Resumen (CSV)", resumen.to_csv(index=False, sep=";").encode("utf-8-sig"), file_name="resumen_%s.csv" % base, mime="text/csv")
    d.download_button("Vegetación (CSV)", tab.to_csv(index=False, sep=";").encode("utf-8-sig"), file_name="vegetacion_%s.csv" % base, mime="text/csv")
    st.dataframe(resumen.astype(str), width="stretch", hide_index=True)
