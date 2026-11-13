# Islas de conservación y plan de rolado

Tablero web (Python + Streamlit) para definir islas de conservación alrededor de los árboles de un lote, estimar la
superficie a rolar y su costo. El lote se dibuja en el mapa, se sube (KML/KMZ/GeoJSON) o se elige de los ejemplos. La
altura de la vegetación se toma del mapa global Meta/WRI Canopy Height Maps v2 (1 m, CC BY 4.0) o de un CHM propio de dron.

## Uso local (Windows)
1. Python 3.10 o superior. 2. `instalar.bat` (una vez). 3. `iniciar_app.bat`.

## Publicar con un enlace (Streamlit Community Cloud)
1. Crear un repositorio en GitHub y subir el contenido de esta carpeta (app.py, nucleo.py, requirements.txt, .streamlit/, fotos/, datos_ejemplo/).
2. Entrar a share.streamlit.io con la cuenta de GitHub, elegir "Create app", el repositorio y `app.py` como archivo principal.
3. Compartir la dirección que entrega el servicio.

En un servidor se ocultan dos opciones que dependen del disco local: la ruta del CHM y el agregado de fotos desde la app
(las fotos se dejan fijas en la carpeta `fotos/` antes de publicar).

## Método y limitaciones
Ver LEEME.txt. Fuente de altura satelital: Meta y World Resources Institute, Canopy Height Maps v2 (arXiv:2603.06382).
