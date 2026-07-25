"""
FITSViewerWidget — GPU-accelerated FITS image viewer (PyQt5).

FitsGLWidget:      OpenGL renderer; STF/stretch applied on GPU via GLSL.
                   Trilinear mipmap LOD for fast, anti-aliased zoom-out
                   (PixInsight-style).  Texture uploaded once as float32;
                   only uniforms change on stretch update → <1 ms response.
FITSViewerWidget:  Composite = FitsGLWidget + compact control bar
                   (mode selector + STF sliders + pixel info).

STF shadow / highlight values are always in **raw pixel units** (same as the
raw FITS data), matching the convention of compute_stf_params().
The sliders provide a normalized (0-1) UI mapped to the current [data_min,
data_max] range; they track automatically after set_data() + auto_stf().

Usage:
    viewer = FITSViewerWidget(parent)
    viewer.set_data(fits_array)         # float32 numpy, 2-D mono or (H,W,3) RGB
    viewer.auto_stf()                   # compute and apply auto-STF
    viewer.fit_in_view()

Overlay markers:
    from apex.gui.widgets.fits_viewer import OverlayMarker
    viewer.set_overlay_markers([OverlayMarker(col=x, row=y, radius=r)])
"""

from __future__ import annotations

import math
from typing import Optional, List, NamedTuple, Tuple

import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QPushButton, QSlider, QOpenGLWidget, QSizePolicy,
)
from PyQt5.QtGui import (
    QOpenGLShaderProgram, QOpenGLShader, QOpenGLBuffer,
    QOpenGLVertexArrayObject, QSurfaceFormat,
    QColor, QPainter, QPen, QBrush, QWheelEvent, QMouseEvent, QKeyEvent, QFont,
)
from PyQt5.QtCore import Qt, pyqtSignal, QPointF, QRectF, QPoint

try:
    from OpenGL import GL as _GL
    _HAS_PYOPENGL = True
except ImportError:
    _HAS_PYOPENGL = False

from apex.utils.stf import compute_stf_params


# ── GLSL shaders ──────────────────────────────────────────────────────────────

_VERT_SRC = """
#version 330 core
layout(location = 0) in vec2 a_pos;
layout(location = 1) in vec2 a_uv;
out vec2 v_uv;
uniform vec2  u_offset;
uniform float u_zoom;
uniform vec2  u_fit_scale;
void main() {
    v_uv = ((a_uv - 0.5) * u_fit_scale) / u_zoom + 0.5 + u_offset;
    gl_Position = vec4(a_pos, 0.0, 1.0);
}
"""

_FRAG_SRC = """
#version 330 core
in vec2 v_uv;
out vec4 frag_color;

uniform sampler2D u_tex;
uniform int   u_channels;
uniform int   u_mode;       // 0=linear 1=stf 2=log 3=asinh 4=sqrt
uniform float u_shadow;
uniform float u_highlight;
uniform float u_mid;
uniform float u_p_low;
uniform float u_p_high;
uniform float u_strength;
uniform float u_lod;
uniform vec3  u_bg;

float mtf(float x) {
    float m = clamp(u_mid, 0.0001, 0.9999);
    float denom = (2.0*m - 1.0)*x - m;
    if (abs(denom) < 1e-7) denom = sign(denom + 1e-10) * 1e-7;
    return clamp((m - 1.0)*x / denom, 0.0, 1.0);
}

float stretch(float raw) {
    if (u_mode == 1) {
        float norm = clamp((raw - u_shadow) / max(u_highlight - u_shadow, 1e-10), 0.0, 1.0);
        return mtf(norm);
    } else if (u_mode == 2) {
        float norm = clamp((raw - u_p_low) / max(u_p_high - u_p_low, 1e-10), 0.0, 1.0);
        float s = max(u_strength, 1e-6);
        return log(1.0 + norm * s) / log(1.0 + s);
    } else if (u_mode == 3) {
        float norm = clamp((raw - u_p_low) / max(u_p_high - u_p_low, 1e-10), 0.0, 1.0);
        float s = max(u_strength, 1e-6);
        return asinh(norm * s) / asinh(s);
    } else if (u_mode == 4) {
        float norm = clamp((raw - u_p_low) / max(u_p_high - u_p_low, 1e-10), 0.0, 1.0);
        return sqrt(norm);
    } else {
        return clamp((raw - u_p_low) / max(u_p_high - u_p_low, 1e-10), 0.0, 1.0);
    }
}

void main() {
    if (v_uv.x < 0.0 || v_uv.x > 1.0 || v_uv.y < 0.0 || v_uv.y > 1.0) {
        frag_color = vec4(u_bg, 1.0);
        return;
    }
    if (u_channels == 3) {
        vec3 rgb = textureLod(u_tex, v_uv, u_lod).rgb;
        frag_color = vec4(stretch(rgb.r), stretch(rgb.g), stretch(rgb.b), 1.0);
    } else {
        float v = stretch(textureLod(u_tex, v_uv, u_lod).r);
        frag_color = vec4(v, v, v, 1.0);
    }
}
"""


# ── Overlay marker ────────────────────────────────────────────────────────────

class OverlayMarker(NamedTuple):
    col:      float
    row:      float
    radius:   float = 8.0
    color:    QColor = QColor(0, 220, 255, 200)
    rejected: bool = False
    label:    str = ""
    member:   bool = False    # if True, draw small filled center dot
    inner_radius: float = 0.0
    outer_radius: float = 0.0
    secondary_color: QColor = QColor(0, 220, 255, 150)
    line_width: float = 1.4
    target_col: float = math.nan
    target_row: float = math.nan
    target_color: QColor = QColor(255, 0, 255, 180)
    label_offset_col: float = 0.0
    label_offset_row: float = 0.0
    label_font_size: int = 8


# ── FitsGLWidget ──────────────────────────────────────────────────────────────

class FitsGLWidget(QOpenGLWidget):
    """
    OpenGL FITS renderer.  STF/stretch runs on GPU via GLSL uniforms;
    no texture re-upload on stretch change.
    Trilinear mipmap provides PixInsight-style zoom-out averaging.
    """

    mouse_moved    = pyqtSignal(float, float, float)   # img_x, img_y, pixel_val
    mouse_pressed  = pyqtSignal(float, float, int)     # img_x, img_y, Qt.MouseButton
    mouse_released = pyqtSignal(float, float, int)     # img_x, img_y, Qt.MouseButton
    zoom_changed   = pyqtSignal(float)
    view_changed   = pyqtSignal()
    selection_changed  = pyqtSignal(float, float, float, float)  # x0, y0, x1, y1 (img coords)
    selection_finished = pyqtSignal(float, float, float, float)

    def __init__(self, parent=None):
        fmt = QSurfaceFormat()
        fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CoreProfile)
        fmt.setSamples(0)
        QSurfaceFormat.setDefaultFormat(fmt)
        super().__init__(parent)

        self._raw_data: Optional[np.ndarray] = None
        self._channels = 1

        # stretch params (all in raw pixel units)
        self._shadow    = 0.0
        self._highlight = 1.0
        self._mid       = 0.5
        self._p_low     = 0.0
        self._p_high    = 1.0
        self._strength  = 10.0
        self._mode      = 1     # 1 = STF default

        # view state
        self._zoom   = 1.0
        self._offset = QPointF(0.0, 0.0)
        self._pan_start:        Optional[QPointF] = None
        self._pan_offset_start: Optional[QPointF] = None
        self._fit_fill = False

        self._overlay_markers: List[OverlayMarker] = []
        self._background_color = QColor("#1e1e2e")

        # Rectangle selection state
        self._selection_mode: bool = False
        self._sel_dragging: bool = False
        self._sel_start_img: Optional[Tuple[float, float]] = None
        self._sel_rect_img: Optional[Tuple[float, float, float, float]] = None  # x0,y0,x1,y1

        # GL objects (created in initializeGL)
        self._program:     Optional[QOpenGLShaderProgram]    = None
        self._vao:         Optional[QOpenGLVertexArrayObject] = None
        self._vbo:         Optional[QOpenGLBuffer]            = None
        self._tex_id:      Optional[int]                      = None
        self._initialized  = False
        self._pending_data: Optional[np.ndarray]             = None

        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.StrongFocus)
        self.setCursor(Qt.OpenHandCursor)

    # ── public API ───────────────────────────────────────────────────────────

    def set_data(self, data: np.ndarray):
        """Upload float32 numpy array (2-D or H×W×3) to GPU texture."""
        self._raw_data = np.asarray(data, dtype=np.float32)
        if not self._initialized:
            self._pending_data = self._raw_data
            return
        self._upload_texture(self._raw_data)
        self.update()

    def set_stretch_mode(self, mode: str):
        """'stf' | 'linear' | 'log' | 'asinh' | 'sqrt'"""
        self._mode = {"linear": 0, "stf": 1, "log": 2, "asinh": 3, "sqrt": 4}.get(
            mode.lower(), 1
        )
        self.update()

    def set_stf(self, shadow: float, highlight: float, mid: float):
        """Raw pixel values for shadow/highlight; mid is 0-1."""
        self._shadow    = float(shadow)
        self._highlight = float(highlight)
        self._mid       = float(mid)
        self.update()

    def set_linear_range(self, p_low: float, p_high: float):
        """Raw pixel values for linear/log/asinh/sqrt black-white point."""
        self._p_low  = float(p_low)
        self._p_high = float(p_high)
        self.update()

    def set_overlay_markers(self, markers: List[OverlayMarker]):
        self._overlay_markers = list(markers)
        self.update()

    def clear_overlay_markers(self):
        self._overlay_markers = []
        self.update()

    def set_background_color(self, color):
        if not isinstance(color, QColor):
            color = QColor(color)
        if color.isValid():
            self._background_color = color
            self.update()

    def fit_in_view(self, fill: bool = False):
        self._fit_fill = bool(fill)
        self._zoom   = 1.0
        self._offset = QPointF(0.0, 0.0)
        self.update()
        self.zoom_changed.emit(self._zoom)

    def set_selection_mode(self, enabled: bool):
        """Toggle rectangle selection mode (left-drag to draw)."""
        self._selection_mode = bool(enabled)
        if not enabled:
            self._sel_dragging = False
        self.setCursor(Qt.CrossCursor if enabled else Qt.OpenHandCursor)

    def set_selection_rect(self, x0: float, y0: float, x1: float, y1: float):
        """Set selection rectangle in image pixel coords."""
        self._sel_rect_img = (float(x0), float(y0), float(x1), float(y1))
        self.update()

    def clear_selection_rect(self):
        self._sel_rect_img = None
        self._sel_dragging = False
        self.update()

    def get_selection_rect(self) -> Optional[Tuple[float, float, float, float]]:
        return self._sel_rect_img

    @property
    def zoom_level(self) -> float:
        return self._zoom

    # ── OpenGL lifecycle ─────────────────────────────────────────────────────

    def initializeGL(self):
        if not _HAS_PYOPENGL:
            return
        from OpenGL import GL as gl

        self._program = QOpenGLShaderProgram(self)
        self._program.addShaderFromSourceCode(QOpenGLShader.Vertex,   _VERT_SRC)
        self._program.addShaderFromSourceCode(QOpenGLShader.Fragment, _FRAG_SRC)
        self._program.link()

        verts = np.array([
            -1, -1,  0, 1,
             1, -1,  1, 1,
             1,  1,  1, 0,
            -1, -1,  0, 1,
             1,  1,  1, 0,
            -1,  1,  0, 0,
        ], dtype=np.float32)

        self._vao = QOpenGLVertexArrayObject(self)
        self._vao.create()
        self._vao.bind()

        self._vbo = QOpenGLBuffer(QOpenGLBuffer.VertexBuffer)
        self._vbo.create()
        self._vbo.bind()
        self._vbo.allocate(verts.tobytes(), verts.nbytes)

        stride = 4 * 4
        self._program.bind()
        self._program.enableAttributeArray(0)
        self._program.setAttributeBuffer(0, gl.GL_FLOAT, 0,     2, stride)
        self._program.enableAttributeArray(1)
        self._program.setAttributeBuffer(1, gl.GL_FLOAT, 2 * 4, 2, stride)

        self._vao.release()
        self._vbo.release()
        self._program.release()

        gl.glClearColor(0.12, 0.12, 0.18, 1.0)
        self._initialized = True

        if self._pending_data is not None:
            self._upload_texture(self._pending_data)
            self._pending_data = None

    def _upload_texture(self, data: np.ndarray):
        if not _HAS_PYOPENGL:
            return
        from OpenGL import GL as gl

        self.makeCurrent()

        if self._tex_id is not None:
            gl.glDeleteTextures([self._tex_id])

        arr = np.asarray(data, dtype=np.float32)
        if arr.ndim == 2:
            self._channels = 1
            h, w = arr.shape
            gl_internal = gl.GL_R32F
            gl_format   = gl.GL_RED
            tex_data    = np.ascontiguousarray(arr)
        else:
            self._channels = 3
            h, w = arr.shape[:2]
            gl_internal = gl.GL_RGB32F
            gl_format   = gl.GL_RGB
            tex_data    = np.ascontiguousarray(arr[:, :, :3])

        tid = gl.glGenTextures(1)
        gl.glBindTexture(gl.GL_TEXTURE_2D, tid)
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MIN_FILTER,
                           gl.GL_LINEAR_MIPMAP_LINEAR)
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST)
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_S, gl.GL_CLAMP_TO_EDGE)
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_T, gl.GL_CLAMP_TO_EDGE)
        gl.glTexImage2D(
            gl.GL_TEXTURE_2D, 0, gl_internal, w, h, 0,
            gl_format, gl.GL_FLOAT, tex_data,
        )
        gl.glGenerateMipmap(gl.GL_TEXTURE_2D)
        gl.glBindTexture(gl.GL_TEXTURE_2D, 0)
        self._tex_id = int(tid)
        self.doneCurrent()

    def paintGL(self):
        if not _HAS_PYOPENGL or not self._initialized:
            return
        from OpenGL import GL as gl

        bg = self._background_color
        gl.glClearColor(bg.redF(), bg.greenF(), bg.blueF(), 1.0)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)

        if self._tex_id is None:
            return

        self._program.bind()
        self._vao.bind()

        gl.glActiveTexture(gl.GL_TEXTURE0)
        gl.glBindTexture(gl.GL_TEXTURE_2D, self._tex_id)

        self._program.setUniformValue("u_tex", 0)

        lod = max(0.0, -math.log2(max(self._zoom, 1e-6)))
        fit_x, fit_y = self._fit_scale()

        self._program.setUniformValue("u_zoom",      float(self._zoom))
        self._program.setUniformValue("u_offset",    float(self._offset.x()), float(self._offset.y()))
        self._program.setUniformValue("u_fit_scale", float(fit_x), float(fit_y))
        self._program.setUniformValue("u_lod",       float(lod))
        self._program.setUniformValue("u_bg",        bg.redF(), bg.greenF(), bg.blueF())
        self._program.setUniformValue("u_channels",  int(self._channels))
        self._program.setUniformValue("u_mode",      int(self._mode))
        self._program.setUniformValue("u_shadow",    float(self._shadow))
        self._program.setUniformValue("u_highlight", float(self._highlight))
        self._program.setUniformValue("u_mid",       float(self._mid))
        self._program.setUniformValue("u_p_low",     float(self._p_low))
        self._program.setUniformValue("u_p_high",    float(self._p_high))
        self._program.setUniformValue("u_strength",  float(self._strength))

        gl.glDrawArrays(gl.GL_TRIANGLES, 0, 6)
        self._vao.release()
        self._program.release()

    def resizeGL(self, w: int, h: int):
        if _HAS_PYOPENGL:
            from OpenGL import GL as gl
            gl.glViewport(0, 0, w, h)

    # ── zoom / pan ───────────────────────────────────────────────────────────

    def wheelEvent(self, event: QWheelEvent):
        factor = 1.15 if event.angleDelta().y() > 0 else 1.0 / 1.15
        self._zoom = max(0.05, min(64.0, self._zoom * factor))
        self.zoom_changed.emit(self._zoom)
        self.update()

    def _img_coords(self, pos: QPointF) -> Tuple[float, float]:
        if self._raw_data is None:
            return (0.0, 0.0)
        uv_x, uv_y = self._widget_pos_to_uv(pos)
        h, w = self._raw_data.shape[:2]
        return float(uv_x * w), float(uv_y * h)

    def mousePressEvent(self, event: QMouseEvent):
        btn = event.button()
        pos = QPointF(event.pos())

        if self._selection_mode and btn == Qt.LeftButton and self._raw_data is not None:
            ix, iy = self._img_coords(pos)
            self._sel_dragging = True
            self._sel_start_img = (ix, iy)
            self._sel_rect_img = (ix, iy, ix, iy)
            self.update()
        elif (self._selection_mode and btn == Qt.RightButton) or \
             (not self._selection_mode and btn == Qt.LeftButton):
            self._pan_start        = pos
            self._pan_offset_start = QPointF(self._offset)
            self.setCursor(Qt.ClosedHandCursor)

        if self._raw_data is not None:
            ix, iy = self._img_coords(pos)
            self.mouse_pressed.emit(ix, iy, int(btn))
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent):
        pos = QPointF(event.pos())

        if self._sel_dragging and self._sel_start_img is not None:
            ix, iy = self._img_coords(pos)
            x0, y0 = self._sel_start_img
            self._sel_rect_img = (
                min(x0, ix), min(y0, iy),
                max(x0, ix), max(y0, iy),
            )
            self.update()
            self.selection_changed.emit(*self._sel_rect_img)
        elif self._pan_start is not None:
            delta = pos - self._pan_start
            self._offset = QPointF(
                self._pan_offset_start.x() - delta.x() / (self.width()  * self._zoom),
                self._pan_offset_start.y() - delta.y() / (self.height() * self._zoom),
            )
            self.update()

        if self._raw_data is not None:
            uv_x, uv_y = self._widget_pos_to_uv(pos)
            h, w = self._raw_data.shape[:2]
            px, py = int(uv_x * w), int(uv_y * h)
            if 0 <= px < w and 0 <= py < h:
                val = float(np.mean(self._raw_data[py, px])) if self._raw_data.ndim == 3 \
                      else float(self._raw_data[py, px])
                self.mouse_moved.emit(float(px), float(py), val)

    def mouseReleaseEvent(self, event: QMouseEvent):
        btn = event.button()
        pos = QPointF(event.pos())

        if self._sel_dragging and btn == Qt.LeftButton:
            self._sel_dragging = False
            if self._sel_rect_img is not None:
                self.selection_finished.emit(*self._sel_rect_img)
        elif self._pan_start is not None and \
             ((self._selection_mode and btn == Qt.RightButton) or
              (not self._selection_mode and btn == Qt.LeftButton)):
            self._pan_start = None
            self.setCursor(Qt.CrossCursor if self._selection_mode else Qt.OpenHandCursor)
            self.view_changed.emit()

        if self._raw_data is not None:
            ix, iy = self._img_coords(pos)
            self.mouse_released.emit(ix, iy, int(btn))

    def keyPressEvent(self, event: QKeyEvent):
        # Let unhandled keys propagate to parent (step window)
        event.ignore()

    # ── overlay painting (QPainter on top of GL) ─────────────────────────────

    def paintEvent(self, event):
        super().paintEvent(event)
        if self._raw_data is None:
            return
        if not self._overlay_markers and self._sel_rect_img is None:
            return

        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setRenderHint(QPainter.TextAntialiasing)
        img_h, img_w = self._raw_data.shape[:2]
        fit_x, fit_y = self._fit_scale()
        if self._sel_rect_img is not None:
            x0, y0, x1, y1 = self._sel_rect_img
            p0 = self._uv_to_widget_pos(x0 / img_w, y0 / img_h)
            p1 = self._uv_to_widget_pos(x1 / img_w, y1 / img_h)
            rect = QRectF(p0, p1).normalized()
            pen = QPen(QColor(255, 60, 60, 230), 2.0, Qt.DashLine)
            painter.setPen(pen)
            painter.setBrush(Qt.NoBrush)
            painter.drawRect(rect)

        for m in self._overlay_markers:
            uv_x = m.col / img_w
            uv_y = m.row / img_h
            wp   = self._uv_to_widget_pos(uv_x, uv_y)
            wx, wy = wp.x(), wp.y()
            scale = min(
                self.width()  / max(1.0, img_w * fit_x) * self._zoom,
                self.height() / max(1.0, img_h * fit_y) * self._zoom,
            )
            r_px = max(4.0, m.radius * scale) if m.radius > 0 else 0.0

            if math.isfinite(m.target_col) and math.isfinite(m.target_row):
                target_wp = self._uv_to_widget_pos(m.target_col / img_w, m.target_row / img_h)
                painter.setPen(QPen(m.target_color, max(1.0, m.line_width)))
                painter.drawLine(wp, target_wp)

            pen_color = QColor(255, 80, 80, 230) if m.rejected else m.color
            if r_px > 0:
                painter.setPen(QPen(pen_color, max(0.5, m.line_width)))
                painter.setBrush(Qt.NoBrush)
                painter.drawEllipse(QRectF(wx - r_px, wy - r_px, 2 * r_px, 2 * r_px))
                for extra_radius in (m.inner_radius, m.outer_radius):
                    if extra_radius > 0:
                        rr_px = max(2.0, extra_radius * scale)
                        painter.setPen(QPen(m.secondary_color, max(0.5, m.line_width * 0.75)))
                        painter.drawEllipse(QRectF(wx - rr_px, wy - rr_px, 2 * rr_px, 2 * rr_px))

            if m.member:
                painter.setPen(Qt.NoPen)
                painter.setBrush(QBrush(pen_color))
                dot_r = max(1.5, r_px * 0.25 if r_px > 0 else 2.0)
                painter.drawEllipse(QRectF(wx - dot_r, wy - dot_r, 2 * dot_r, 2 * dot_r))

            if m.label:
                label_font = QFont("Arial", max(6, int(m.label_font_size)), QFont.Bold)
                painter.setFont(label_font)
                if m.label_offset_col or m.label_offset_row:
                    tx = wx + m.label_offset_col * scale
                    ty = wy + m.label_offset_row * scale
                else:
                    tx, ty = wx - r_px - 2, wy - r_px - 2
                # Dark outline pass
                painter.setPen(QPen(QColor(0, 0, 0, 180), 2.5))
                for ox, oy in ((-1, -1), (1, -1), (-1, 1), (1, 1)):
                    painter.drawText(QPointF(tx + ox, ty + oy), m.label)
                # Foreground pass
                painter.setPen(QPen(pen_color, 1.0))
                painter.drawText(QPointF(tx, ty), m.label)

        painter.end()

    # ── helpers ──────────────────────────────────────────────────────────────

    def _fit_scale(self):
        if self._raw_data is None or self.width() <= 0 or self.height() <= 0:
            return 1.0, 1.0
        img_h, img_w = self._raw_data.shape[:2]
        if img_h <= 0 or img_w <= 0:
            return 1.0, 1.0
        img_asp  = float(img_w) / float(img_h)
        view_asp = float(self.width()) / float(self.height())
        if self._fit_fill:
            return (1.0, img_asp / view_asp) if view_asp > img_asp \
                   else (view_asp / img_asp, 1.0)
        return (view_asp / img_asp, 1.0) if view_asp > img_asp \
               else (1.0, img_asp / view_asp)

    def _widget_pos_to_uv(self, pos: QPointF):
        fx, fy = self._fit_scale()
        ux = ((pos.x() / max(1, self.width()))  - 0.5) * fx / self._zoom + 0.5 + self._offset.x()
        uy = ((pos.y() / max(1, self.height())) - 0.5) * fy / self._zoom + 0.5 + self._offset.y()
        return ux, uy

    def _uv_to_widget_pos(self, uv_x: float, uv_y: float) -> QPointF:
        fx, fy = self._fit_scale()
        qx = ((uv_x - self._offset.x() - 0.5) * self._zoom / fx) + 0.5
        qy = ((uv_y - self._offset.y() - 0.5) * self._zoom / fy) + 0.5
        return QPointF(qx * self.width(), qy * self.height())


# ── FITSViewerWidget ─────────────────────────────────────────────────────────

_MODES    = ["STF", "Linear", "Log", "Asinh", "Sqrt"]
_MODE_KEYS = ["stf", "linear", "log", "asinh", "sqrt"]

def _ctrl_ss() -> str:
    """Control-bar QSS from the live theme tokens (built at construction so a
    theme switch applies to newly opened viewers)."""
    from apex.gui.theme import Tokens as T
    return (
        f"QWidget#FitsCtrlBar {{ background: {T.SURFACE_ALT};"
        f" border-bottom: 1px solid {T.BORDER_STRONG}; }}"
        f"QWidget#FitsCtrlBar QLabel {{ color: {T.TEXT_SUB}; font-size: 11px; }}"
        # Compact paddings: the themed 6/14px button+input padding clips the
        # fixed-width compact controls in this 30px bar.
        f"QWidget#FitsCtrlBar QPushButton {{ padding: 2px 8px; }}"
        f"QWidget#FitsCtrlBar QComboBox {{ padding: 2px 6px; }}"
    )


class FITSViewerWidget(QWidget):
    """
    Composite FITS viewer: GPU GL canvas + compact stretch control bar.

    STF shadow/highlight parameters are in **raw pixel units** throughout.
    Sliders show a normalized 0-1 view of the current data range.

    Signals
    -------
    mouse_moved(x, y, val)         — image pixel coords + raw value
    mouse_pressed(x, y, qt_button)  — any button press in image coords
    mouse_released(x, y, qt_button) — any button release in image coords
    zoom_changed(zoom)             — current zoom factor
    """

    mouse_moved        = pyqtSignal(float, float, float)
    mouse_pressed      = pyqtSignal(float, float, int)
    mouse_released     = pyqtSignal(float, float, int)
    zoom_changed       = pyqtSignal(float)
    selection_changed  = pyqtSignal(float, float, float, float)
    selection_finished = pyqtSignal(float, float, float, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # A floor, not a target: the Expanding policy above already claims the
        # spare room. 360 was more than a 1280x704 laptop screen could give
        # Step 4 once the plot and table below took their share, which pushed
        # the whole window past the monitor.
        self.setMinimumHeight(240)

        # data range cache (raw pixel units) — updated in set_data()
        self._data_min: float = 0.0
        self._data_max: float = 1.0

        # current STF params in raw pixel units
        self._shadow_raw:    float = 0.0
        self._highlight_raw: float = 1.0
        self._mid_val:       float = 0.5

        # current linear range in raw pixel units
        self._lin_low_raw:  float = 0.0
        self._lin_high_raw: float = 1.0

        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ── control bar ──────────────────────────────────────────────────────
        ctrl = QWidget()
        ctrl.setObjectName("FitsCtrlBar")
        ctrl.setFixedHeight(30)
        ctrl.setStyleSheet(_ctrl_ss())
        bar = QHBoxLayout(ctrl)
        bar.setContentsMargins(6, 2, 6, 2)
        bar.setSpacing(6)

        self._mode_combo = QComboBox()
        self._mode_combo.addItems(_MODES)
        # 80px clipped "Linear"/"STF" under the themed 6/8px input padding.
        self._mode_combo.setFixedWidth(96)
        self._mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        bar.addWidget(self._mode_combo)

        # STF sliders panel (shadow / mid / highlight + Auto)
        self._stf_panel = QWidget()
        sp = QHBoxLayout(self._stf_panel)
        sp.setContentsMargins(0, 0, 0, 0)
        sp.setSpacing(4)

        sp.addWidget(QLabel("Lo"))
        self._sl_shadow = QSlider(Qt.Horizontal)
        self._sl_shadow.setRange(0, 1000)
        self._sl_shadow.setValue(0)
        self._sl_shadow.setFixedWidth(80)
        self._sl_shadow.valueChanged.connect(self._on_stf_slider)
        sp.addWidget(self._sl_shadow)

        sp.addWidget(QLabel("Mid"))
        self._sl_mid = QSlider(Qt.Horizontal)
        self._sl_mid.setRange(1, 999)
        self._sl_mid.setValue(500)
        self._sl_mid.setFixedWidth(80)
        self._sl_mid.valueChanged.connect(self._on_stf_slider)
        sp.addWidget(self._sl_mid)

        sp.addWidget(QLabel("Hi"))
        self._sl_highlight = QSlider(Qt.Horizontal)
        self._sl_highlight.setRange(0, 1000)
        self._sl_highlight.setValue(1000)
        self._sl_highlight.setFixedWidth(80)
        self._sl_highlight.valueChanged.connect(self._on_stf_slider)
        sp.addWidget(self._sl_highlight)

        self._btn_auto = QPushButton("Auto")
        self._btn_auto.setFixedWidth(64)  # 44px clipped the label to "ut"
        self._btn_auto.clicked.connect(self.auto_stf)
        sp.addWidget(self._btn_auto)

        bar.addWidget(self._stf_panel)

        # Linear range sliders panel
        self._lin_panel = QWidget()
        lp = QHBoxLayout(self._lin_panel)
        lp.setContentsMargins(0, 0, 0, 0)
        lp.setSpacing(4)
        lp.addWidget(QLabel("Lo"))
        self._sl_lin_lo = QSlider(Qt.Horizontal)
        self._sl_lin_lo.setRange(0, 500)
        self._sl_lin_lo.setValue(0)
        self._sl_lin_lo.setFixedWidth(80)
        self._sl_lin_lo.valueChanged.connect(self._on_lin_slider)
        lp.addWidget(self._sl_lin_lo)
        lp.addWidget(QLabel("Hi"))
        self._sl_lin_hi = QSlider(Qt.Horizontal)
        self._sl_lin_hi.setRange(500, 1000)
        self._sl_lin_hi.setValue(1000)
        self._sl_lin_hi.setFixedWidth(80)
        self._sl_lin_hi.valueChanged.connect(self._on_lin_slider)
        lp.addWidget(self._sl_lin_hi)
        self._lin_panel.setVisible(False)
        bar.addWidget(self._lin_panel)

        bar.addStretch()

        self._lbl_pixel = QLabel("X:— Y:— Val:—")
        self._lbl_pixel.setFixedWidth(170)
        bar.addWidget(self._lbl_pixel)

        self._lbl_zoom = QLabel("100%")
        self._lbl_zoom.setFixedWidth(50)
        bar.addWidget(self._lbl_zoom)

        root.addWidget(ctrl)

        # ── GL canvas ────────────────────────────────────────────────────────
        self._gl = FitsGLWidget(self)
        self._gl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self._gl.mouse_moved.connect(self._on_mouse_moved)
        self._gl.mouse_pressed.connect(self.mouse_pressed)
        self._gl.mouse_released.connect(self.mouse_released)
        self._gl.zoom_changed.connect(self._on_zoom_changed)
        self._gl.selection_changed.connect(self.selection_changed)
        self._gl.selection_finished.connect(self.selection_finished)
        root.addWidget(self._gl)

    # ── public API ───────────────────────────────────────────────────────────

    @property
    def gl(self) -> "FitsGLWidget":
        return self._gl

    def set_selection_mode(self, enabled: bool):
        self._gl.set_selection_mode(enabled)

    def set_selection_rect(self, x0: float, y0: float, x1: float, y1: float):
        self._gl.set_selection_rect(x0, y0, x1, y1)

    def clear_selection_rect(self):
        self._gl.clear_selection_rect()

    def get_selection_rect(self):
        return self._gl.get_selection_rect()

    def set_data(self, data: np.ndarray):
        """Upload float32 numpy array. Auto-initializes linear range to 0.1/99.9 percentile."""
        arr = np.asarray(data, dtype=np.float32)
        finite = arr[np.isfinite(arr)].ravel()
        if finite.size:
            self._data_min = float(finite.min())
            self._data_max = float(finite.max())
            step = max(1, finite.size // 65536)
            sample = finite[::step][:65536]
            lo = float(np.percentile(sample, 0.1))
            hi = float(np.percentile(sample, 99.9))
            if hi <= lo:
                hi = lo + 1.0
            self._lin_low_raw  = lo
            self._lin_high_raw = hi
            self._gl.set_linear_range(lo, hi)
            self._update_lin_sliders()
        self._gl.set_data(arr)

    def set_data_auto_stf(self, data: np.ndarray):
        """Upload data and apply auto-STF immediately."""
        self.set_data(data)
        self.set_stretch_mode("stf")
        self.auto_stf()

    def fit_in_view(self):
        self._gl.fit_in_view()

    def set_stretch_mode(self, mode: str):
        idx = _MODE_KEYS.index(mode.lower()) if mode.lower() in _MODE_KEYS else 0
        self._mode_combo.setCurrentIndex(idx)

    def set_stf(self, shadow: float, highlight: float, mid: float):
        """Set STF params (raw pixel units for shadow/highlight, 0-1 for mid)."""
        self._shadow_raw    = float(shadow)
        self._highlight_raw = float(highlight)
        self._mid_val       = float(mid)
        self._gl.set_stf(shadow, highlight, mid)
        self._update_stf_sliders()

    def set_shadow_highlight(self, shadow: float, highlight: float):
        """Update shadow/highlight only (keep current mid). Raw pixel units."""
        self.set_stf(shadow, highlight, self._mid_val)

    def set_linear_range(self, low: float, high: float):
        """Set linear/log/asinh/sqrt black-white point in raw pixel units."""
        self._lin_low_raw = float(low)
        self._lin_high_raw = float(high)
        self._gl.set_linear_range(low, high)
        self._update_lin_sliders()

    def auto_stf(self):
        """Compute auto-STF from current data and apply."""
        data = self._gl._raw_data
        if data is None:
            return
        s, h, m = compute_stf_params(data)
        self.set_stf(s, h, m)

    def set_overlay_markers(self, markers):
        self._gl.set_overlay_markers(markers)

    def clear_overlay_markers(self):
        self._gl.clear_overlay_markers()

    def get_stf_params(self) -> Tuple[float, float, float]:
        """Return (shadow_raw, highlight_raw, mid) — raw pixel units."""
        return self._shadow_raw, self._highlight_raw, self._mid_val

    def get_data_range(self) -> Tuple[float, float]:
        """Return (data_min, data_max) in raw pixel units."""
        return self._data_min, self._data_max

    # ── internal ─────────────────────────────────────────────────────────────

    def _on_mode_changed(self, idx: int):
        key = _MODE_KEYS[idx]
        self._gl.set_stretch_mode(key)
        self._stf_panel.setVisible(key == "stf")
        self._lin_panel.setVisible(key in ("linear", "log", "asinh", "sqrt"))

    def _raw_to_norm(self, raw: float) -> float:
        """Map a raw pixel value to 0-1 slider fraction."""
        rng = self._data_max - self._data_min
        if rng <= 0:
            return 0.0
        return max(0.0, min(1.0, (raw - self._data_min) / rng))

    def _norm_to_raw(self, norm: float) -> float:
        """Map a 0-1 slider fraction to raw pixel value."""
        return self._data_min + norm * (self._data_max - self._data_min)

    def _update_stf_sliders(self):
        for sl, val in [
            (self._sl_shadow,    self._raw_to_norm(self._shadow_raw)),
            (self._sl_mid,       self._mid_val),
            (self._sl_highlight, self._raw_to_norm(self._highlight_raw)),
        ]:
            sl.blockSignals(True)
            sl.setValue(int(val * 1000))
            sl.blockSignals(False)

    def _update_lin_sliders(self):
        lo_n = self._raw_to_norm(self._lin_low_raw)
        hi_n = self._raw_to_norm(self._lin_high_raw)
        for sl, val in [(self._sl_lin_lo, lo_n), (self._sl_lin_hi, hi_n)]:
            sl.blockSignals(True)
            sl.setValue(int(val * 1000))
            sl.blockSignals(False)

    def _on_stf_slider(self):
        shadow_raw    = self._norm_to_raw(self._sl_shadow.value()    / 1000.0)
        mid           =                   self._sl_mid.value()       / 1000.0
        highlight_raw = self._norm_to_raw(self._sl_highlight.value() / 1000.0)
        if highlight_raw <= shadow_raw:
            highlight_raw = shadow_raw + max(1.0, (self._data_max - self._data_min) * 0.001)
        self._shadow_raw    = shadow_raw
        self._highlight_raw = highlight_raw
        self._mid_val       = mid
        self._gl.set_stf(shadow_raw, highlight_raw, mid)

    def _on_lin_slider(self):
        lo_raw = self._norm_to_raw(self._sl_lin_lo.value() / 1000.0)
        hi_raw = self._norm_to_raw(self._sl_lin_hi.value() / 1000.0)
        self._lin_low_raw  = lo_raw
        self._lin_high_raw = hi_raw
        self._gl.set_linear_range(lo_raw, hi_raw)

    def _on_mouse_moved(self, x: float, y: float, val: float):
        self._lbl_pixel.setText(f"X:{int(x)} Y:{int(y)} Val:{val:.1f}")
        self.mouse_moved.emit(x, y, val)

    def _on_zoom_changed(self, z: float):
        self._lbl_zoom.setText(f"{z * 100:.0f}%")
        self.zoom_changed.emit(z)
