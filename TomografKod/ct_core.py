from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


def bresenham_line(x0: int, y0: int, x1: int, y1: int) -> np.ndarray:
    """
    Zwraca tablicę punktów (x, y) na odcinku bez "dziur".
    Implementacja klasycznego algorytmu Bresenhama.
    """
    points: list[tuple[int, int]] = []
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    err = dx - dy

    x, y = x0, y0
    while True:
        points.append((x, y))
        if x == x1 and y == y1:
            break
        e2 = 2 * err
        if e2 > -dy:
            err -= dy
            x += sx
        if e2 < dx:
            err += dx
            y += sy

    return np.asarray(points, dtype=np.int32)


@dataclass(frozen=True)
class ParallelGeometry:
    """
    Model równoległy: dla kąta alpha (rad) mamy zestaw n promieni równoległych,
    o przesunięciach t równomiernie w zakresie [-l/2, l/2] w jednostkach pikseli.
    """

    n_detectors: int
    span_px: float

    def detector_offsets(self) -> np.ndarray:
        if self.n_detectors < 2:
            return np.array([0.0], dtype=np.float32)
        return np.linspace(-self.span_px / 2.0, self.span_px / 2.0, self.n_detectors, dtype=np.float32)


def _clip_line_to_box(p0: np.ndarray, p1: np.ndarray, w: int, h: int) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Przycina odcinek do prostokąta [0,w-1]x[0,h-1] metodą Liang–Barsky.
    Zwraca (q0,q1) w układzie (x,y) lub None gdy brak przecięcia.
    """
    x0, y0 = float(p0[0]), float(p0[1])
    x1, y1 = float(p1[0]), float(p1[1])
    dx = x1 - x0
    dy = y1 - y0

    p = [-dx, dx, -dy, dy]
    q = [x0 - 0.0, (w - 1.0) - x0, y0 - 0.0, (h - 1.0) - y0]

    u1, u2 = 0.0, 1.0
    for pi, qi in zip(p, q):
        if pi == 0.0:
            if qi < 0.0:
                return None
            continue
        t = qi / pi
        if pi < 0.0:
            if t > u2:
                return None
            if t > u1:
                u1 = t
        else:
            if t < u1:
                return None
            if t < u2:
                u2 = t

    q0 = np.array([x0 + u1 * dx, y0 + u1 * dy], dtype=np.float32)
    q1 = np.array([x0 + u2 * dx, y0 + u2 * dy], dtype=np.float32)
    return q0, q1


def _ray_endpoints_parallel(alpha: float, t: float, w: int, h: int) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Wyznacza końce promienia dla modelu równoległego.
    Linia opisana jest równaniem: (x-cx)*cos(alpha) + (y-cy)*sin(alpha) = t
    Kierunek promienia: d = (-sin(alpha), cos(alpha)) (prostopadły do normalnej).
    """
    cx = (w - 1) / 2.0
    cy = (h - 1) / 2.0

    n = np.array([np.cos(alpha), np.sin(alpha)], dtype=np.float32)
    d = np.array([-np.sin(alpha), np.cos(alpha)], dtype=np.float32)
    p = np.array([cx, cy], dtype=np.float32) + n * float(t)

    diag = float(np.hypot(w, h))
    p0 = p - d * (diag * 2.0)
    p1 = p + d * (diag * 2.0)
    clipped = _clip_line_to_box(p0, p1, w, h)
    if clipped is None:
        return None
    return clipped


def radon_transform_parallel(
    image: np.ndarray,
    delta_alpha_deg: float,
    geom: ParallelGeometry,
    progress_angles: int | None = None,
) -> np.ndarray:
    """
    Transformata Radona liczona przez śledzenie promieni (Bresenham).
    Zwraca sinogram o kształcie [n_angles, n_detectors].

    Absorpcja addytywna: jako wartość projekcji bierzemy ŚREDNIĄ jasność na promieniu,
    żeby uniknąć artefaktów związanych z różną "długością dyskretną" odcinków.
    """
    if image.ndim != 2:
        raise ValueError("image must be 2D (grayscale)")
    if delta_alpha_deg <= 0:
        raise ValueError("delta_alpha_deg must be > 0")

    img = image.astype(np.float32, copy=False)
    h, w = img.shape

    n_angles = int(np.ceil(180.0 / float(delta_alpha_deg)))
    if progress_angles is not None:
        n_angles = max(0, min(n_angles, int(progress_angles)))

    offsets = geom.detector_offsets()
    sino = np.zeros((n_angles, offsets.shape[0]), dtype=np.float32)

    for ai in range(n_angles):
        alpha = np.deg2rad(ai * float(delta_alpha_deg))
        for di, t in enumerate(offsets):
            endpoints = _ray_endpoints_parallel(alpha, float(t), w, h)
            if endpoints is None:
                continue
            p0, p1 = endpoints
            x0, y0 = int(round(float(p0[0]))), int(round(float(p0[1])))
            x1, y1 = int(round(float(p1[0]))), int(round(float(p1[1])))

            pts = bresenham_line(x0, y0, x1, y1)
            xs = pts[:, 0]
            ys = pts[:, 1]
            valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
            xs = xs[valid]
            ys = ys[valid]
            if xs.size == 0:
                continue
            sino[ai, di] = float(np.mean(img[ys, xs]))

    return sino


def ramp_filter_fft(projection: np.ndarray) -> np.ndarray:
    """
    Filtr typu ramp (Ram-Lak) zastosowany do 1D projekcji.
    Implementacja przez FFT, bez gotowych filtrów CT.
    """
    p = projection.astype(np.float32, copy=False)
    n = p.shape[0]
    pad = int(2 ** int(np.ceil(np.log2(max(1, 2 * n)))))
    P = np.fft.rfft(p, n=pad)
    freqs = np.fft.rfftfreq(pad).astype(np.float32)
    H = np.abs(freqs)  # ramp
    out = np.fft.irfft(P * H, n=pad)[:n]
    return out.astype(np.float32, copy=False)


def filter_sinogram(sino: np.ndarray) -> np.ndarray:
    """
    Filtruje każdy "view" (wiersz) sinogramu filtrem ramp.
    """
    if sino.ndim != 2:
        raise ValueError("sinogram must be 2D")
    out = np.empty_like(sino, dtype=np.float32)
    for i in range(sino.shape[0]):
        out[i, :] = ramp_filter_fft(sino[i, :])
    return out


def backproject_parallel(
    sino: np.ndarray,
    delta_alpha_deg: float,
    out_shape: tuple[int, int],
    geom: ParallelGeometry,
    filtered: bool = True,
    progress_angles: int | None = None,
) -> np.ndarray:
    """
    Odwrotna transformata przez projekcję wsteczną.
    Implementacja "ray-driven": dla każdego promienia rozprowadza wartość po pikselach Bresenhamem.
    """
    if sino.ndim != 2:
        raise ValueError("sinogram must be 2D")
    if delta_alpha_deg <= 0:
        raise ValueError("delta_alpha_deg must be > 0")

    n_angles_total = sino.shape[0]
    n_angles = n_angles_total
    if progress_angles is not None:
        n_angles = max(0, min(n_angles_total, int(progress_angles)))

    offsets = geom.detector_offsets()
    if offsets.shape[0] != sino.shape[1]:
        raise ValueError("sinogram detector count != geometry n_detectors")

    use = sino[:n_angles, :]
    if filtered:
        use = filter_sinogram(use)

    h, w = out_shape
    out = np.zeros((h, w), dtype=np.float32)

    for ai in range(n_angles):
        alpha = np.deg2rad(ai * float(delta_alpha_deg))
        for di, t in enumerate(offsets):
            val = float(use[ai, di])
            endpoints = _ray_endpoints_parallel(alpha, float(t), w, h)
            if endpoints is None:
                continue
            p0, p1 = endpoints
            x0, y0 = int(round(float(p0[0]))), int(round(float(p0[1])))
            x1, y1 = int(round(float(p1[0]))), int(round(float(p1[1])))

            pts = bresenham_line(x0, y0, x1, y1)
            xs = pts[:, 0]
            ys = pts[:, 1]
            valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
            xs = xs[valid]
            ys = ys[valid]
            if xs.size == 0:
                continue
            out[ys, xs] += val

    if n_angles > 0:
        out /= float(n_angles)

    return out


def normalize01(img: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    a = img.astype(np.float32, copy=False)
    mn = float(np.min(a))
    mx = float(np.max(a))
    if mx - mn < eps:
        return np.zeros_like(a, dtype=np.float32)
    return (a - mn) / (mx - mn)


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    a01 = normalize01(a)
    b01 = normalize01(b)
    return float(np.sqrt(np.mean((a01 - b01) ** 2)))

