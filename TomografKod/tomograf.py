from __future__ import annotations
import pydicom
import os
import threading
import dicom_support
import tkinter as tk
from dataclasses import dataclass
from tkinter import filedialog, messagebox, ttk

import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from PIL import Image

from ct_core import ParallelGeometry, backproject_parallel, normalize01, radon_transform_parallel, rmse


matplotlib.use("TkAgg")


def _to_grayscale_np(path: str) -> np.ndarray:
    img = Image.open(path).convert("L")
    arr = np.asarray(img, dtype=np.float32)
    return arr / 255.0


def _safe_float(v: str, default: float) -> float:
    try:
        return float(v)
    except Exception:
        return default


def _safe_int(v: str, default: int) -> int:
    try:
        return int(float(v))
    except Exception:
        return default


@dataclass
class AppState:
    dicom_file: pydicom.dataset.FileDataset|None =None
    image_in: np.ndarray | None = None
    sinogram: np.ndarray | None = None
    image_out: np.ndarray | None = None
    last_path: str | None = None


class TomografApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("Symulator CT (2D) – Radon + rekonstrukcja (Bresenham)")
        self.state = AppState()

        self._build_ui()
        self._render_all()

    def _build_ui(self) -> None:
        self.root.rowconfigure(0, weight=1)
        self.root.columnconfigure(0, weight=1)

        main = ttk.Frame(self.root, padding=10)
        main.grid(row=0, column=0, sticky="nsew")
        main.columnconfigure(0, weight=0)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        controls = ttk.LabelFrame(main, text="Sterowanie", padding=10)
        controls.grid(row=0, column=0, sticky="nsw", padx=(0, 10))

        btn_load = ttk.Button(controls, text="Wczytaj BMP", command=self.on_load)
        
        btn_load.grid(row=0, column=0, sticky="ew")
        
        ttk.Separator(controls).grid(row=1, column=0, sticky="ew", pady=8)

        form = ttk.Frame(controls)
        form.grid(row=2, column=0, sticky="ew")
        form.columnconfigure(1, weight=1)

        ttk.Label(form, text="∆α [deg]").grid(row=0, column=0, sticky="w")
        self.var_delta = tk.StringVar(value="2.0")
        ttk.Entry(form, textvariable=self.var_delta, width=10).grid(row=0, column=1, sticky="ew", padx=(8, 0))

        ttk.Label(form, text="n detektorów").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.var_n = tk.StringVar(value="180")
        ttk.Entry(form, textvariable=self.var_n, width=10).grid(row=1, column=1, sticky="ew", padx=(8, 0), pady=(6, 0))

        ttk.Label(form, text="l [px]").grid(row=2, column=0, sticky="w", pady=(6, 0))
        self.var_l = tk.StringVar(value="256")
        ttk.Entry(form, textvariable=self.var_l, width=10).grid(row=2, column=1, sticky="ew", padx=(8, 0), pady=(6, 0))

        self.var_iter = tk.BooleanVar(value=False)
        ttk.Checkbutton(controls, text="Tryb iteracyjny (suwak postępu)", variable=self.var_iter, command=self._render_all).grid(
            row=3, column=0, sticky="w", pady=(10, 0)
        )

        self.var_filtered = tk.BooleanVar(value=True)
        ttk.Checkbutton(controls, text="Filtrowanie (Ramp / Ram-Lak)", variable=self.var_filtered).grid(
            row=4, column=0, sticky="w", pady=(6, 0)
        )

        ttk.Label(controls, text="Postęp [liczba kątów]").grid(row=5, column=0, sticky="w", pady=(10, 0))
        self.var_progress = tk.IntVar(value=0)
        self.slider = ttk.Scale(controls, from_=0, to=1, orient="horizontal", command=self.on_progress_slider)
        self.slider.grid(row=6, column=0, sticky="ew")

        ttk.Separator(controls).grid(row=7, column=0, sticky="ew", pady=10)

        ttk.Button(controls, text="Generuj sinogram", command=self.on_generate_sinogram).grid(row=8, column=0, sticky="ew")
        ttk.Button(controls, text="Rekonstrukcja", command=self.on_reconstruct).grid(row=9, column=0, sticky="ew", pady=(6, 0))
        ttk.Button(controls, text="RMSE (wejście vs wynik)", command=self.on_rmse).grid(row=10, column=0, sticky="ew", pady=(6, 0))
        dicom_button= ttk.Button(controls, text="Zapisz jako Dicom", command= lambda: dicom_support.DicomEditor(tk.Toplevel(),self.state.image_out,self.state.dicom_file))
        dicom_button.grid(row=11, column=0, sticky="ew")
        self.status = tk.StringVar(value="Gotowe.")
        ttk.Label(controls, textvariable=self.status, wraplength=260).grid(row=12, column=0, sticky="w", pady=(10, 0))

        for i in range(12):
            controls.rowconfigure(i, pad=2)
        controls.columnconfigure(0, weight=1)

        view = ttk.Frame(main)
        view.grid(row=0, column=1, sticky="nsew")
        view.rowconfigure(0, weight=1)
        view.columnconfigure(0, weight=1)

        self.fig = Figure(figsize=(10, 4), dpi=100)
        self.ax_in = self.fig.add_subplot(1, 3, 1)
        self.ax_sino = self.fig.add_subplot(1, 3, 2)
        self.ax_out = self.fig.add_subplot(1, 3, 3)

        self.canvas = FigureCanvasTkAgg(self.fig, master=view)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    def _current_params(self) -> tuple[float, ParallelGeometry]:
        delta = _safe_float(self.var_delta.get(), 2.0)
        n = max(1, _safe_int(self.var_n.get(), 180))
        l = max(1.0, _safe_float(self.var_l.get(), 256.0))
        return delta, ParallelGeometry(n_detectors=n, span_px=l)

    def _max_angles(self) -> int:
        delta, _ = self._current_params()
        return int(np.ceil(180.0 / float(delta)))

    def _sync_slider_range(self) -> None:
        mx = max(1, self._max_angles())
        self.slider.configure(from_=0, to=mx)
        if float(self.slider.get()) > mx:
            self.slider.set(mx)

    def set_status(self, text: str) -> None:
        self.status.set(text)
        self.root.update_idletasks()

    def on_progress_slider(self, _v: str) -> None:
        if self.var_iter.get():
            self._render_all()

    def on_load(self) -> None:
        path = filedialog.askopenfilename(
            title="Wybierz plik",
            filetypes=[("Bitmapy", "*.bmp;*.dib"), ("Wszystkie pliki", "*.*")],
        )
        if not path:
            return
        try:
            try:
                dicomFile=dicom_support.load(path)
                self.state.dicom_file = dicomFile     
                img=dicomFile.pixel_array
            except:
                self.state.dicom_file = None
                img = _to_grayscale_np(path)
        except Exception as e:
            messagebox.showerror("Błąd", f"Nie udało się wczytać obrazu.\n\n{e}")
            return
        
        self.state.image_in = img
        self.state.sinogram = None
        self.state.image_out = None
        self.state.last_path = path
        
        h, w = img.shape
        self.var_l.set(str(int(0.9 * min(w, h))))
        self._sync_slider_range()
        self.slider.set(self._max_angles() if self.var_iter.get() else 0)

        self.set_status(f"Wczytano: {os.path.basename(path)} ({w}×{h}).")
        self._render_all()

    def _run_bg(self, fn, done_text: str) -> None:
        def runner():
            try:
                fn()
                self.root.after(0, lambda: self.set_status(done_text))
                self.root.after(0, self._render_all)
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Błąd", str(e)))
                self.root.after(0, lambda: self.set_status("Błąd obliczeń."))

        threading.Thread(target=runner, daemon=True).start()

    def on_generate_sinogram(self) -> None:
        if self.state.image_in is None:
            messagebox.showinfo("Info", "Najpierw wczytaj obraz BMP.")
            return

        self._sync_slider_range()
        delta, geom = self._current_params()
        progress = int(round(float(self.slider.get()))) if self.var_iter.get() else None

        def work():
            self.root.after(0, lambda: self.set_status("Liczę sinogram..."))
            self.state.sinogram = radon_transform_parallel(self.state.image_in, delta, geom, progress_angles=progress)
            self.state.image_out = None

        self._run_bg(work, "Sinogram gotowy.")

    def on_reconstruct(self) -> None:
        if self.state.image_in is None:
            messagebox.showinfo("Info", "Najpierw wczytaj obraz BMP.")
            return

        self._sync_slider_range()
        delta, geom = self._current_params()
        progress = int(round(float(self.slider.get()))) if self.var_iter.get() else None
        filtered = bool(self.var_filtered.get())

        def work():
            self.root.after(0, lambda: self.set_status("Rekonstrukcja..."))
            if self.state.sinogram is None:
                self.state.sinogram = radon_transform_parallel(self.state.image_in, delta, geom, progress_angles=progress)
            out = backproject_parallel(
                self.state.sinogram,
                delta,
                out_shape=self.state.image_in.shape,
                geom=geom,
                filtered=filtered,
                progress_angles=progress,
            )
            self.state.image_out = normalize01(out)

        self._run_bg(work, "Rekonstrukcja gotowa.")

    def on_rmse(self) -> None:
        if self.state.image_in is None or self.state.image_out is None:
            messagebox.showinfo("Info", "Wczytaj obraz i wykonaj rekonstrukcję.")
            return
        val = rmse(self.state.image_in, self.state.image_out)
        messagebox.showinfo("RMSE", f"RMSE = {val:.6f} (po normalizacji do 0..1)")

    def _render_all(self) -> None:
        self._sync_slider_range()

        self.ax_in.clear()
        self.ax_sino.clear()
        self.ax_out.clear()

        self.ax_in.set_title("Obraz wejściowy")
        if self.state.image_in is not None:
            self.ax_in.imshow(self.state.image_in, cmap="gray", vmin=0, vmax=1)
        self.ax_in.axis("off")

        self.ax_sino.set_title("Sinogram")
        if self.state.sinogram is not None:
            self.ax_sino.imshow(self.state.sinogram, cmap="gray", aspect="auto")
        self.ax_sino.set_xlabel("detektor")
        self.ax_sino.set_ylabel("kąt")

        self.ax_out.set_title("Obraz wyjściowy")
        if self.state.image_out is not None:
            self.ax_out.imshow(self.state.image_out, cmap="gray", vmin=0, vmax=1)
        self.ax_out.axis("off")

        if self.var_iter.get():
            self.ax_sino.text(
                0.5,
                -0.18,
                f"Postęp: {int(round(float(self.slider.get())))} / {self._max_angles()} kątów",
                transform=self.ax_sino.transAxes,
                ha="center",
                va="top",
                fontsize=9,
            )

        self.fig.tight_layout()
        self.canvas.draw_idle()


def main() -> None:
    root = tk.Tk()
    try:
        root.call("tk", "scaling", 1.2)
    except Exception:
        pass
    TomografApp(root)
    root.minsize(1000, 420)
    root.mainloop()


if __name__ == "__main__":
    main()
