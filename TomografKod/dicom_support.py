import pydicom
def load(path):
    ds = pydicom.dcmread(path)
    return ds
def edit(df):
    if df==None:
        print("Brak pliku")
        return
    print("hejo")
import tkinter as tk
from tkinter import filedialog, messagebox
import pydicom
from pydicom.dataset import Dataset, FileDataset
from PIL import Image, ImageTk
import numpy as np
import datetime
from pydicom.uid import generate_uid
class DicomEditor:
    def __init__(self, root,img,ds):
        self.root = root
        self.root.title("DICOM Patient Data Editor")
        self.root.geometry("900x600")

        self.current_ds = None  
        self.img_label = None
        self.out_image = img
        self.current_ds=ds
        self.setup_ui()

    def setup_ui(self):

        form_frame = tk.Frame(self.root, padx=20, pady=20)
        form_frame.pack(side=tk.LEFT, fill=tk.Y)

        tk.Label(form_frame, text="Dane Pacjenta", font=('Arial', 14, 'bold')).pack(pady=10)

        self.fields = {}
        labels = [
            ("Imię i Nazwisko:", "PatientName"),
            ("ID Pacjenta:", "PatientID"),
            ("Data badania (YYYYMMDD):", "ContentDate"),
            ("Komentarz:", "ImageComments")
        ]

        for label_text, attr in labels:
            tk.Label(form_frame, text=label_text).pack(anchor=tk.W)
            entry = tk.Entry(form_frame, width=30)
            entry.pack(pady=5)
            self.fields[attr] = entry

        tk.Button(form_frame, text="Wczytaj DICOM", command=self.load_file, bg="#e1e1e1").pack(fill=tk.X, pady=10)
        tk.Button(form_frame, text="Zapisz Zmiany", command=self.save_file, bg="#4CAF50", fg="white").pack(fill=tk.X)

        self.image_frame = tk.Frame(self.root, bg="black")
        self.image_frame.pack(side=tk.RIGHT, expand=True, fill=tk.BOTH)
        
        self.canvas = tk.Canvas(self.image_frame, bg="black")
        self.canvas.pack(expand=True, fill=tk.BOTH)
        self.display_image()
    def load_file(self):
        path = filedialog.askopenfilename(filetypes=[("DICOM files", "*.dcm"), ("All files", "*.*")])
        if not path:
            return

        try:
            self.current_ds = pydicom.dcmread(path)

            self.update_entries()

            self.display_image()
            
        except Exception as e:
            messagebox.showerror("Błąd", f"Nie udało się odczytać pliku: {e}")

    def update_entries(self):
        # Czyścimy i wpisujemy nowe dane
        for attr, entry in self.fields.items():
            entry.delete(0, tk.END)
            val = getattr(self.current_ds, attr, "")
            entry.insert(0, str(val))

    def display_image(self):
        

        pix = self.out_image
        

        pix = ((pix - np.min(pix)) / (np.max(pix) - np.min(pix)) * 255).astype(np.uint8)

        img = Image.fromarray(pix)
        print("hello")
       

        img.thumbnail((500, 500))
        self.tk_img = ImageTk.PhotoImage(image=img)
        
        self.canvas.create_image(250, 250, image=self.tk_img, anchor=tk.CENTER)

    def save_file(self):

        # 1. Tworzenie metadanych pliku (nagłówek)
        file_meta = Dataset()
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.7' 
        file_meta.MediaStorageSOPInstanceUID = generate_uid()
        file_meta.ImplementationClassUID = generate_uid()

        file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian 

        # 2. Inicjalizacja pustego pliku DICOM
        filename = "nowy_obraz_pacjenta.dcm"

        ds = FileDataset(filename, {}, file_meta=file_meta, preamble=b"\0" * 128)

        # 3. Dodawanie tagów medycznych i identyfikacyjnych
        ds.PatientName = self.fields['PatientName'].get()
        ds.PatientID = self.fields['PatientID'].get()
        ds.ImageComments = self.fields['ImageComments'].get()
        
        
        dt = datetime.datetime.now()
        ds.ContentDate = self.fields['ContentDate'].get() or dt.strftime('%Y%m%d')
        ds.ContentTime = dt.strftime('%H%M%S')
        ds.StudyDate = dt.strftime('%Y%m%d')

        
        ds.StudyInstanceUID = generate_uid()
        ds.SeriesInstanceUID = generate_uid()
        ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
        ds.SOPClassUID = file_meta.MediaStorageSOPClassUID
        ds.Modality = "OT" 

        img_array = self.out_image.copy() 
        
        # Wymiary obrazu
        ds.Rows = img_array.shape[0]
        ds.Columns = img_array.shape[1]
        
        # Konfiguracja sposobu wyświetlania pikseli
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.PixelRepresentation = 0
        ds.BitsAllocated = 8
        ds.BitsStored = 8
        ds.HighBit = 7
        ds.WindowCenter = 127
        ds.WindowWidth = 255
        ds.RescaleIntercept = "0"
        ds.RescaleSlope = "1"

        # 5. Zapis pikseli - Normalizacja i rzutowanie
        if img_array.dtype != np.uint8 or np.max(img_array) <= 1.0:
            if np.max(img_array) != np.min(img_array):
                img_array = ((img_array - np.min(img_array)) / (np.max(img_array) - np.min(img_array)) * 255.0)
            img_array = np.clip(img_array, 0, 255).astype(np.uint8)
        ds.PixelData = img_array.tobytes()



        # 6. Zapis pliku na dysk
        ds.is_little_endian = True
        ds.is_implicit_VR = False 


        save_path = filedialog.asksaveasfilename(defaultextension=".dcm")
        if save_path:
            ds.save_as(save_path)
            messagebox.showinfo("Sukces", "Plik zapisany pomyślnie!")
if __name__ == "__main__":
    root = tk.Tk()
    app = DicomEditor(root)
    root.mainloop()