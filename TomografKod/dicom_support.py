import pydicom
def load(path):
    ds = pydicom.dcmread(path)
    return ds
def edit(df):
    if df==None:
        print("Brak pliku")
        return
    print("hejo")