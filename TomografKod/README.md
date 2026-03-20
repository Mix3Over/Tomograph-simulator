# Symulator tomografu komputerowego (2D) – projekt na 3.0

To jest prosta aplikacja okienkowa w Pythonie, która:
- wczytuje obraz wejściowy **bitmapa (BMP), skala szarości**
- wykonuje **transformatę Radona** (model **równoległy**) → **sinogram**
- wykonuje **odwrotną transformatę** przez **filtrowaną projekcję wsteczną** → obraz wynikowy
- wizualizuje: obraz wejściowy, sinogram, obraz wyjściowy
- pozwala liczyć „na raz” albo iteracyjnie (suwak postępu)

Ważne: obliczenia transformaty są zaimplementowane samodzielnie; do przejścia po pikselach promienia używany jest **algorytm Bresenhama**. Nie ma rotowania obrazu w celu symulacji obrotu układu.

## Instalacja

W katalogu projektu:

```bash
python -m pip install -r requirements.txt
```

## Uruchomienie

```bash
python tomograf.py
```

## Jak przetestować (krok po kroku)

1. Uruchom aplikację.
2. Kliknij **Wczytaj BMP** i wybierz prosty obraz w skali szarości:
   - najlepiej BMP 8-bit (albo RGB – aplikacja i tak przekonwertuje do skali szarości),
   - obraz prostokątny jest dozwolony (aplikacja działa na prostokątnych).
3. Ustaw parametry:
   - **∆α** – krok kątowy w stopniach (np. 2.0 → 180/2=90 projekcji),
   - **n** – liczba detektorów (np. 180),
   - **l** – rozpiętość detektorów w pikselach (np. min(w,h) * 0.9).
4. Kliknij **Generuj sinogram**:
   - w trybie „na raz” sinogram policzy się od razu,
   - w trybie iteracyjnym włącz **Tryb iteracyjny** i ruszaj suwakiem **Postęp** – sinogram będzie liczony do wybranego kąta.
5. Kliknij **Rekonstrukcja**:
   - opcja **Filtrowanie** powinna być włączona (daje dużo lepszy wynik),
   - w trybie iteracyjnym suwak kontroluje liczbę użytych projekcji.
6. (Opcjonalnie) kliknij **RMSE** – policzy błąd między wejściem i wynikiem (wspólny zakres 0..1).

### Zalecane wartości do demo
- ∆α: 2.0
- n: 180
- l: 0.9 * min(w, h)
- filtrowanie: włączone

## Ograniczenia (świadome, pod ocenę 3.0)
- Model: **równoległy** (parallel beam).
- Symulacja pochłaniania: **addytywna** (średnia jasność na promieniu).
- Rekonstrukcja: filtrowana projekcja wsteczna z prostym filtrem typu „ramp” w dziedzinie FFT.

