import subprocess


def run_matlab_script(color, filename, iteration):
    # Generowanie nazwy pliku wynikowego
    result_filename = f"{filename}_{color}_iteration-{iteration}.jpg"

    # Wybieranie funkcji MATLAB-a na podstawie koloru
    if color == "red":
        matlab_function = "bitmap_h_red"
    elif color == "green":
        matlab_function = "bitmap_h_green"
    elif color == "blue":
        matlab_function = "bitmap_h_blue"
    else:
        matlab_function = "bitmap_h"  # Dla ogólnego przypadku lub koloru innego niż RGB

    # Przygotowanie komendy do uruchomienia
    matlab_command = f"{matlab_function}('{filename}', 4, 4, 1, {iteration}, true, '{result_filename}')"

    # Uruchamianie skryptu MATLAB-a
    subprocess.run(
        ["matlab", "-batch", matlab_command],
        check=True,
    )


def main():
    filename = "aurora.jpg"  # Zamień na nazwę pliku, z którym chcesz pracować

    # Uruchamianie skryptów dla kolorów i iteracji
    for color in ["all", "red", "green", "blue"]:
        for iteration in range(1, 11):
            run_matlab_script(color, filename, iteration)


if __name__ == "__main__":
    main()
