import imageio

# Parametry
output_file = "simulation.gif"  # Nazwa wynikowego pliku GIF
frame_files = [f"../out_{i}.png" for i in range(101)]  # Lista plików
frame_rate = 10  # Liczba klatek na sekundę

# Tworzenie GIF-a
with imageio.get_writer(output_file, mode="I", duration=1/frame_rate) as writer:
    for frame_file in frame_files:
        writer.append_data(imageio.imread(frame_file))

print(f"GIF zapisany jako {output_file}")
